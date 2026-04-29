#!/usr/bin/env Rscript
# Stage 2: GaGa / MiGaGa2 / Ga / limma fits on MAQC, with per-fit caching so
# we can resume after a 45s bash timeout.
args <- commandArgs(FALSE)
this_file <- sub("^--file=", "", args[grepl("^--file=", args)][1])
src_dir <- dirname(normalizePath(this_file, mustWork = TRUE))
source(file.path(src_dir, "benchmark_lib.R"))
bench <- setup_benchmark()
require_pkgs(c("gaga", "EBarrays", "limma"))

cache_file <- result_file("maqc_loaded.rds")
if (!file.exists(cache_file))
  stop("Run Claude/07a_maqc_cache.R first to cache the GEO matrices.")
loaded <- readRDS(cache_file)
affy_meta <- loaded$affy_meta; affy <- loaded$affy
qpcr_meta <- loaded$qpcr_meta; qpcr <- loaded$qpcr
gpl570 <- loaded$gpl570;       gpl4097 <- loaded$gpl4097

n_genes <- env_int("MAQC_N_GENES", 0)
gaga_method <- Sys.getenv("MAQC_GAGA_METHOD", unset = "EM")

first_site_filter <- function(meta) {
  grepl("^MAQC_AFX_1_[ABCD][1-5]$", meta$title)
}
group_from_title <- function(title) {
  sub("^MAQC_[A-Z]+_[0-9]+_([ABCD]).*$", "\\1", title)
}

keep <- first_site_filter(affy_meta)
affy_meta_1 <- affy_meta[keep, , drop = FALSE]
affy_1 <- affy[, affy_meta_1$gsm, drop = FALSE]
colnames(affy_1) <- affy_meta_1$title
groups_affy <- group_from_title(affy_meta_1$title)

if (n_genes > 0 && n_genes < nrow(affy_1)) {
  log_affy_1 <- log2(pmax(affy_1, .Machine$double.eps))
  vars <- apply(log_affy_1, 1, stats::var)
  keep_genes <- order(vars, decreasing = TRUE)[seq_len(n_genes)]
  affy_fit <- affy_1[keep_genes, , drop = FALSE]
} else {
  affy_fit <- affy_1
}

maqc_patterns <- function() {
  pat <- matrix(c(
    0, 0, 0, 0,   # all equal
    0, 1, 1, 1,   # A vs (C, D, B)
    0, 0, 1, 1,   # (A, C) vs (D, B)
    0, 0, 0, 1,   # (A, C, D) vs B
    0, 1, 2, 3    # all distinct (titration)
  ), ncol = 4, byrow = TRUE)
  colnames(pat) <- c("A", "C", "D", "B"); pat
}
maqc_eb_patterns <- function(groups) {
  vals <- list(rep(1, length(groups)),
               ifelse(groups == "A", 1, 2),
               ifelse(groups %in% c("A", "C"), 1, 2),
               ifelse(groups == "B", 2, 1),
               match(groups, c("A", "C", "D", "B")))
  EBarrays::ebPatterns(vapply(vals, paste, collapse = " ", character(1)))
}
fit_gaga_maqc <- function(x, groups, nclust, method) {
  for (fm in unique(c(method, "quickEM", "EM"))) {
    fit <- tryCatch(
      gaga::fitGG(x, groups, patterns = maqc_patterns(), equalcv = TRUE,
                  nclust = nclust, method = fm, trace = FALSE),
      error = function(e) NULL)
    if (is.null(fit)) next
    fit <- tryCatch(
      gaga::parest(fit, x = x, groups = groups, alpha = 0.05),
      error = function(e) NULL)
    if (!is.null(fit)) { fit$method_used <- fm; return(fit) }
  }
  stop("GaGa fit failed.")
}
fit_ga_maqc <- function(x, groups) {
  for (nm in c("raw", "log2")) {
    xx <- if (nm == "raw") x else log2(pmax(x, .Machine$double.eps))
    fit <- tryCatch(
      EBarrays::emfit(data = xx, family = "GG",
                      hypotheses = maqc_eb_patterns(groups),
                      num.iter = 20, verbose = FALSE),
      error = function(e) NULL)
    if (!is.null(fit)) {
      post <- EBarrays::postprob(fit, xx)$pattern
      return(list(fit = fit, pp = post,
                  method_used = paste0("EBarrays_GG_", nm)))
    }
  }
  stop("Ga fit failed.")
}

# ---------------------------------------------------------------------------
# Per-method caching
# ---------------------------------------------------------------------------
cache_dir <- result_file("maqc_fits"); dir.create(cache_dir, recursive = TRUE,
                                                  showWarnings = FALSE)
fit_or_cache <- function(name, fn) {
  f <- file.path(cache_dir, paste0(name, ".rds"))
  if (file.exists(f)) { cat("[cache] loading", name, "\n"); return(readRDS(f)) }
  cat("[fit ] running", name, "\n")
  el <- system.time(out <- fn())[["elapsed"]]
  out$elapsed <- el
  saveRDS(out, f); out
}

target <- Sys.getenv("MAQC_RUN_STAGE", unset = "all")
gg_fit <- if (target %in% c("all", "GaGa")) fit_or_cache("GaGa",
            function() fit_gaga_maqc(affy_fit, groups_affy, 1, gaga_method))
mg_fit <- if (target %in% c("all", "MiGaGa2")) fit_or_cache("MiGaGa2",
            function() fit_gaga_maqc(affy_fit, groups_affy, 2, gaga_method))
ga_fit <- if (target %in% c("all", "Ga")) tryCatch(fit_or_cache("Ga",
            function() fit_ga_maqc(affy_fit, groups_affy)),
            error = function(e) { warning("Ga skipped: ", conditionMessage(e)); NULL })

if (target != "all") {
  cat("Stage '", target, "' done. Re-run with MAQC_RUN_STAGE=all (or other)\n",
      sep = "")
  quit(status = 0)
}

# ---------------------------------------------------------------------------
# Build curves once all three fits exist
# ---------------------------------------------------------------------------
limma_qpcr_validation <- function(qpcr, groups) {
  groups <- factor(groups); design <- stats::model.matrix(~ groups)
  fit <- limma::eBayes(limma::lmFit(
    log2(pmax(qpcr, .Machine$double.eps)), design))
  tt <- limma::topTable(fit, coef = seq_len(ncol(design))[-1],
                        number = Inf, sort.by = "none")
  p <- tt$P.Value; names(p) <- rownames(tt); p
}
limma_f_score <- function(x, groups) {
  groups <- factor(groups); design <- stats::model.matrix(~ groups)
  fit <- limma::eBayes(limma::lmFit(
    log2(pmax(x, .Machine$double.eps)), design))
  tt <- limma::topTable(fit, coef = seq_len(ncol(design))[-1],
                        number = Inf, sort.by = "none")
  p <- tt$P.Value; names(p) <- rownames(tt); -log10(pmax(p, 1e-300))
}
split_symbols <- function(x) {
  x <- toupper(trimws(x)); x <- strsplit(x, "///", fixed = TRUE)
  lapply(x, function(z) {
    z <- trimws(z); z[nzchar(z) & z != "---" & z != "NA"]
  })
}
make_affy_symbol_map <- function(g) {
  syms <- split_symbols(g[["Gene symbol"]])
  rows <- lapply(seq_along(syms), function(i) {
    if (!length(syms[[i]])) return(NULL)
    data.frame(probe_id = g$ID[i], symbol = syms[[i]],
               stringsAsFactors = FALSE)
  })
  unique(do.call(rbind, rows))
}
score_by_symbol <- function(scores, affy_map) {
  idx <- match(affy_map$probe_id, names(scores))
  map <- affy_map[!is.na(idx), , drop = FALSE]
  map$score <- scores[idx[!is.na(idx)]]
  stats::aggregate(score ~ symbol, data = map, FUN = max, na.rm = TRUE)
}
validation_curve <- function(score, validated) {
  ok <- is.finite(score) & !is.na(validated)
  score <- score[ok]; validated <- as.logical(validated[ok])
  ord <- order(score, decreasing = TRUE); validated <- validated[ord]
  n <- length(validated)
  x <- c(0, cumsum(!validated) / n); y <- c(0, cumsum(validated) / n)
  auc <- sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
  std_auc <- if (sum(!validated) == 0 || sum(validated) == 0) NA_real_
             else auc / ((sum(!validated) / n) * (sum(validated) / n))
  data.frame(rank = seq_along(x) - 1, prop_nonvalidated = x,
             prop_validated = y, auc_scaled = auc, auc_standard = std_auc)
}

qpcr_groups <- group_from_title(qpcr_meta$title)
colnames(qpcr) <- qpcr_meta$title
qpcr_p <- limma_qpcr_validation(qpcr, qpcr_groups)
qpcr_annot <- gpl4097[, c("ID", "ORF", "AssayID"), drop = FALSE]
names(qpcr_annot) <- c("qpcr_id", "symbol", "assay_id")
qpcr_annot$symbol <- toupper(trimws(qpcr_annot$symbol))
qpcr_info <- merge(
  data.frame(qpcr_id = rownames(qpcr), qpcr_p = qpcr_p,
             stringsAsFactors = FALSE),
  qpcr_annot, by = "qpcr_id", all.x = TRUE, sort = FALSE)
qpcr_info$validated <- qpcr_info$qpcr_p < 0.05
affy_map <- make_affy_symbol_map(gpl570)

gaga_score <- 1 - gg_fit$pp[, 1]; names(gaga_score) <- rownames(affy_fit)
migaga_score <- 1 - mg_fit$pp[, 1]; names(migaga_score) <- rownames(affy_fit)
limma_score <- limma_f_score(affy_fit, groups_affy)

method_scores <- list(
  GaGa = score_by_symbol(gaga_score, affy_map),
  MiGaGa2 = score_by_symbol(migaga_score, affy_map),
  limma_BH = score_by_symbol(limma_score, affy_map))
if (!is.null(ga_fit)) {
  ga_score <- 1 - ga_fit$pp[, 1]; names(ga_score) <- rownames(affy_fit)
  method_scores <- c(list(Ga = score_by_symbol(ga_score, affy_map)),
                     method_scores)
}

curve_rows <- list(); summary_rows <- list()
for (m in names(method_scores)) {
  dat <- merge(qpcr_info, method_scores[[m]],
               by = "symbol", all.x = FALSE, sort = FALSE)
  curve <- validation_curve(dat$score, dat$validated)
  curve$method <- m; curve_rows[[m]] <- curve
  fin <- curve[nrow(curve), ]
  summary_rows[[m]] <- data.frame(
    method = m, n_qpcr_assays = nrow(qpcr_info),
    n_mapped_assays = nrow(dat),
    n_validated_mapped = sum(dat$validated),
    auc_scaled = fin$auc_scaled, auc_standard = fin$auc_standard)
}
curve_all <- do.call(rbind, curve_rows)
roc_summary <- do.call(rbind, summary_rows)

paper_counts <- data.frame(
  method = rep(c("GaGa", "MiGaGa2"), each = 5), pattern = rep(0:4, 2),
  paper_count = c(20272, 0, 1429, 3935, 29039,
                  16328, 0, 1323, 3697, 33327))
pattern_rows <- list()
if (!is.null(ga_fit))
  pattern_rows[["Ga"]] <- data.frame(method = "Ga", pattern = 0:4,
    count = as.integer(tabulate(max.col(ga_fit$pp, ties.method = "first"),
                                nbins = 5)))
pattern_counts <- do.call(rbind, c(pattern_rows, list(
  data.frame(method = "GaGa", pattern = 0:4,
             count = as.integer(tabulate(max.col(gg_fit$pp,
                                                ties.method = "first"),
                                         nbins = 5))),
  data.frame(method = "MiGaGa2", pattern = 0:4,
             count = as.integer(tabulate(max.col(mg_fit$pp,
                                                ties.method = "first"),
                                         nbins = 5))))))
pattern_counts <- merge(pattern_counts, paper_counts,
                        by = c("method", "pattern"), all.x = TRUE)
pattern_counts$diff_from_paper <- pattern_counts$count - pattern_counts$paper_count

write.csv(pattern_counts, result_file("maqc_pattern_counts.csv"),
          row.names = FALSE)
write.csv(roc_summary, result_file("maqc_roc_summary.csv"),
          row.names = FALSE)
write.csv(curve_all, result_file("maqc_roc_curves.csv"), row.names = FALSE)

plot_maqc <- function(file, device = c("pdf", "png")) {
  device <- match.arg(device)
  if (device == "pdf") grDevices::pdf(file, width = 6.5, height = 5)
  else grDevices::png(file, width = 1300, height = 1000, res = 200)
  plot(NA,
       xlim = c(0, max(curve_all$prop_nonvalidated, na.rm = TRUE)),
       ylim = c(0, max(curve_all$prop_validated, na.rm = TRUE)),
       xlab = "Proportion nonvalidated qPCR assays",
       ylab = "Proportion validated qPCR assays",
       main = "MAQC validation curves")
  cols <- c(Ga = "#666666", GaGa = "#1B9E77",
            MiGaGa2 = "#D95F02", limma_BH = "#7570B3")
  for (m in names(method_scores)) {
    cc <- curve_all[curve_all$method == m, ]
    lines(cc$prop_nonvalidated, cc$prop_validated, col = cols[[m]], lwd = 2)
  }
  leg <- merge(data.frame(method = names(method_scores)), roc_summary,
               by = "method", sort = FALSE)
  legend("bottomright",
         legend = paste0(leg$method, " AUC=",
                         sprintf("%.4f", leg$auc_scaled)),
         col = cols[leg$method], lwd = 2, bty = "n")
  grDevices::dev.off()
}
plot_maqc(result_file("maqc_roc_curves.pdf"), "pdf")
plot_maqc(result_file("maqc_roc_curves.png"), "png")

cat("\n=== MAQC pattern counts vs Rossell 2009 (Tab 3) ===\n")
print(pattern_counts, row.names = FALSE)
cat("\n=== MAQC ROC AUC summary (Fig 2(b)-style) ===\n")
print(roc_summary, row.names = FALSE)
