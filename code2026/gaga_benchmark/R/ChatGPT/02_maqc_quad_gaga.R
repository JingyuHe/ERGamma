args <- commandArgs(FALSE)
this_file <- sub("^--file=", "", args[grepl("^--file=", args)][1])
source(file.path(dirname(normalizePath(this_file, mustWork = TRUE)), "quad_gaga_lib.R"))

bench <- setup_paths()
require_pkgs(c("gaga", "Biobase"))

fit_gaga_fallback <- function(x, groups, patterns, method) {
  methods <- unique(c(method, "quickEM", "EM"))
  errors <- character()
  for (m in methods) {
    fit <- tryCatch({
      fit_gaga_original(x, groups, patterns, nclust = 1, method = m, fdr = 0.05)
    }, error = function(e) {
      errors <<- c(errors, paste0(m, ": ", conditionMessage(e)))
      NULL
    })
    if (!is.null(fit)) {
      fit$method_used <- m
      return(fit)
    }
  }
  stop("GaGa fit failed: ", paste(errors, collapse = " | "))
}

plot_validation_curves <- function(curves, out_pdf, out_png) {
  methods <- unique(curves$method)
  cols <- c(GaGa = "#1B9E77", Quad_GaGa = "#D95F02", limma_F = "#7570B3")
  draw <- function() {
    plot(NA, NA, xlim = c(0, max(curves$prop_nonvalidated, na.rm = TRUE)),
         ylim = c(0, max(curves$prop_validated, na.rm = TRUE)),
         xlab = "Proportion non-validated qPCR assays",
         ylab = "Proportion validated qPCR assays")
    abline(0, 1, col = "grey80", lty = 2)
    for (m in methods) {
      cur <- curves[curves$method == m, ]
      lines(cur$prop_nonvalidated, cur$prop_validated,
            col = cols[[m]], lwd = 2)
    }
    legend("bottomright", legend = methods, col = cols[methods], lwd = 2, bty = "n")
  }
  grDevices::pdf(out_pdf, width = 6, height = 5)
  draw()
  grDevices::dev.off()
  grDevices::png(out_png, width = 1500, height = 1200, res = 220)
  draw()
  grDevices::dev.off()
}

method <- Sys.getenv("QUAD_MAQC_FIT_METHOD", unset = "EM")
max_probes <- env_int("QUAD_MAQC_MAX_PROBES", 0)
fdr <- env_num("QUAD_MAQC_FDR", 0.05)
progress_every <- env_int("QUAD_MAQC_PROGRESS_EVERY", 200)
seed <- env_int("QUAD_MAQC_SEED", 2026)
set.seed(seed)

data_dir <- file.path(bench, "data")
affy_file <- file.path(data_dir, "GSE5350-GPL570_series_matrix.txt.gz")
qpcr_file <- file.path(data_dir, "GSE5350-GPL4097_TaqMan_series_matrix.txt.gz")
gpl570_file <- file.path(data_dir, "GPL570.annot.gz")
gpl4097_file <- file.path(data_dir, "GPL4097_family.soft.gz")

cat("Reading MAQC matrices and annotations\n")
affy_meta <- read_series_metadata(affy_file)
affy <- read_series_matrix(affy_file)
qpcr_meta <- read_series_metadata(qpcr_file)
qpcr <- read_series_matrix(qpcr_file)
gpl570 <- read_platform_table(gpl570_file, gz = TRUE)
gpl4097 <- read_platform_table(gpl4097_file, gz = TRUE)

keep <- first_site_filter(affy_meta)
affy_meta_1 <- affy_meta[keep, , drop = FALSE]
affy_1 <- affy[, affy_meta_1$gsm, drop = FALSE]
colnames(affy_1) <- affy_meta_1$title
groups_affy <- factor(group_from_title(affy_meta_1$title), levels = c("A", "C", "D", "B"))

qpcr_groups <- factor(group_from_title(qpcr_meta$title), levels = c("A", "C", "D", "B"))
colnames(qpcr) <- qpcr_meta$title
qpcr_p <- limma_qpcr_validation(qpcr, qpcr_groups)
qpcr_annot <- gpl4097[, c("ID", "ORF", "AssayID"), drop = FALSE]
names(qpcr_annot) <- c("qpcr_id", "symbol", "assay_id")
qpcr_annot$symbol <- toupper(trimws(qpcr_annot$symbol))
qpcr_info <- merge(
  data.frame(qpcr_id = rownames(qpcr), qpcr_p = qpcr_p, stringsAsFactors = FALSE),
  qpcr_annot,
  by = "qpcr_id",
  all.x = TRUE,
  sort = FALSE
)
qpcr_info$validated <- qpcr_info$qpcr_p < 0.05
qpcr_symbols <- unique(qpcr_info$symbol[nzchar(qpcr_info$symbol)])

affy_map <- make_affy_symbol_map(gpl570)
mapped_probes <- unique(affy_map$probe_id[affy_map$symbol %in% qpcr_symbols])
mapped_probes <- intersect(mapped_probes, rownames(affy_1))
affy_fit <- affy_1[mapped_probes, , drop = FALSE]

if (max_probes > 0 && max_probes < nrow(affy_fit)) {
  vars <- apply(log2(pmax(affy_fit, .Machine$double.eps)), 1, stats::var)
  keep_probes <- order(vars, decreasing = TRUE)[seq_len(max_probes)]
  affy_fit <- affy_fit[keep_probes, , drop = FALSE]
}

affy_map_fit <- affy_map[affy_map$probe_id %in% rownames(affy_fit), , drop = FALSE]
patterns <- maqc_patterns()

cat("MAQC Quad-GaGa benchmark\n")
cat("  probes=", nrow(affy_fit), " arrays=", ncol(affy_fit),
    " mapped_qpcr_symbols=", length(qpcr_symbols), "\n", sep = "")
cat("  qPCR validation p-values from ",
    if (requireNamespace("limma", quietly = TRUE)) "limma" else "base F test",
    "\n", sep = "")

elapsed_gaga <- system.time({
  gg_fit <- fit_gaga_fallback(affy_fit, groups_affy, patterns, method = method)
})[["elapsed"]]

elapsed_quad <- system.time({
  quad <- quad_gaga_pp(affy_fit, groups_affy, patterns, gg_fit,
                       progress_every = progress_every)
})[["elapsed"]]

gaga_score <- 1 - gg_fit$pp[, 1]
quad_score <- 1 - quad$pp[, 1]
limma_p <- limma_f_pvalue(affy_fit, groups_affy)
limma_score <- -log10(pmax(limma_p, 1e-300))
names(gaga_score) <- rownames(affy_fit)
names(quad_score) <- rownames(affy_fit)

quad_calls <- find_quad_calls(quad$pp, fdr = fdr)
gaga_calls <- as.logical(gaga::findgenes(gg_fit, affy_fit, as.character(groups_affy),
                                         fdrmax = fdr, parametric = TRUE)$d)

probe_scores <- data.frame(
  probe_id = rownames(affy_fit),
  gaga_score = gaga_score,
  quad_score = quad_score,
  limma_score = limma_score,
  gaga_pp_null = gg_fit$pp[, 1],
  quad_pp_null = quad$pp[, 1],
  gaga_call = gaga_calls,
  quad_call = quad_calls,
  stringsAsFactors = FALSE
)

pattern_counts <- rbind(
  data.frame(method = "GaGa", pattern = 0:4,
             count = as.integer(tabulate(max.col(gg_fit$pp, ties.method = "first"), nbins = 5))),
  data.frame(method = "Quad_GaGa", pattern = 0:4,
             count = as.integer(tabulate(max.col(quad$pp, ties.method = "first"), nbins = 5)))
)

method_scores <- list(
  GaGa = score_by_symbol(gaga_score, affy_map_fit),
  Quad_GaGa = score_by_symbol(quad_score, affy_map_fit),
  limma_F = score_by_symbol(limma_score, affy_map_fit)
)

curve_rows <- list()
score_rows <- list()
summary_rows <- list()
for (m in names(method_scores)) {
  dat <- merge(qpcr_info, method_scores[[m]], by = "symbol", all.x = FALSE, sort = FALSE)
  curve <- validation_curve(dat$score, dat$validated)
  curve$method <- m
  curve_rows[[m]] <- curve
  dat$method <- m
  score_rows[[m]] <- dat
  final <- curve[nrow(curve), ]
  summary_rows[[m]] <- data.frame(
    method = m,
    n_fit_probes = nrow(affy_fit),
    n_qpcr_assays = nrow(qpcr_info),
    n_mapped_assays = nrow(dat),
    n_validated_mapped = sum(dat$validated),
    auc_scaled = final$auc_scaled,
    auc_standard = final$auc_standard,
    spearman_probe_score_vs_gaga = if (m == "GaGa") 1 else
      stats::cor(probe_scores$gaga_score,
                 if (m == "Quad_GaGa") probe_scores$quad_score else probe_scores$limma_score,
                 method = "spearman", use = "pairwise.complete.obs"),
    n_de_at_fdr = if (m == "GaGa") sum(gaga_calls) else if (m == "Quad_GaGa") sum(quad_calls) else NA_integer_,
    stringsAsFactors = FALSE
  )
}

curves <- do.call(rbind, curve_rows)
validation_scores <- do.call(rbind, score_rows)
summary <- do.call(rbind, summary_rows)
summary$gaga_fit_method <- gg_fit$method_used
summary$elapsed_gaga_seconds <- elapsed_gaga
summary$elapsed_quad_seconds <- elapsed_quad
summary$qPCR_validation_engine <- if (requireNamespace("limma", quietly = TRUE)) "limma" else "base_F"

write.csv(summary, chatgpt_result_file("maqc_quad_gaga_validation_summary.csv"),
          row.names = FALSE)
write.csv(curves, chatgpt_result_file("maqc_quad_gaga_validation_curve.csv"),
          row.names = FALSE)
write.csv(validation_scores, chatgpt_result_file("maqc_quad_gaga_validation_scores.csv"),
          row.names = FALSE)
write.csv(probe_scores, chatgpt_result_file("maqc_quad_gaga_probe_scores.csv"),
          row.names = FALSE)
write.csv(pattern_counts, chatgpt_result_file("maqc_quad_gaga_pattern_counts.csv"),
          row.names = FALSE)
saveRDS(list(config = list(max_probes = max_probes, fdr = fdr,
                           fit_method = method, seed = seed),
             summary = summary, curves = curves, validation_scores = validation_scores,
             probe_scores = probe_scores, pattern_counts = pattern_counts,
             fit = gg_fit, quad = quad),
        chatgpt_result_file("maqc_quad_gaga.rds"))

plot_validation_curves(
  curves,
  out_pdf = chatgpt_result_file("maqc_quad_gaga_validation_curve.pdf"),
  out_png = chatgpt_result_file("maqc_quad_gaga_validation_curve.png")
)

print(summary)
print(pattern_counts)
