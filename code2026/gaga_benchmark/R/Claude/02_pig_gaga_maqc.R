#!/usr/bin/env Rscript
# =============================================================================
# 02_pig_gaga_maqc.R
#
# Apply the PIG-augmented Bayesian GaGa sampler to the first-site Affymetrix
# arm of MAQC and compare the qPCR-validation ROC against:
#   * gaga::fitGG (Stirling-EM, Rossell 2009)
#   * limma BH F-test
#
# The MAQC analysis uses 5 titration patterns (all-equal, A vs CDB, AC vs DB,
# B alone, all-different) over 4 pools A, C, D, B in titration-level order.
# =============================================================================

suppressMessages({
  library(gaga)
  library(EBarrays)
  library(limma)
})
src_dir <- normalizePath(dirname(sub("^--file=", "",
  commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][1])))
source(file.path(src_dir, "pig_gaga_lib.R"))

bench <- normalizePath(file.path(src_dir, "..", ".."))
results_dir <- file.path(src_dir, "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
result_file <- function(x) file.path(results_dir, x)

cache_file <- file.path(bench, "Claude", "results", "maqc_loaded.rds")
if (!file.exists(cache_file))
  stop("Run gaga_benchmark/Claude/07a_maqc_cache.R first to cache GEO data.")
loaded <- readRDS(cache_file)

# Subset to high-variance probes for tractable PIG MCMC
n_genes <- as.integer(Sys.getenv("MAQC_N_GENES", unset = "5000"))
n_iter  <- as.integer(Sys.getenv("PIG_GAGA_ITER", unset = "60"))
n_burn  <- as.integer(Sys.getenv("PIG_GAGA_BURN", unset = "20"))
N_pig   <- as.integer(Sys.getenv("PIG_GAGA_NPIG", unset = "60"))
seed    <- as.integer(Sys.getenv("PIG_GAGA_SEED", unset = "42"))

affy_meta <- loaded$affy_meta; affy <- loaded$affy
qpcr_meta <- loaded$qpcr_meta; qpcr <- loaded$qpcr
gpl570 <- loaded$gpl570; gpl4097 <- loaded$gpl4097

first_site <- grepl("^MAQC_AFX_1_[ABCD][1-5]$", affy_meta$title)
affy_meta_1 <- affy_meta[first_site, , drop = FALSE]
affy_1 <- affy[, affy_meta_1$gsm, drop = FALSE]
colnames(affy_1) <- affy_meta_1$title
groups_affy <- factor(sub("^MAQC_[A-Z]+_[0-9]+_([ABCD]).*$", "\\1",
                          affy_meta_1$title),
                      levels = c("A", "C", "D", "B"))

# Subsample to a mix of high-variance and random probes so the qPCR-mapped
# subset includes both validated and non-validated genes.  All-top-variance
# subsetting tends to keep only validated genes (since high-variance probes
# tend to be the differentially-expressed ones), which kills the ROC.
log2_affy <- log2(pmax(affy_1, 1))
set.seed(seed)
if (n_genes >= nrow(affy_1)) {
  keep <- seq_len(nrow(affy_1))
} else {
  vars <- apply(log2_affy, 1, stats::var)
  n_top <- floor(n_genes / 2)
  n_rand <- n_genes - n_top
  top_idx <- order(vars, decreasing = TRUE)[seq_len(n_top)]
  rand_idx <- sample(setdiff(seq_len(nrow(affy_1)), top_idx), n_rand)
  keep <- sort(c(top_idx, rand_idx))
}
affy_fit <- affy_1[keep, , drop = FALSE]
affy_fit <- pmax(affy_fit, 1)
log2_fit <- log2_affy[keep, , drop = FALSE]

cat(sprintf("MAQC PIG-GaGa: n_probes=%d  n_arrays=%d (5 per pool)\n",
            nrow(affy_fit), ncol(affy_fit)))

# Five patterns, columns A, C, D, B
patterns <- matrix(c(0,0,0,0,
                     0,1,1,1,
                     0,0,1,1,
                     0,0,0,1,
                     0,1,2,3),
                   ncol = 4, byrow = TRUE)
colnames(patterns) <- levels(groups_affy)

# ---------- gaga::fitGG ----------------------------------------------------
cat("\n[ref] Fitting gaga::fitGG (Stirling EM) ...\n")
# gaga expects patterns colnames to match groups levels in factor order
groups_for_gaga <- factor(as.character(groups_affy))
patterns_gg <- patterns
colnames(patterns_gg) <- levels(groups_for_gaga)
gg_fit <- gaga::fitGG(affy_fit, groups_for_gaga, patterns = patterns_gg,
                      equalcv = TRUE, nclust = 1,
                      method = "quickEM", trace = FALSE)
gg_fit <- gaga::parest(gg_fit, x = affy_fit, groups = groups_for_gaga, alpha = 0.05)
gg_pp <- gg_fit$pp
gg_score <- 1 - gg_pp[, 1]
gg_par <- gaga::getpar(gg_fit)
cat(sprintf("[ref] gaga::fitGG  alpha=%.3f  alpha0=%.3f  nu=%.3f\n",
            gg_par$balpha, gg_par$a0, gg_par$nu))

# ---------- limma F ---------------------------------------------------------
cat("\n[ref] limma F-test ...\n")
design <- stats::model.matrix(~ groups_affy)
fit_l <- limma::eBayes(limma::lmFit(log2_fit, design))
tt <- limma::topTable(fit_l, coef = 2:ncol(design), number = Inf,
                      sort.by = "none")
limma_score <- -log10(pmax(tt$P.Value, 1e-300))
names(limma_score) <- rownames(tt)

# ---------- PIG-GaGa --------------------------------------------------------
cat(sprintf("\n[PIG] PIG-GaGa Gibbs (iter=%d  burn=%d  N_pig=%d) ...\n",
            n_iter, n_burn, N_pig))
set.seed(seed)
init <- list(alpha = gg_par$balpha, alpha0 = gg_par$a0,
             lambda0 = gg_par$nu,
             pi = colMeans(gg_pp))
mu_pri <- max(gg_par$balpha + 1.5, 2)
mu0_pri <- max(gg_par$a0 + 1.5, 2)

t0 <- Sys.time()
out <- pig_gaga_gibbs(
  X = affy_fit, group = groups_affy, patterns = patterns,
  pi_prior = rep(1, nrow(patterns)),
  delta = 1, mu = mu_pri, eta_prior = 0,
  delta0 = 1, mu0 = mu0_pri, eta_prior0 = 0,
  init = init, n_iter = n_iter, n_burn = n_burn,
  N_pig = N_pig, verbose = TRUE)
elapsed_pig <- as.numeric(Sys.time() - t0, units = "secs")
pig_score <- 1 - out$pp[, 1]
cat(sprintf("[PIG]  alpha=%.3f  alpha0=%.3f  lambda0=%.3f  elapsed=%.1fs\n",
            mean(tail(out$alpha, n_iter - n_burn)),
            mean(tail(out$alpha0, n_iter - n_burn)),
            mean(tail(out$lambda0, n_iter - n_burn)),
            elapsed_pig))

# ---------- qPCR validation -------------------------------------------------
qpcr_groups <- factor(sub("^MAQC_[A-Z]+_[0-9]+_([ABCD]).*$", "\\1",
                          qpcr_meta$title), levels = c("A", "C", "D", "B"))
colnames(qpcr) <- qpcr_meta$title
fit_q <- limma::eBayes(limma::lmFit(
  log2(pmax(qpcr, .Machine$double.eps)),
  stats::model.matrix(~ qpcr_groups)))
tt_q <- limma::topTable(fit_q, coef = 2:ncol(stats::model.matrix(~ qpcr_groups)),
                        number = Inf, sort.by = "none")
qpcr_p <- tt_q$P.Value; names(qpcr_p) <- rownames(tt_q)

qpcr_annot <- gpl4097[, c("ID", "ORF", "AssayID"), drop = FALSE]
names(qpcr_annot) <- c("qpcr_id", "symbol", "assay_id")
qpcr_annot$symbol <- toupper(trimws(qpcr_annot$symbol))
qpcr_info <- merge(
  data.frame(qpcr_id = rownames(qpcr), qpcr_p = qpcr_p,
             stringsAsFactors = FALSE),
  qpcr_annot, by = "qpcr_id", all.x = TRUE, sort = FALSE)
qpcr_info$validated <- qpcr_info$qpcr_p < 0.05

# Affy probe → gene symbol map
split_symbols <- function(x) {
  x <- toupper(trimws(x)); x <- strsplit(x, "///", fixed = TRUE)
  lapply(x, function(z) {
    z <- trimws(z); z[nzchar(z) & z != "---" & z != "NA"]
  })
}
syms <- split_symbols(gpl570[["Gene symbol"]])
rows <- lapply(seq_along(syms), function(i) {
  if (!length(syms[[i]])) return(NULL)
  data.frame(probe_id = gpl570$ID[i], symbol = syms[[i]],
             stringsAsFactors = FALSE)
})
affy_map <- unique(do.call(rbind, rows))

score_by_symbol <- function(scores) {
  names(scores) <- rownames(affy_fit)
  idx <- match(affy_map$probe_id, names(scores))
  m <- affy_map[!is.na(idx), , drop = FALSE]
  m$score <- scores[idx[!is.na(idx)]]
  stats::aggregate(score ~ symbol, data = m, FUN = max, na.rm = TRUE)
}

methods <- list(
  PIG_GaGa = pig_score,
  gaga_fitGG = gg_score,
  limma_F = limma_score
)

validation_curve <- function(score, validated) {
  ok <- is.finite(score) & !is.na(validated)
  score <- score[ok]; validated <- as.logical(validated[ok])
  ord <- order(score, decreasing = TRUE); validated <- validated[ord]
  n <- length(validated)
  x <- c(0, cumsum(!validated) / n); y <- c(0, cumsum(validated) / n)
  auc <- sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
  std_auc <- if (sum(!validated) == 0 || sum(validated) == 0) NA_real_
             else auc / ((sum(!validated) / n) * (sum(validated) / n))
  data.frame(auc_scaled = auc, auc_standard = std_auc)
}
roc_summary <- do.call(rbind, lapply(names(methods), function(m) {
  s_by_sym <- score_by_symbol(methods[[m]])
  dat <- merge(qpcr_info, s_by_sym, by = "symbol", sort = FALSE)
  cu <- validation_curve(dat$score, dat$validated)
  data.frame(method = m, n_mapped_assays = nrow(dat),
             n_validated = sum(dat$validated), cu)
}))
print(roc_summary, row.names = FALSE)

write.csv(roc_summary, result_file("pig_maqc_roc_summary.csv"),
          row.names = FALSE)
saveRDS(list(out = out, gg_pp = gg_pp, limma_score = limma_score,
             roc = roc_summary),
        result_file("pig_maqc.rds"))
cat("\nDone. Outputs in", results_dir, "\n")
