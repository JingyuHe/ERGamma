#!/usr/bin/env Rscript
# =============================================================================
# 03_pig_pergene_armstrong.R
#
# Per-gene PIG-GaGa for Armstrong using IDENTICAL hyperparameters and model
# specification as gaga::fitGG (per-gene α_i with Ga(b_α, b_α/ν_α) prior).
# The only difference is that we replace gaga's Stirling-based marginal-
# likelihood approximation with PIG-EM + Laplace (PIG paper §3.3).
# =============================================================================

suppressMessages({
  library(gaga)
  library(limma)
})
src_dir <- normalizePath(dirname(sub("^--file=", "",
  commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][1])))
source(file.path(src_dir, "pig_gaga_pergene.R"))

bench <- normalizePath(file.path(src_dir, "..", ".."))
results_dir <- file.path(src_dir, "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
result_file <- function(x) file.path(results_dir, x)

data_choice <- Sys.getenv("ARMSTRONG_DATA", unset = "schliep_filtered")
em_iter <- as.integer(Sys.getenv("PIG_EM_ITER", unset = "30"))

# ---------- data loaders ----------------------------------------------------
load_armstrong_schliep <- function(data_file) {
  raw <- read.delim(data_file, check.names = FALSE, stringsAsFactors = FALSE)
  groups <- as.character(unlist(raw[1, -1]))
  x <- as.matrix(data.frame(lapply(raw[-1, -1, drop = FALSE], as.numeric),
                            check.names = FALSE))
  rownames(x) <- raw[-1, 1]; list(x = x, groups = groups)
}
load_armstrong_orange <- function(data_file) {
  raw <- read.delim(data_file, check.names = FALSE, stringsAsFactors = FALSE)
  groups <- raw$class[-c(1, 2)]
  xdf <- raw[-c(1, 2), -1, drop = FALSE]
  x <- t(apply(xdf, 2, as.numeric))
  rownames(x) <- colnames(raw)[-1]; list(x = x, groups = groups)
}

if (data_choice == "schliep_filtered") {
  dat <- load_armstrong_schliep(file.path(bench, "data",
                                          "armstrong-2002-v2_database.txt"))
  transform_choice <- Sys.getenv("ARMSTRONG_TRANSFORM", unset = "log")
} else {
  dat <- load_armstrong_orange(file.path(bench, "data",
                                         "MLL_orange_12533.tab"))
  transform_choice <- Sys.getenv("ARMSTRONG_TRANSFORM", unset = "log2_floor20")
}
keep <- dat$groups %in% c("ALL", "MLL")
mll_idx <- which(dat$groups == "MLL"); keep[tail(mll_idx, 2)] <- FALSE
X_raw <- dat$x[, keep, drop = FALSE]
groups <- factor(dat$groups[keep], levels = c("ALL", "MLL"))

log_t <- function(X, t) {
  if (t == "log")           log(pmax(X, 1))
  else if (t == "log2_floor20") log2(pmax(X, 20))
  else if (t == "log2_clip")    log2(pmax(X, 1) + 1)
  else stop("unknown transform: ", t)
}
# Both gaga and PIG-pergene fit the same Gamma model on the SAME log-scale
# matrix.  This is what gaga::fitGG expects and what every Rossell-style
# benchmark uses.  Passing a different scale to PIG would make the comparison
# meaningless (the user pointed this out).
X_log <- log_t(X_raw, transform_choice)
X     <- X_log

cat(sprintf("Armstrong per-gene PIG-GaGa: data=%s transform=%s genes=%d samples=%d\n",
            data_choice, transform_choice, nrow(X), ncol(X)))
cat(sprintf("groups: ALL=%d MLL=%d\n",
            sum(groups == "ALL"), sum(groups == "MLL")))

patterns <- matrix(c(0, 0, 0, 1), 2, 2)
colnames(patterns) <- levels(groups)

# ---------- gaga::fitGG to get hyperparameters AND reference DE calls -------
cat("\n[ref] gaga::fitGG (Stirling EM) ...\n")
gg_fit <- gaga::fitGG(X_log, groups, patterns = patterns, equalcv = TRUE,
                      nclust = 1, method = "quickEM", trace = FALSE)
gg_fit <- gaga::parest(gg_fit, x = X_log, groups = groups, alpha = 0.05)
gg_par <- gaga::getpar(gg_fit)
gg_genes <- gaga::findgenes(gg_fit, X_log, groups, fdrmax = 0.05,
                            parametric = TRUE)
gg_calls <- as.logical(gg_genes$d != 0)
gg_score <- 1 - gg_fit$pp[, 1]
cat(sprintf("[ref] hyperparams: alpha_0=%.3f  lambda_0=%.3f  ",
            gg_par$a0, gg_par$nu),
    sprintf("b_alpha=%.3f  nu_alpha=%.3f  pi=(%.3f,%.3f)\n",
            gg_par$balpha, gg_par$nualpha,
            gg_par$probpat[1], gg_par$probpat[2]))
cat(sprintf("[ref] gaga #DE = %d\n", sum(gg_calls)))

# ---------- limma_BH ------------------------------------------------------
fit_l <- limma::eBayes(limma::lmFit(X_log, stats::model.matrix(~ groups)))
limma_p <- fit_l$p.value[, 2]
limma_calls <- p.adjust(limma_p, "BH") <= 0.05
limma_score <- -log10(pmax(limma_p, 1e-300))
cat(sprintf("[ref] limma_BH #DE = %d\n", sum(limma_calls)))

# ---------- per-gene PIG-EM + Laplace -------------------------------------
hyper <- list(alpha0 = gg_par$a0, nu = gg_par$nu,
              b_alpha = gg_par$balpha, nu_alpha = gg_par$nualpha)
pi_v <- as.numeric(gg_par$probpat)

cat(sprintf("\n[PIG] per-gene PIG-EM + Laplace (em_iter=%d) ...\n", em_iter))
t0 <- Sys.time()
pig <- pig_gaga_pergene(X = X, group = groups, patterns = patterns,
                        hyper = hyper, pi_v = pi_v,
                        em_iter = em_iter, verbose = TRUE)
elapsed_pig <- as.numeric(Sys.time() - t0, units = "secs")
pig_score <- 1 - pig$pp[, 1]

# Same FDR rule that gaga::findgenes uses (Newton & Bonferroni-style):
fdr_calls <- function(pp, fdr = 0.05) {
  null_prob <- pp[, 1]
  ord <- order(null_prob)            # most likely DE first
  cum_fdr <- cumsum(null_prob[ord]) / seq_along(ord)
  cutoff_rank <- max(c(0, which(cum_fdr <= fdr)))
  calls <- rep(FALSE, nrow(pp))
  if (cutoff_rank > 0) calls[ord[seq_len(cutoff_rank)]] <- TRUE
  calls
}
pig_calls <- fdr_calls(pig$pp, 0.05)
cat(sprintf("[PIG]  PIG-pergene #DE=%d   elapsed=%.1fs\n",
            sum(pig_calls), elapsed_pig))

# ---------- agreement -------------------------------------------------------
agree <- function(a, b, na, nb) {
  data.frame(method_pair = paste0(na, " vs ", nb),
             both_DE = sum(a & b), only_first = sum(a & !b),
             only_second = sum(!a & b), both_EE = sum(!a & !b),
             jaccard = sum(a & b) / max(sum(a | b), 1))
}
agreement <- rbind(
  agree(pig_calls, gg_calls,    "PIG_pergene", "gaga_fitGG"),
  agree(pig_calls, limma_calls, "PIG_pergene", "limma_BH"),
  agree(gg_calls,  limma_calls, "gaga_fitGG",  "limma_BH"))

# Posterior-prob correlation between PIG-pergene and gaga
cor_score <- cor(pig_score, gg_score, method = "spearman")
cor_pp    <- cor(1 - pig$pp[, 1], 1 - gg_fit$pp[, 1], method = "spearman")

summary_df <- data.frame(
  method = c("PIG_pergene", "gaga_fitGG", "limma_BH"),
  n_DE = c(sum(pig_calls), sum(gg_calls), sum(limma_calls)),
  elapsed_sec = c(elapsed_pig, NA, NA),
  spearman_with_gaga_score = c(cor_score, 1, NA),
  stringsAsFactors = FALSE)

# Save outputs
suffix <- paste("pergene", data_choice, transform_choice, sep = "_")
write.csv(summary_df,
  result_file(paste0("pig_armstrong_summary_", suffix, ".csv")),
  row.names = FALSE)
write.csv(agreement,
  result_file(paste0("pig_armstrong_agreement_", suffix, ".csv")),
  row.names = FALSE)
saveRDS(list(pig = pig, gg_calls = gg_calls, gg_score = gg_score,
             limma_calls = limma_calls, limma_score = limma_score,
             pig_calls = pig_calls, pig_score = pig_score,
             groups = groups, hyper = hyper),
        result_file(paste0("pig_armstrong_", suffix, ".rds")))

cat("\n=== Per-gene PIG-GaGa Armstrong summary ===\n")
print(summary_df, row.names = FALSE)
cat("\n=== Agreement (#DE overlap) ===\n")
print(agreement, row.names = FALSE)
cat(sprintf("\nSpearman(PIG-pergene posterior, gaga posterior) = %.4f\n",
            cor_score))
cat("\nOutputs saved to ", results_dir, " with suffix ", suffix, "\n",
    sep = "")
