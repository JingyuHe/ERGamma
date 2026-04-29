#!/usr/bin/env Rscript
# =============================================================================
# 01_pig_gaga_armstrong.R
#
# Apply the PIG-augmented Bayesian GaGa sampler to the Armstrong (ALL vs MLL)
# leukemia data and compare the resulting DE calls to:
#   * gaga::fitGG (Stirling-EM, Rossell 2009),
#   * limma BH on the same log-transformed matrix.
#
# Outputs go to gaga_benchmark/R/Claude/results/pig_armstrong_*.csv.
# Default uses the schliep_filtered mirror (2194 probes) for speed.  Set
# ARMSTRONG_DATA=orange_full to run on the 12533-probe matrix instead.
# =============================================================================

suppressMessages({
  library(gaga)
  library(limma)
})
src_dir <- normalizePath(dirname(sub("^--file=", "",
  commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][1])))
source(file.path(src_dir, "pig_gaga_lib.R"))

bench <- normalizePath(file.path(src_dir, "..", ".."))
results_dir <- file.path(src_dir, "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
result_file <- function(x) file.path(results_dir, x)

data_choice <- Sys.getenv("ARMSTRONG_DATA", unset = "schliep_filtered")
n_iter <- as.integer(Sys.getenv("PIG_GAGA_ITER", unset = "150"))
n_burn <- as.integer(Sys.getenv("PIG_GAGA_BURN", unset = "50"))
N_pig  <- as.integer(Sys.getenv("PIG_GAGA_NPIG", unset = "80"))
seed   <- as.integer(Sys.getenv("PIG_GAGA_SEED", unset = "42"))
delta_pri  <- as.numeric(Sys.getenv("PIG_GAGA_DELTA",  unset = "1"))
delta0_pri <- as.numeric(Sys.getenv("PIG_GAGA_DELTA0", unset = "1"))

# ---------- data loaders ----------------------------------------------------
load_armstrong_schliep <- function(data_file) {
  raw <- read.delim(data_file, check.names = FALSE, stringsAsFactors = FALSE)
  groups <- as.character(unlist(raw[1, -1]))
  x <- as.matrix(data.frame(lapply(raw[-1, -1, drop = FALSE], as.numeric),
                            check.names = FALSE))
  rownames(x) <- raw[-1, 1]
  list(x = x, groups = groups)
}
load_armstrong_orange <- function(data_file) {
  raw <- read.delim(data_file, check.names = FALSE, stringsAsFactors = FALSE)
  groups <- raw$class[-c(1, 2)]
  xdf <- raw[-c(1, 2), -1, drop = FALSE]
  x <- t(apply(xdf, 2, as.numeric))
  rownames(x) <- colnames(raw)[-1]
  list(x = x, groups = groups)
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

# Restrict to ALL + MLL, drop the last 2 MLL (the U95Av2 ones)
keep <- dat$groups %in% c("ALL", "MLL")
mll_idx <- which(dat$groups == "MLL")
keep[tail(mll_idx, 2)] <- FALSE
X_raw <- dat$x[, keep, drop = FALSE]
groups <- factor(dat$groups[keep], levels = c("ALL", "MLL"))

apply_t <- function(X, t) {
  if (t == "log")           list(x_log = log(pmax(X, 1)),         x_pos = pmax(X, 1))
  else if (t == "log2_floor20") list(x_log = log2(pmax(X, 20)),       x_pos = pmax(X, 20))
  else if (t == "log2_clip")    list(x_log = log2(pmax(X, 1) + 1),    x_pos = pmax(X, 1) + 1)
  else stop("unknown transform: ", t)
}
tx <- apply_t(X_raw, transform_choice)
# Use the positive-scale version for GaGa (gamma sampling needs positive values)
X <- tx$x_pos

cat(sprintf("Armstrong PIG-GaGa: data=%s transform=%s genes=%d samples=%d\n",
            data_choice, transform_choice, nrow(X), ncol(X)))
cat(sprintf("groups: ALL=%d MLL=%d\n",
            sum(groups == "ALL"), sum(groups == "MLL")))

# ---------- run gaga::fitGG (Stirling EM) for reference ---------------------
patterns <- matrix(c(0, 0, 0, 1), 2, 2)
colnames(patterns) <- levels(groups)
cat("\n[ref] Fitting gaga::fitGG (Stirling EM) ...\n")
gg_fit <- gaga::fitGG(tx$x_log, groups, patterns = patterns,
                      equalcv = TRUE, nclust = 1,
                      method = "quickEM", trace = FALSE)
gg_fit <- gaga::parest(gg_fit, x = tx$x_log, groups = groups, alpha = 0.05)
gg_genes <- gaga::findgenes(gg_fit, tx$x_log, groups, fdrmax = 0.05,
                            parametric = TRUE)
gg_calls <- as.logical(gg_genes$d != 0)
gg_score <- 1 - gg_fit$pp[, 1]
gg_par <- gaga::getpar(gg_fit)
cat(sprintf("[ref] gaga::fitGG  alpha=%.3f  alpha0=%.3f  nu=%.3f  #DE=%d\n",
            gg_par$balpha, gg_par$a0, gg_par$nu, sum(gg_calls)))

# ---------- run limma_BH ----------------------------------------------------
cat("\n[ref] Running limma_BH on log-scale data ...\n")
design <- stats::model.matrix(~ groups)
fit_l <- limma::eBayes(limma::lmFit(tx$x_log, design))
limma_p <- fit_l$p.value[, 2]
limma_calls <- p.adjust(limma_p, "BH") <= 0.05
limma_score <- -log10(pmax(limma_p, 1e-300))
cat(sprintf("[ref] limma_BH  #DE=%d\n", sum(limma_calls)))

# ---------- run PIG-GaGa Gibbs ---------------------------------------------
cat(sprintf("\n[PIG] Running PIG-GaGa Gibbs (iter=%d burn=%d N_pig=%d) ...\n",
            n_iter, n_burn, N_pig))
set.seed(seed)
# Use Stirling-EM hyperparameter estimates as Damsleth prior centres so the
# PIG sampler starts in a sensible region.
mu_pri  <- max(gg_par$balpha + 1.5, 2)
mu0_pri <- max(gg_par$a0 + 1.5, 2)
init <- list(alpha = gg_par$balpha, alpha0 = gg_par$a0,
             lambda0 = gg_par$nu, pi = c(1 - gg_par$probpat[2],
                                         gg_par$probpat[2]))

t0 <- Sys.time()
out <- pig_gaga_gibbs(
  X = X, group = groups, patterns = patterns,
  pi_prior = c(1, 1),
  delta = delta_pri,  mu = mu_pri,  eta_prior = 0,
  delta0 = delta0_pri, mu0 = mu0_pri, eta_prior0 = 0,
  init = init,
  n_iter = n_iter, n_burn = n_burn,
  N_pig = N_pig, verbose = TRUE,
  fix_alpha0 = FALSE, fix_lambda0 = FALSE)
elapsed_pig <- as.numeric(Sys.time() - t0, units = "secs")

pig_score <- de_score_from_pp(out$pp)
pig_calls <- fdr_calls(out$pp, 0.05)
cat(sprintf("\n[PIG] PIG-GaGa  alpha=%.3f (mean post burn-in)  alpha0=%.3f  ",
            mean(tail(out$alpha, n_iter - n_burn)),
            mean(tail(out$alpha0, n_iter - n_burn))),
    sprintf("lambda0=%.3f  #DE=%d  elapsed=%.1fs\n",
            mean(tail(out$lambda0, n_iter - n_burn)),
            sum(pig_calls), elapsed_pig))

# ---------- agreement between methods --------------------------------------
agree_table <- function(a, b, na, nb) {
  data.frame(method_pair = paste0(na, " vs ", nb),
             both_DE     = sum( a &  b),
             only_first  = sum( a & !b),
             only_second = sum(!a &  b),
             both_EE     = sum(!a & !b),
             jaccard     = sum( a &  b) / max(sum(a | b), 1))
}
agree <- rbind(
  agree_table(pig_calls, gg_calls,    "PIG_GaGa", "gaga_fitGG"),
  agree_table(pig_calls, limma_calls, "PIG_GaGa", "limma_BH"),
  agree_table(gg_calls,  limma_calls, "gaga_fitGG", "limma_BH"))

summary_df <- data.frame(
  method = c("PIG_GaGa", "gaga_fitGG", "limma_BH"),
  n_DE = c(sum(pig_calls), sum(gg_calls), sum(limma_calls)),
  elapsed_sec = c(elapsed_pig, NA, NA),
  stringsAsFactors = FALSE)
trace_df <- data.frame(
  iter = seq_len(n_iter), alpha = out$alpha,
  alpha0 = out$alpha0, lambda0 = out$lambda0)

# ---------- save outputs ----------------------------------------------------
suffix <- paste(data_choice, transform_choice, sep = "_")
write.csv(summary_df,
  result_file(paste0("pig_armstrong_summary_", suffix, ".csv")),
  row.names = FALSE)
write.csv(agree,
  result_file(paste0("pig_armstrong_agreement_", suffix, ".csv")),
  row.names = FALSE)
write.csv(trace_df,
  result_file(paste0("pig_armstrong_trace_", suffix, ".csv")),
  row.names = FALSE)
saveRDS(list(out = out, gg_calls = gg_calls, gg_score = gg_score,
             limma_calls = limma_calls, limma_score = limma_score,
             pig_calls = pig_calls, pig_score = pig_score,
             groups = groups, transform = transform_choice),
        result_file(paste0("pig_armstrong_", suffix, ".rds")))

cat("\n=== PIG-GaGa Armstrong summary ===\n")
print(summary_df, row.names = FALSE)
cat("\n=== Method-pair agreement (#DE overlap) ===\n")
print(agree, row.names = FALSE)
cat("\nOutputs saved to ", results_dir, " with suffix ", suffix, "\n", sep = "")
