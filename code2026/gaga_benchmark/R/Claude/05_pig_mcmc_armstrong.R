#!/usr/bin/env Rscript
# =============================================================================
# 05_pig_mcmc_armstrong.R
#
# Per-gene PIG-augmented MCMC for Armstrong, comparing posterior pattern
# probabilities to gaga::fitGG (Stirling) and to limma_BH.
# Uses the exact target chain from PIG_GAGA_MCMC_FORMULAS.md plus a log-α
# RW fallback for mixing.
# =============================================================================

suppressMessages({ library(gaga); library(limma) })
src_dir <- normalizePath(dirname(sub("^--file=", "",
  commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][1])))
source(file.path(src_dir, "pig_gaga_mcmc_lib.R"))

bench <- normalizePath(file.path(src_dir, "..", ".."))
results_dir <- file.path(src_dir, "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
result_file <- function(x) file.path(results_dir, x)

data_choice <- Sys.getenv("ARMSTRONG_DATA", unset = "schliep_filtered")
n_iter <- as.integer(Sys.getenv("PIG_MCMC_ITER", unset = "300"))
n_burn <- as.integer(Sys.getenv("PIG_MCMC_BURN", unset = "100"))
trunc  <- as.integer(Sys.getenv("PIG_MCMC_TRUNC", unset = "30"))
rw_sd  <- as.numeric(Sys.getenv("PIG_MCMC_RWSD", unset = "0.2"))
seed   <- as.integer(Sys.getenv("PIG_MCMC_SEED", unset = "42"))

# loaders --------------------------------------------------------------------
load_armstrong_schliep <- function(f) {
  raw <- read.delim(f, check.names = FALSE, stringsAsFactors = FALSE)
  groups <- as.character(unlist(raw[1, -1]))
  x <- as.matrix(data.frame(lapply(raw[-1, -1, drop = FALSE], as.numeric),
                            check.names = FALSE))
  rownames(x) <- raw[-1, 1]; list(x = x, groups = groups)
}
load_armstrong_orange <- function(f) {
  raw <- read.delim(f, check.names = FALSE, stringsAsFactors = FALSE)
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
  if (t == "log")               log(pmax(X, 1))
  else if (t == "log2_floor20") log2(pmax(X, 20))
  else if (t == "log2_clip")    log2(pmax(X, 1) + 1)
  else stop("unknown transform")
}
X_log <- log_t(X_raw, transform_choice)
X     <- X_log

cat(sprintf("Armstrong PIG-MCMC: data=%s transform=%s genes=%d samples=%d  (n_iter=%d burn=%d)\n",
            data_choice, transform_choice, nrow(X), ncol(X), n_iter, n_burn))

patterns <- matrix(c(0, 0, 0, 1), 2, 2)
colnames(patterns) <- levels(groups)

# gaga reference --------------------------------------------------------------
cat("\n[ref] gaga::fitGG (Stirling) ...\n")
gg_fit <- gaga::fitGG(X_log, groups, patterns = patterns, equalcv = TRUE,
                      nclust = 1, method = "quickEM", trace = FALSE)
gg_fit <- gaga::parest(gg_fit, x = X_log, groups = groups, alpha = 0.05)
gg_par <- gaga::getpar(gg_fit)
gg_genes <- gaga::findgenes(gg_fit, X_log, groups, fdrmax = 0.05,
                            parametric = TRUE)
gg_calls <- as.logical(gg_genes$d != 0)
gg_score <- 1 - gg_fit$pp[, 1]
cat(sprintf("[ref] hyper a0=%.3f nu=%.3f b_alpha=%.3f nu_alpha=%.3g  pi=(%.3f,%.3f)  #DE=%d\n",
            gg_par$a0, gg_par$nu, gg_par$balpha, gg_par$nualpha,
            gg_par$probpat[1], gg_par$probpat[2], sum(gg_calls)))

# limma_BH --------------------------------------------------------------------
fit_l <- limma::eBayes(limma::lmFit(X_log, stats::model.matrix(~ groups)))
limma_p <- fit_l$p.value[, 2]
limma_calls <- p.adjust(limma_p, "BH") <= 0.05
limma_score <- -log10(pmax(limma_p, 1e-300))
cat(sprintf("[ref] limma_BH #DE = %d\n", sum(limma_calls)))

# PIG-MCMC --------------------------------------------------------------------
hyper <- list(alpha0 = gg_par$a0, nu = gg_par$nu,
              b_alpha = gg_par$balpha, nu_alpha = gg_par$nualpha)
pi_v <- as.numeric(gg_par$probpat)

set.seed(seed)
ckpt <- result_file(paste0("pig_mcmc_armstrong_",
                           data_choice, "_", transform_choice,
                           "_ckpt.rds"))
cat(sprintf("\n[PIG] MCMC (n_iter=%d burn=%d trunc=%d rw_sd=%.2f) ...\n",
            n_iter, n_burn, trunc, rw_sd))
t0 <- Sys.time()
out <- pig_gaga_mcmc(X = X, group = groups, patterns = patterns,
                     hyper = hyper, pi_v = pi_v,
                     n_iter = n_iter, n_burn = n_burn,
                     trunc = trunc, rw_sd = rw_sd,
                     init_em_iter = 15,
                     checkpoint_file = ckpt,
                     checkpoint_every = as.integer(Sys.getenv("PIG_MCMC_CHKPT", unset = "10")),
                     verbose = TRUE)
elapsed_pig <- as.numeric(Sys.time() - t0, units = "secs")
pig_score <- 1 - out$pp[, 1]
pig_calls <- fdr_calls_pp(out$pp, 0.05)
cat(sprintf("[PIG]  elapsed=%.1fs  acc=%.3f  median alpha=%.2f  #DE=%d\n",
            elapsed_pig, mean(out$acc_rate), median(out$alpha_mean),
            sum(pig_calls)))

# agreement ------------------------------------------------------------------
agree <- function(a, b, na, nb) data.frame(
  pair = paste0(na, " vs ", nb),
  both_DE = sum(a & b), only_first = sum(a & !b),
  only_second = sum(!a & b), both_EE = sum(!a & !b),
  jaccard = sum(a & b) / max(sum(a | b), 1))
agreement <- rbind(
  agree(pig_calls, gg_calls, "PIG_MCMC", "gaga_fitGG"),
  agree(pig_calls, limma_calls, "PIG_MCMC", "limma_BH"),
  agree(gg_calls, limma_calls, "gaga_fitGG", "limma_BH"))

summary_df <- data.frame(
  method = c("PIG_MCMC", "gaga_fitGG", "limma_BH"),
  n_DE = c(sum(pig_calls), sum(gg_calls), sum(limma_calls)),
  elapsed_sec = c(elapsed_pig, NA, NA),
  spearman_with_gaga_score = c(
    cor(pig_score, gg_score, method = "spearman"),
    1.0,
    cor(limma_score, gg_score, method = "spearman")),
  stringsAsFactors = FALSE)

suffix <- paste("mcmc", data_choice, transform_choice, sep = "_")
write.csv(summary_df,
  result_file(paste0("pig_mcmc_armstrong_summary_", suffix, ".csv")),
  row.names = FALSE)
write.csv(agreement,
  result_file(paste0("pig_mcmc_armstrong_agreement_", suffix, ".csv")),
  row.names = FALSE)
saveRDS(list(out = out, gg_calls = gg_calls, gg_score = gg_score,
             pig_calls = pig_calls, pig_score = pig_score,
             limma_calls = limma_calls, limma_score = limma_score,
             groups = groups, hyper = hyper, n_iter = n_iter, n_burn = n_burn),
        result_file(paste0("pig_mcmc_armstrong_", suffix, ".rds")))

cat("\n=== PIG-MCMC Armstrong summary ===\n")
print(summary_df, row.names = FALSE)
cat("\n=== Agreement ===\n")
print(agreement, row.names = FALSE)
cat(sprintf("\nSpearman(PIG-MCMC, gaga) on raw posterior probs = %.4f\n",
            cor(pig_score, gg_score, method = "spearman")))
cat("\nOutputs in", results_dir, "with suffix", suffix, "\n")
