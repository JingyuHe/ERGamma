#!/usr/bin/env Rscript
# =============================================================================
# 06_pig_mcmc_maqc.R
#
# Per-gene PIG-augmented MCMC for MAQC (5 titration patterns over 4 pools),
# comparing to gaga::fitGG (Stirling) and limma F.
# =============================================================================

suppressMessages({ library(gaga); library(EBarrays); library(limma) })
src_dir <- normalizePath(dirname(sub("^--file=", "",
  commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][1])))
source(file.path(src_dir, "pig_gaga_mcmc_lib.R"))

bench <- normalizePath(file.path(src_dir, "..", ".."))
results_dir <- file.path(src_dir, "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
result_file <- function(x) file.path(results_dir, x)

cache_file <- file.path(bench, "Claude", "results", "maqc_loaded.rds")
if (!file.exists(cache_file))
  stop("Run gaga_benchmark/Claude/07a_maqc_cache.R first.")
loaded <- readRDS(cache_file)

n_genes <- as.integer(Sys.getenv("MAQC_N_GENES", unset = "5000"))
n_iter  <- as.integer(Sys.getenv("PIG_MCMC_ITER", unset = "200"))
n_burn  <- as.integer(Sys.getenv("PIG_MCMC_BURN", unset = "60"))
trunc   <- as.integer(Sys.getenv("PIG_MCMC_TRUNC", unset = "20"))
rw_sd   <- as.numeric(Sys.getenv("PIG_MCMC_RWSD", unset = "0.2"))
seed    <- as.integer(Sys.getenv("PIG_MCMC_SEED", unset = "42"))
chkpt_every <- as.integer(Sys.getenv("PIG_MCMC_CHKPT", unset = "5"))

affy_meta <- loaded$affy_meta; affy <- loaded$affy
qpcr_meta <- loaded$qpcr_meta; qpcr <- loaded$qpcr
gpl570 <- loaded$gpl570; gpl4097 <- loaded$gpl4097

extract_pool <- function(t) {
  m <- regmatches(t, regexec("^MAQC_[A-Z]+_[0-9]+_([ABCD])[0-9]+$", t))
  vapply(m, function(x) if (length(x) >= 2) x[2] else NA_character_,
         character(1))
}

first_site <- grepl("^MAQC_AFX_1_[ABCD][1-5]$", affy_meta$title)
am1 <- affy_meta[first_site, , drop = FALSE]
A1 <- affy[, am1$gsm, drop = FALSE]; colnames(A1) <- am1$title
groups_affy <- factor(extract_pool(am1$title), levels = c("A","C","D","B"))

log2_A <- log2(pmax(A1, 1))
set.seed(seed)
if (n_genes >= nrow(A1)) {
  keep <- seq_len(nrow(A1))
} else {
  vars <- apply(log2_A, 1, stats::var)
  n_top <- floor(n_genes / 2); n_rand <- n_genes - n_top
  top_idx <- order(vars, decreasing = TRUE)[seq_len(n_top)]
  rand_idx <- sample(setdiff(seq_len(nrow(A1)), top_idx), n_rand)
  keep <- sort(c(top_idx, rand_idx))
}
log2_fit <- log2_A[keep, , drop = FALSE]
log2_fit <- pmax(log2_fit, 1e-6)

cat(sprintf("MAQC PIG-MCMC: probes=%d arrays=%d  iter=%d burn=%d\n",
            nrow(log2_fit), ncol(log2_fit), n_iter, n_burn))

patterns <- matrix(c(
  0,0,0,0, 0,1,1,1, 0,0,1,1, 0,0,0,1, 0,1,2,3),
  ncol = 4, byrow = TRUE)
colnames(patterns) <- levels(groups_affy)

cat("\n[ref] gaga::fitGG ...\n")
groups_chr <- as.character(groups_affy)
gg_fit <- gaga::fitGG(log2_fit, groups_chr, patterns = patterns, equalcv = TRUE,
                      nclust = 1, method = "quickEM", trace = FALSE)
gg_fit <- gaga::parest(gg_fit, x = log2_fit, groups = groups_chr, alpha = 0.05)
gg_par <- gaga::getpar(gg_fit)
gg_pp <- gg_fit$pp
gg_score <- 1 - gg_pp[, 1]
cat(sprintf("[ref] hyper a0=%.3f nu=%.3f b_alpha=%.3f nu_alpha=%.3g\n",
            gg_par$a0, gg_par$nu, gg_par$balpha, gg_par$nualpha))

fit_l <- limma::eBayes(limma::lmFit(log2_fit,
                                    stats::model.matrix(~ groups_affy)))
tt <- limma::topTable(fit_l, coef = 2:ncol(stats::model.matrix(~groups_affy)),
                      number = Inf, sort.by = "none")
limma_score <- -log10(pmax(tt$P.Value, 1e-300))
names(limma_score) <- rownames(tt)

hyper <- list(alpha0 = gg_par$a0, nu = gg_par$nu,
              b_alpha = gg_par$balpha, nu_alpha = gg_par$nualpha)
pi_v <- as.numeric(gg_par$probpat)

set.seed(seed)
ckpt <- result_file("pig_mcmc_maqc_ckpt.rds")
cat(sprintf("\n[PIG] MCMC ...\n"))
t0 <- Sys.time()
out <- pig_gaga_mcmc(X = log2_fit, group = groups_affy, patterns = patterns,
                     hyper = hyper, pi_v = pi_v,
                     n_iter = n_iter, n_burn = n_burn,
                     trunc = trunc, rw_sd = rw_sd,
                     init_em_iter = 15,
                     checkpoint_file = ckpt,
                     checkpoint_every = chkpt_every,
                     verbose = TRUE)
elapsed_pig <- as.numeric(Sys.time() - t0, units = "secs")
pig_pp <- out$pp; pig_score <- 1 - pig_pp[, 1]
cat(sprintf("[PIG] elapsed=%.1fs  acc=%.3f  median alpha=%.2f  Spearman(PIG, gaga)=%.4f\n",
            elapsed_pig, mean(out$acc_rate), median(out$alpha_mean),
            cor(pig_score, gg_score, method = "spearman")))

# qPCR validation -----------------------------------------------------------
qpcr_groups <- factor(extract_pool(qpcr_meta$title), levels = c("A","C","D","B"))
colnames(qpcr) <- qpcr_meta$title
qpcr_p <- {
  fit_q <- limma::eBayes(limma::lmFit(
    log2(pmax(qpcr, .Machine$double.eps)),
    stats::model.matrix(~ qpcr_groups)))
  tt_q <- limma::topTable(fit_q,
    coef = 2:ncol(stats::model.matrix(~qpcr_groups)),
    number = Inf, sort.by = "none")
  p <- tt_q$P.Value; names(p) <- rownames(tt_q); p
}
qpcr_annot <- gpl4097[, c("ID", "ORF", "AssayID"), drop = FALSE]
names(qpcr_annot) <- c("qpcr_id", "symbol", "assay_id")
qpcr_annot$symbol <- toupper(trimws(qpcr_annot$symbol))
qpcr_info <- merge(
  data.frame(qpcr_id = rownames(qpcr), qpcr_p = qpcr_p,
             stringsAsFactors = FALSE),
  qpcr_annot, by = "qpcr_id", all.x = TRUE, sort = FALSE)
qpcr_info$validated <- qpcr_info$qpcr_p < 0.05

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
  names(scores) <- rownames(log2_fit)
  idx <- match(affy_map$probe_id, names(scores))
  m <- affy_map[!is.na(idx), , drop = FALSE]
  m$score <- scores[idx[!is.na(idx)]]
  stats::aggregate(score ~ symbol, data = m, FUN = max, na.rm = TRUE)
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
  data.frame(auc_scaled = auc, auc_standard = std_auc)
}
methods <- list(PIG_MCMC = pig_score, gaga_fitGG = gg_score,
                limma_F = limma_score)
roc_summary <- do.call(rbind, lapply(names(methods), function(m) {
  s <- score_by_symbol(methods[[m]])
  dat <- merge(qpcr_info, s, by = "symbol", sort = FALSE)
  cu <- validation_curve(dat$score, dat$validated)
  data.frame(method = m, n_mapped = nrow(dat),
             n_validated = sum(dat$validated), cu)
}))

write.csv(roc_summary, result_file("pig_mcmc_maqc_roc_summary.csv"),
          row.names = FALSE)
saveRDS(list(out = out, gg_pp = gg_pp, limma_score = limma_score,
             roc = roc_summary, hyper = hyper),
        result_file("pig_mcmc_maqc.rds"))

cat("\n=== MAQC PIG-MCMC ROC summary ===\n")
print(roc_summary, row.names = FALSE)
cat("\nSpearman correlations (probe-level):\n")
cat(sprintf("  PIG-MCMC vs gaga: %.4f\n",
            cor(pig_score, gg_score, method = "spearman")))
cat(sprintf("  PIG-MCMC vs limma: %.4f\n",
            cor(pig_score, limma_score, method = "spearman")))
cat(sprintf("  gaga vs limma: %.4f\n",
            cor(gg_score, limma_score, method = "spearman")))
