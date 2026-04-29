#!/usr/bin/env Rscript
# =============================================================================
# 04_pig_pergene_maqc.R
#
# Per-gene PIG-EM applied to MAQC Section 7.  Identical model and same
# log2-scale data as gaga::fitGG; only the per-gene marginal-likelihood
# computation switches from Stirling to PIG-EM + Laplace.
#
# Bug fixes vs 02_pig_gaga_maqc.R:
#   * groups factor has levels in titration order A, C, D, B; pattern matrix
#     columns are in the SAME order.  We pass the factor unchanged to gaga
#     so the hypotheses align.
#   * Inverse-Gamma prior on λ uses rate α_0/ν (matching gaga::simGG source).
#   * Both gaga and PIG see the same log2-scale matrix.
# =============================================================================

suppressMessages({
  library(gaga)
  library(EBarrays)
  library(limma)
})
src_dir <- normalizePath(dirname(sub("^--file=", "",
  commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][1])))
source(file.path(src_dir, "pig_gaga_pergene.R"))

bench <- normalizePath(file.path(src_dir, "..", ".."))
results_dir <- file.path(src_dir, "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
result_file <- function(x) file.path(results_dir, x)

cache_file <- file.path(bench, "Claude", "results", "maqc_loaded.rds")
if (!file.exists(cache_file))
  stop("Run gaga_benchmark/Claude/07a_maqc_cache.R first.")
loaded <- readRDS(cache_file)

n_genes <- as.integer(Sys.getenv("MAQC_N_GENES", unset = "8000"))
em_iter <- as.integer(Sys.getenv("PIG_EM_ITER", unset = "20"))
seed    <- as.integer(Sys.getenv("PIG_GAGA_SEED", unset = "42"))

affy_meta <- loaded$affy_meta; affy <- loaded$affy
qpcr_meta <- loaded$qpcr_meta; qpcr <- loaded$qpcr
gpl570 <- loaded$gpl570; gpl4097 <- loaded$gpl4097

first_site <- grepl("^MAQC_AFX_1_[ABCD][1-5]$", affy_meta$title)
affy_meta_1 <- affy_meta[first_site, , drop = FALSE]
affy_1 <- affy[, affy_meta_1$gsm, drop = FALSE]
colnames(affy_1) <- affy_meta_1$title
extract_pool <- function(titles) {
  m <- regmatches(titles, regexec("^MAQC_[A-Z]+_[0-9]+_([ABCD])[0-9]+$",
                                  titles))
  vapply(m, function(x) if (length(x) >= 2) x[2] else NA_character_,
         character(1))
}
groups_affy <- factor(extract_pool(affy_meta_1$title),
                      levels = c("A", "C", "D", "B"))   # titration order

# Mix top-variance + random sample of probes
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
log2_fit <- log2_affy[keep, , drop = FALSE]
log2_fit <- pmax(log2_fit, 1e-6)  # gaga needs strictly positive

cat(sprintf("MAQC per-gene PIG-EM: n_probes=%d  arrays=%d (%s)\n",
            nrow(log2_fit), ncol(log2_fit),
            paste(table(groups_affy), collapse = "+")))

# Patterns: rows are 5 hypotheses, columns A/C/D/B in titration order
patterns <- matrix(c(
  0,0,0,0,    # all equal
  0,1,1,1,    # A vs (C,D,B)
  0,0,1,1,    # (A,C) vs (D,B)
  0,0,0,1,    # (A,C,D) vs B
  0,1,2,3),   # all different
  ncol = 4, byrow = TRUE)
colnames(patterns) <- levels(groups_affy)   # A, C, D, B

# ---------- gaga::fitGG (reference) ----------------------------------------
cat("\n[ref] gaga::fitGG (Stirling EM) ...\n")
# gaga::groups2int compares colnames(patterns) to unique(groups) as a set;
# pass groups as a character vector with values matching colnames(patterns).
groups_chr <- as.character(groups_affy)
gg_fit <- gaga::fitGG(log2_fit, groups_chr, patterns = patterns,
                      equalcv = TRUE, nclust = 1,
                      method = "quickEM", trace = FALSE)
gg_fit <- gaga::parest(gg_fit, x = log2_fit, groups = groups_chr,
                       alpha = 0.05)
gg_par <- gaga::getpar(gg_fit)
gg_pp <- gg_fit$pp
gg_score <- 1 - gg_pp[, 1]
cat(sprintf("[ref] hyper: a0=%.3f  nu=%.3f  b_alpha=%.3f  nu_alpha=%.3g\n",
            gg_par$a0, gg_par$nu, gg_par$balpha, gg_par$nualpha))

# limma F
cat("[ref] limma F ...\n")
fit_l <- limma::eBayes(limma::lmFit(log2_fit,
                                    stats::model.matrix(~ groups_affy)))
tt <- limma::topTable(fit_l, coef = 2:ncol(stats::model.matrix(~groups_affy)),
                      number = Inf, sort.by = "none")
limma_score <- -log10(pmax(tt$P.Value, 1e-300))
names(limma_score) <- rownames(tt)

# ---------- per-gene PIG-EM + Laplace --------------------------------------
hyper <- list(alpha0 = gg_par$a0, nu = gg_par$nu,
              b_alpha = gg_par$balpha, nu_alpha = gg_par$nualpha)
pi_v <- as.numeric(gg_par$probpat)

cat(sprintf("\n[PIG] per-gene PIG-EM + Laplace (em_iter=%d) ...\n", em_iter))
t0 <- Sys.time()
pig <- pig_gaga_pergene(X = log2_fit, group = groups_affy,
                        patterns = patterns, hyper = hyper, pi_v = pi_v,
                        em_iter = em_iter, verbose = TRUE)
elapsed_pig <- as.numeric(Sys.time() - t0, units = "secs")
pig_pp <- pig$pp
pig_score <- 1 - pig_pp[, 1]

cat(sprintf("[PIG]  elapsed=%.1fs   Spearman(PIG, gaga) score = %.4f\n",
            elapsed_pig, cor(pig_score, gg_score, method = "spearman")))

# ---------- qPCR ROC --------------------------------------------------------
qpcr_groups <- factor(extract_pool(qpcr_meta$title),
                      levels = c("A", "C", "D", "B"))
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

# probe → symbol map
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

methods <- list(PIG_pergene = pig_score,
                gaga_fitGG  = gg_score,
                limma_F     = limma_score)
roc_summary <- do.call(rbind, lapply(names(methods), function(m) {
  s <- score_by_symbol(methods[[m]])
  dat <- merge(qpcr_info, s, by = "symbol", sort = FALSE)
  cu <- validation_curve(dat$score, dat$validated)
  data.frame(method = m,
             n_mapped = nrow(dat),
             n_validated = sum(dat$validated),
             cu)
}))

write.csv(roc_summary, result_file("pig_pergene_maqc_roc_summary.csv"),
          row.names = FALSE)
saveRDS(list(pig = pig, gg_pp = gg_pp, limma_score = limma_score,
             roc = roc_summary, hyper = hyper),
        result_file("pig_pergene_maqc.rds"))

cat("\n=== MAQC per-gene PIG-EM ROC summary ===\n")
print(roc_summary, row.names = FALSE)
cat("\nSpearman correlations between methods (probe-level scores):\n")
cat(sprintf("  PIG vs gaga: %.4f\n",
            cor(pig_score, gg_score, method = "spearman")))
cat(sprintf("  PIG vs limma: %.4f\n",
            cor(pig_score, limma_score, method = "spearman")))
cat(sprintf("  gaga vs limma: %.4f\n",
            cor(gg_score, limma_score, method = "spearman")))
