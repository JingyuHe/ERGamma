#!/usr/bin/env Rscript
# Aggregate the existing PIG-slice empirical runs (Armstrong + MAQC) into a
# single comparison table for the paper.

src_dir <- normalizePath(dirname(sub("^--file=", "",
  commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][1])))
results_dir_pig <- normalizePath(file.path(src_dir, "..", "results"))
agg_dir <- file.path(src_dir, "summary")
dir.create(agg_dir, recursive = TRUE, showWarnings = FALSE)
af <- function(f) file.path(agg_dir, f)

# ---------------- Armstrong: Schliep ------------------------------------
arm_schliep <- readRDS(file.path(results_dir_pig,
  "pig_slice_armstrong_slice_schliep_filtered_log.rds"))
arm_orange  <- readRDS(file.path(results_dir_pig,
  "pig_slice_armstrong_slice_orange_full_log2_floor20.rds"))
maqc        <- readRDS(file.path(results_dir_pig, "pig_slice_maqc.rds"))

armstrong_row <- function(label, obj) {
  pp_pig  <- obj$out$pp[, 1]
  pp_gaga <- obj$gg_pp[, 1]
  data.frame(
    dataset = label,
    n_genes = nrow(obj$out$pp),
    n_arrays = length(obj$groups),
    nDE_gaga = sum(obj$gg_calls),
    nDE_pig = sum(obj$pig_calls),
    nDE_limma = sum(obj$limma_calls),
    overlap_pig_gaga = sum(obj$gg_calls & obj$pig_calls),
    only_pig = sum(obj$pig_calls & !obj$gg_calls),
    only_gaga = sum(!obj$pig_calls & obj$gg_calls),
    jaccard_pig_gaga = sum(obj$gg_calls & obj$pig_calls) /
                       max(sum(obj$gg_calls | obj$pig_calls), 1),
    spearman_pp_pig_gaga = stats::cor(pp_pig, pp_gaga, method = "spearman"),
    mean_abs_pp_diff = mean(abs(pp_pig - pp_gaga)),
    max_abs_pp_diff = max(abs(pp_pig - pp_gaga)),
    median_alpha_pig = median(obj$out$alpha_mean),
    stringsAsFactors = FALSE
  )
}

armstrong_tbl <- rbind(
  armstrong_row("Armstrong-Schliep", arm_schliep),
  armstrong_row("Armstrong-Orange", arm_orange))
write.csv(armstrong_tbl, af("empirical_armstrong.csv"), row.names = FALSE)

# ---------------- MAQC ROC summary --------------------------------------
maqc_roc <- maqc$roc
write.csv(maqc_roc, af("empirical_maqc_roc.csv"), row.names = FALSE)

# Score correlations on MAQC
maqc_pig_score <- 1 - maqc$out$pp[, 1]
maqc_gaga_score <- 1 - maqc$gg_pp[, 1]
maqc_extra <- data.frame(
  dataset = "MAQC",
  n_probes = nrow(maqc$out$pp),
  spearman_pp_pig_gaga =
    stats::cor(maqc_pig_score, maqc_gaga_score, method = "spearman"),
  mean_abs_pp_diff = mean(abs(maqc$out$pp[, 1] - maqc$gg_pp[, 1])),
  AUC_PIG  = maqc_roc$auc_standard[maqc_roc$method == "PIG_slice"],
  AUC_gaga = maqc_roc$auc_standard[maqc_roc$method == "gaga_fitGG"],
  AUC_limma = maqc_roc$auc_standard[maqc_roc$method == "limma_F"],
  stringsAsFactors = FALSE)
write.csv(maqc_extra, af("empirical_maqc_score.csv"), row.names = FALSE)

# Combined headline empirical table
empirical_headline <- data.frame(
  dataset = c("Armstrong (Schliep mirror)", "Armstrong (Orange mirror)", "MAQC"),
  n_features = c(nrow(arm_schliep$out$pp),
                 nrow(arm_orange$out$pp),
                 nrow(maqc$out$pp)),
  metric_used = c("#DE @ FDR=0.05",
                  "#DE @ FDR=0.05",
                  "AUC vs qPCR validation"),
  gaga_value = c(as.character(sum(arm_schliep$gg_calls)),
                 as.character(sum(arm_orange$gg_calls)),
                 sprintf("%.4f", maqc_roc$auc_standard[maqc_roc$method == "gaga_fitGG"])),
  pig_value = c(as.character(sum(arm_schliep$pig_calls)),
                as.character(sum(arm_orange$pig_calls)),
                sprintf("%.4f", maqc_roc$auc_standard[maqc_roc$method == "PIG_slice"])),
  spearman_pig_vs_gaga = sprintf("%.4f", c(
    stats::cor(arm_schliep$out$pp[, 1], arm_schliep$gg_pp[, 1], method = "spearman"),
    stats::cor(arm_orange$out$pp[, 1],  arm_orange$gg_pp[, 1], method = "spearman"),
    stats::cor(maqc$out$pp[, 1], maqc$gg_pp[, 1], method = "spearman"))),
  mean_abs_pp_diff = sprintf("%.4f", c(
    mean(abs(arm_schliep$out$pp[, 1] - arm_schliep$gg_pp[, 1])),
    mean(abs(arm_orange$out$pp[, 1] - arm_orange$gg_pp[, 1])),
    mean(abs(maqc$out$pp[, 1] - maqc$gg_pp[, 1])))),
  stringsAsFactors = FALSE
)
write.csv(empirical_headline, af("empirical_headline.csv"), row.names = FALSE)

cat("\n=== EMPIRICAL HEADLINE ===\n")
print(empirical_headline, row.names = FALSE)
cat("\n=== ARMSTRONG details ===\n")
print(armstrong_tbl, row.names = FALSE)
cat("\n=== MAQC details ===\n")
print(maqc_roc, row.names = FALSE)
print(maqc_extra, row.names = FALSE)
