args <- commandArgs(FALSE)
this_file <- sub("^--file=", "", args[grepl("^--file=", args)][1])
src_dir <- dirname(normalizePath(this_file, mustWork = TRUE))
source(file.path(src_dir, "quad_gaga_lib.R"))

bench <- setup_paths()
out_dir <- file.path(src_dir, "results")
root_result_dir <- file.path(bench, "results")

read_csv_if <- function(path) {
  if (!file.exists(path)) {
    stop("Missing required result file: ", path)
  }
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

fmt <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", formatC(as.numeric(x), format = "f", digits = digits))
}

md_table <- function(df) {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  if (!nrow(df)) {
    return("")
  }
  df[] <- lapply(df, as.character)
  header <- paste0("| ", paste(names(df), collapse = " | "), " |")
  sep <- paste0("| ", paste(rep("---", ncol(df)), collapse = " | "), " |")
  rows <- apply(df, 1, function(z) paste0("| ", paste(z, collapse = " | "), " |"))
  paste(c(header, sep, rows), collapse = "\n")
}

copy_plot_path <- function(name) {
  file.path(out_dir, name)
}

plot_armstrong_scatter <- function(scores, summary, pdf_file, png_file) {
  catg <- ifelse(scores$gaga_call & scores$pig_mcmc_call, "both",
                 ifelse(scores$gaga_call, "GaGa only",
                        ifelse(scores$pig_mcmc_call, "PIG only", "neither")))
  cols <- c("both" = "#1B9E77", "GaGa only" = "#D95F02",
            "PIG only" = "#7570B3", "neither" = "grey75")
  draw <- function() {
    plot(
      scores$gaga_score,
      scores$pig_mcmc_score,
      pch = 20,
      cex = 0.45,
      col = cols[catg],
      xlab = "Original GaGa DE score, 1 - pp(null)",
      ylab = "PIG-MCMC GaGa DE score, 1 - pp(null)",
      main = "Armstrong ALL vs MLL: posterior scores"
    )
    abline(0, 1, col = "black", lty = 2)
    legend("bottomright", legend = names(cols), col = cols, pch = 20, bty = "n")
    mtext(
      paste0(
        "Jaccard = ", fmt(summary$jaccard_with_other[1], 3),
        ", Spearman = ", fmt(summary$spearman_score_vs_other[1], 3)
      ),
      side = 3,
      line = 0.1,
      cex = 0.8
    )
  }
  grDevices::pdf(pdf_file, width = 6.5, height = 6)
  draw()
  grDevices::dev.off()
  grDevices::png(png_file, width = 1500, height = 1400, res = 220)
  draw()
  grDevices::dev.off()
}

plot_maqc_patterns <- function(pattern_counts, pdf_file, png_file) {
  counts <- xtabs(count ~ method + pattern, data = pattern_counts)
  draw <- function() {
    barplot(
      counts,
      beside = TRUE,
      col = c("#1B9E77", "#D95F02"),
      xlab = "Pattern",
      ylab = "Number of probes",
      main = "MAQC qPCR-mapped probes: most probable pattern"
    )
    legend("topright", legend = rownames(counts),
           fill = c("#1B9E77", "#D95F02"), bty = "n")
  }
  grDevices::pdf(pdf_file, width = 7, height = 5)
  draw()
  grDevices::dev.off()
  grDevices::png(png_file, width = 1500, height = 1100, res = 220)
  draw()
  grDevices::dev.off()
}

plot_maqc_curve <- function(curves, pdf_file, png_file) {
  cols <- c("GaGa" = "#1B9E77", "PIG_MCMC_GaGa" = "#D95F02",
            "limma_F" = "#7570B3")
  draw <- function() {
    plot(
      NA, NA,
      xlim = c(0, max(curves$prop_nonvalidated, na.rm = TRUE)),
      ylim = c(0, max(curves$prop_validated, na.rm = TRUE)),
      xlab = "Proportion non-validated qPCR assays",
      ylab = "Proportion validated qPCR assays",
      main = "MAQC validation curve"
    )
    abline(0, 1, col = "grey80", lty = 2)
    for (m in unique(curves$method)) {
      cur <- curves[curves$method == m, ]
      lines(cur$prop_nonvalidated, cur$prop_validated, col = cols[[m]], lwd = 2)
    }
    legend("bottomright", legend = unique(curves$method),
           col = cols[unique(curves$method)], lwd = 2, bty = "n")
  }
  grDevices::pdf(pdf_file, width = 6, height = 5)
  draw()
  grDevices::dev.off()
  grDevices::png(png_file, width = 1500, height = 1200, res = 220)
  draw()
  grDevices::dev.off()
}

arm_table_path <- file.path(
  bench,
  "Claude", "results",
  "armstrong_table1_summary_schliep_filtered_log_equalcv_exclu95av2.csv"
)
if (!file.exists(arm_table_path)) {
  arm_table_path <- file.path(
    root_result_dir,
    "armstrong_table1_summary_schliep_filtered_log2_clip_exclude2mll.csv"
  )
}
arm_table_schliep <- read_csv_if(arm_table_path)
arm_fig_schliep <- read_csv_if(file.path(
  root_result_dir,
  "armstrong_figure1_summary_schliep_filtered_log.csv"
))
arm_pig_summary <- read_csv_if(file.path(out_dir, "armstrong_pig_mcmc_gaga_summary.csv"))
arm_pig_scores <- read_csv_if(file.path(out_dir, "armstrong_pig_mcmc_gaga_scores.csv"))

maqc_orig_summary <- read_csv_if(file.path(root_result_dir, "maqc_roc_summary.csv"))
maqc_orig_patterns <- read_csv_if(file.path(root_result_dir, "maqc_pattern_counts.csv"))
maqc_pig_summary <- read_csv_if(file.path(out_dir, "maqc_pig_mcmc_gaga_validation_summary.csv"))
maqc_pig_patterns <- read_csv_if(file.path(out_dir, "maqc_pig_mcmc_gaga_pattern_counts.csv"))
maqc_pig_curves <- read_csv_if(file.path(out_dir, "maqc_pig_mcmc_gaga_validation_curve.csv"))

plot_armstrong_scatter(
  arm_pig_scores,
  arm_pig_summary,
  copy_plot_path("empirical_armstrong_gaga_vs_pig_scores.pdf"),
  copy_plot_path("empirical_armstrong_gaga_vs_pig_scores.png")
)
plot_maqc_patterns(
  maqc_pig_patterns,
  copy_plot_path("empirical_maqc_pattern_counts_gaga_vs_pig.pdf"),
  copy_plot_path("empirical_maqc_pattern_counts_gaga_vs_pig.png")
)
plot_maqc_curve(
  maqc_pig_curves,
  copy_plot_path("empirical_maqc_validation_curve_gaga_vs_pig.pdf"),
  copy_plot_path("empirical_maqc_validation_curve_gaga_vs_pig.png")
)

arm_table_show <- arm_table_schliep[
  arm_table_schliep$method %in% c("GaGa", "MiGaGa2", "limma_BH") &
    arm_table_schliep$n_per_group %in% c("5", "10", "15", "All data"),
  intersect(
    c("method", "n_per_group", "n_de", "paper_n_de",
      "reproducibility", "paper_repro", "n_success"),
    names(arm_table_schliep)
  )
]
arm_table_show$n_de <- fmt(arm_table_show$n_de, 1)
if ("paper_n_de" %in% names(arm_table_show)) {
  arm_table_show$paper_n_de <- fmt(arm_table_show$paper_n_de, 1)
}
arm_table_show$reproducibility <- fmt(arm_table_show$reproducibility, 3)
if ("paper_repro" %in% names(arm_table_show)) {
  arm_table_show$paper_repro <- fmt(arm_table_show$paper_repro, 3)
}

arm_pig_show <- arm_pig_summary[, c(
  "method", "n_de", "overlap_with_other", "jaccard_with_other",
  "spearman_score_vs_other", "mean_abs_score_diff", "elapsed_seconds",
  "n_iter", "burnin", "thin", "trunc"
)]
arm_pig_show$jaccard_with_other <- fmt(arm_pig_show$jaccard_with_other, 3)
arm_pig_show$spearman_score_vs_other <- fmt(arm_pig_show$spearman_score_vs_other, 3)
arm_pig_show$mean_abs_score_diff <- fmt(arm_pig_show$mean_abs_score_diff, 4)
arm_pig_show$elapsed_seconds <- fmt(arm_pig_show$elapsed_seconds, 1)

maqc_orig_show <- maqc_orig_summary[, c(
  "method", "n_mapped_assays", "n_validated_mapped", "auc_standard"
)]
maqc_orig_show$auc_standard <- fmt(maqc_orig_show$auc_standard, 3)

maqc_pig_show <- maqc_pig_summary[, c(
  "method", "n_fit_probes", "n_mapped_assays", "n_validated_mapped",
  "auc_standard", "spearman_probe_score_vs_gaga", "n_de_at_fdr",
  "elapsed_pig_seconds", "n_iter", "burnin", "thin", "trunc"
)]
maqc_pig_show$auc_standard <- fmt(maqc_pig_show$auc_standard, 3)
maqc_pig_show$spearman_probe_score_vs_gaga <- fmt(maqc_pig_show$spearman_probe_score_vs_gaga, 3)
maqc_pig_show$elapsed_pig_seconds <- fmt(maqc_pig_show$elapsed_pig_seconds, 1)

maqc_patterns_show <- maqc_pig_patterns
maqc_patterns_show$count <- as.integer(maqc_patterns_show$count)

paper_pattern_show <- maqc_orig_patterns[
  maqc_orig_patterns$method %in% c("GaGa", "MiGaGa2"),
  c("method", "pattern", "count", "paper_count", "diff_from_paper")
]

lines <- c(
  "# Empirical Results: Original GaGa and PIG-MCMC GaGa",
  "",
  "This report is organized in the same order as the empirical part of the GaGa paper: Armstrong leukemia data first, then MAQC validation.",
  "",
  "## 1. Armstrong ALL vs MLL",
  "",
  "### 1.1 Original GaGa Reproduction",
  "",
  paste0(
    "Data: `", arm_fig_schliep$data[1], "`, transform: `",
    arm_fig_schliep$transform[1], "`, genes: ", arm_fig_schliep$n_genes[1],
    ", samples: ", arm_fig_schliep$n_samples[1], " (ALL ",
    arm_fig_schliep$n_all[1], ", MLL ", arm_fig_schliep$n_mll[1], ")."
  ),
  "",
  "Figure 1 reproduction:",
  "",
  paste0("![Armstrong Figure 1 reproduction](", normalizePath(file.path(root_result_dir, "armstrong_figure1_reproduction_schliep_filtered_log.png"), mustWork = FALSE), ")"),
  "",
  "Table 1-style reproducibility summary. When paper columns are present, they are the values reported by Rossell/GaGa and the remaining columns are our reproduction:",
  "",
  md_table(arm_table_show),
  "",
  "### 1.2 PIG-MCMC GaGa on Armstrong",
  "",
  "PIG-MCMC uses the original GaGa empirical-Bayes hyperparameters and replaces the posterior scoring layer with the PIG-augmented sampler.",
  "",
  md_table(arm_pig_show),
  "",
  paste0("![Armstrong score scatter](", normalizePath(copy_plot_path("empirical_armstrong_gaga_vs_pig_scores.png"), mustWork = FALSE), ")"),
  "",
  "## 2. MAQC Validation",
  "",
  "### 2.1 Original GaGa Reproduction",
  "",
  "Original MAQC validation curve summary from the reproduced GaGa pipeline:",
  "",
  md_table(maqc_orig_show),
  "",
  "Pattern counts versus the counts reported in the GaGa paper:",
  "",
  md_table(paper_pattern_show),
  "",
  paste0("![Original MAQC ROC reproduction](", normalizePath(file.path(root_result_dir, "maqc_roc_curves.png"), mustWork = FALSE), ")"),
  "",
  "### 2.2 PIG-MCMC GaGa on qPCR-Mapped MAQC Probes",
  "",
  "For PIG-MCMC, GaGa hyperparameters are fitted on all first-site Affymetrix probes, matching the original empirical-Bayes step; the PIG scorer is then run on all qPCR-mapped probes used by the validation analysis.",
  "",
  md_table(maqc_pig_show),
  "",
  "Most probable pattern counts on this qPCR-mapped probe set:",
  "",
  md_table(maqc_patterns_show),
  "",
  paste0("![MAQC validation curve](", normalizePath(copy_plot_path("empirical_maqc_validation_curve_gaga_vs_pig.png"), mustWork = FALSE), ")"),
  "",
  paste0("![MAQC pattern counts](", normalizePath(copy_plot_path("empirical_maqc_pattern_counts_gaga_vs_pig.png"), mustWork = FALSE), ")"),
  "",
  "## Interpretation",
  "",
  "- Armstrong: PIG-MCMC and original GaGa make essentially the same calls on this data set; the selected sets have high overlap and nearly identical posterior score ranks.",
  "- MAQC: PIG-MCMC remains close to original GaGa but is slightly more conservative on the qPCR-mapped probe set. The validation AUC is modestly lower than original GaGa and limma in this run.",
  "- Together with the oracle simulation, this supports a focused claim: PIG-MCMC fixes the low-shape/Stirling failure mode without materially changing ordinary empirical examples where the original approximation is already stable.",
  ""
)

report_file <- file.path(out_dir, "EMPIRICAL_RESULTS.md")
writeLines(lines, report_file)

cat("Wrote empirical report:\n")
cat(report_file, "\n")
cat("Key outputs:\n")
cat(copy_plot_path("empirical_armstrong_gaga_vs_pig_scores.png"), "\n")
cat(copy_plot_path("empirical_maqc_validation_curve_gaga_vs_pig.png"), "\n")
cat(copy_plot_path("empirical_maqc_pattern_counts_gaga_vs_pig.png"), "\n")
