args <- commandArgs(FALSE)
this_file <- sub("^--file=", "", args[grepl("^--file=", args)][1])
src_dir <- dirname(normalizePath(this_file, mustWork = TRUE))
source(file.path(src_dir, "quad_gaga_lib.R"))

setup_paths()

result_dir <- file.path(src_dir, "results")
prefixes <- strsplit(
  Sys.getenv(
    "ORACLE_PAPER_PREFIXES",
    "oracle_paper_S0,oracle_paper_S1,oracle_paper_S2,oracle_paper_S3,oracle_paper_S4,oracle_paper_S5"
  ),
  ",",
  fixed = TRUE
)[[1]]
prefixes <- trimws(prefixes)
combined_prefix <- Sys.getenv("ORACLE_PAPER_COMBINED_PREFIX",
                              unset = "oracle_paper_combined")

read_kind <- function(kind) {
  files <- file.path(result_dir, paste0(prefixes, "_", kind, ".csv"))
  missing <- files[!file.exists(files)]
  if (length(missing)) {
    stop("Missing result files: ", paste(missing, collapse = ", "))
  }
  do.call(rbind, lapply(files, read.csv, check.names = FALSE))
}

write_kind <- function(x, kind) {
  write.csv(x, file.path(result_dir, paste0(combined_prefix, "_", kind, ".csv")),
            row.names = FALSE)
}

paper_summary <- function(data, by, vars) {
  key <- interaction(data[by], drop = TRUE, lex.order = TRUE, sep = "\r")
  pieces <- lapply(vars, function(v) {
    rows <- lapply(split(seq_len(nrow(data)), key), function(idx) {
      z <- data[[v]][idx]
      z <- z[is.finite(z)]
      data.frame(
        data[idx[1], by, drop = FALSE],
        n = length(z),
        mean = if (length(z)) mean(z) else NA_real_,
        sd = if (length(z) > 1) stats::sd(z) else NA_real_,
        se = if (length(z) > 1) stats::sd(z) / sqrt(length(z)) else NA_real_,
        row.names = NULL
      )
    })
    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    data.frame(out[, by, drop = FALSE], metric = v,
               out[, c("n", "mean", "sd", "se"), drop = FALSE],
               stringsAsFactors = FALSE)
  })
  do.call(rbind, pieces)
}

mean_summary <- function(data, formula, vars) {
  agg <- stats::aggregate(
    stats::as.formula(paste("cbind(", paste(vars, collapse = ","), ") ~ ",
                            formula)),
    data,
    function(z) mean(z, na.rm = TRUE)
  )
  names(agg)[-(seq_along(strsplit(formula, "\\+")[[1]]))] <-
    paste0(names(agg)[-(seq_along(strsplit(formula, "\\+")[[1]]))], "_mean")
  agg
}

approx <- read_kind("approx_by_rep")
decision <- read_kind("decision_by_rep")
jaccard <- read_kind("jaccard_by_rep")
alpha_bin <- read_kind("alpha_bin_by_rep")
config <- read_kind("config_by_rep")
timing <- read_kind("timing_by_rep")

write_kind(approx, "approx_by_rep")
write_kind(decision, "decision_by_rep")
write_kind(jaccard, "jaccard_by_rep")
write_kind(alpha_bin, "alpha_bin_by_rep")
write_kind(config, "config_by_rep")
write_kind(timing, "timing_by_rep")

approx_paper <- paper_summary(
  approx,
  by = c("scenario", "method"),
  vars = c("mae_pp0", "rmse_pp0", "q95_abs_pp0", "max_abs_pp0",
           "spearman_score", "pearson_score")
)
decision_paper <- paper_summary(
  decision,
  by = c("scenario", "method", "fdr_target"),
  vars = c("n_selected", "empirical_fdr", "power", "auc", "avg_precision")
)
jaccard_paper <- paper_summary(
  jaccard,
  by = c("scenario", "method", "fdr_target"),
  vars = c("jaccard_vs_exact", "exact_only", "method_only")
)
alpha_paper <- paper_summary(
  alpha_bin,
  by = c("scenario", "method", "alpha_bin"),
  vars = c("mae_pp0", "q95_abs_pp0")
)

write_kind(approx_paper, "approx_paper_summary")
write_kind(decision_paper, "decision_paper_summary")
write_kind(jaccard_paper, "jaccard_paper_summary")
write_kind(alpha_paper, "alpha_bin_paper_summary")

approx_mean <- stats::aggregate(
  cbind(mae_pp0, rmse_pp0, q95_abs_pp0, max_abs_pp0,
        spearman_score, pearson_score) ~ scenario + method,
  approx,
  function(z) mean(z, na.rm = TRUE)
)
decision_mean <- stats::aggregate(
  cbind(n_selected, empirical_fdr, power, auc, avg_precision) ~
    scenario + method + fdr_target,
  decision,
  function(z) mean(z, na.rm = TRUE)
)
jaccard_mean <- stats::aggregate(
  cbind(jaccard_vs_exact, exact_only, method_only) ~
    scenario + method + fdr_target,
  jaccard,
  function(z) mean(z, na.rm = TRUE)
)
alpha_mean <- stats::aggregate(
  cbind(n, mae_pp0, q95_abs_pp0) ~ scenario + method + alpha_bin,
  alpha_bin,
  function(z) mean(z, na.rm = TRUE)
)
timing_mean <- stats::aggregate(
  cbind(exact_seconds, stirling_seconds, pig_seconds) ~ scenario,
  timing,
  function(z) mean(z, na.rm = TRUE)
)

write_kind(approx_mean, "approx_summary")
write_kind(decision_mean, "decision_summary")
write_kind(jaccard_mean, "jaccard_summary")
write_kind(alpha_mean, "alpha_bin_summary")
write_kind(timing_mean, "timing_summary")

fmt <- function(x, digits = 4) formatC(x, format = "f", digits = digits)

approx_key <- approx_paper[approx_paper$metric %in%
                             c("mae_pp0", "q95_abs_pp0", "spearman_score"),
                           , drop = FALSE]
approx_key$value <- paste0(fmt(approx_key$mean), " (", fmt(approx_key$se), ")")
approx_wide <- reshape(
  approx_key[, c("scenario", "method", "metric", "value")],
  idvar = c("scenario", "method"),
  timevar = "metric",
  direction = "wide"
)
names(approx_wide) <- sub("^value\\.", "", names(approx_wide))

decision05 <- decision_paper[
  abs(decision_paper$fdr_target - 0.05) < 1e-12 &
    decision_paper$metric %in% c("n_selected", "empirical_fdr", "power", "auc"),
  , drop = FALSE
]
decision05$value <- paste0(fmt(decision05$mean), " (", fmt(decision05$se), ")")
decision_wide <- reshape(
  decision05[, c("scenario", "method", "metric", "value")],
  idvar = c("scenario", "method"),
  timevar = "metric",
  direction = "wide"
)
names(decision_wide) <- sub("^value\\.", "", names(decision_wide))

manuscript_table <- merge(approx_wide, decision_wide,
                          by = c("scenario", "method"), all = TRUE)
write_kind(manuscript_table, "manuscript_table_fdr005")

pdf(file.path(result_dir, paste0(combined_prefix, "_summary_plots.pdf")),
    width = 10, height = 7)
op <- par(no.readonly = TRUE)
on.exit(par(op), add = TRUE)
par(mfrow = c(2, 2), mar = c(8, 4, 3, 1))

mat_mae <- xtabs(mae_pp0 ~ method + scenario, data = approx_mean)
barplot(mat_mae, beside = TRUE, las = 2, ylab = "mean abs error in pp0",
        main = "Posterior null probability error vs exact",
        col = c("steelblue", "gray55"))
legend("topright", legend = rownames(mat_mae),
       fill = c("steelblue", "gray55"), bty = "n")

mat_spear <- xtabs(spearman_score ~ method + scenario, data = approx_mean)
barplot(mat_spear, beside = TRUE, las = 2, ylab = "Spearman score correlation",
        main = "DE score rank agreement vs exact",
        col = c("steelblue", "gray55"), ylim = c(0.85, 1.0))
legend("bottomleft", legend = rownames(mat_spear),
       fill = c("steelblue", "gray55"), bty = "n")

dec05 <- decision_mean[abs(decision_mean$fdr_target - 0.05) < 1e-12, ]
mat_fdr <- xtabs(empirical_fdr ~ method + scenario, data = dec05)
barplot(mat_fdr, beside = TRUE, las = 2, ylab = "empirical FDR",
        main = "Empirical FDR at posterior FDR 0.05",
        col = c("black", "steelblue", "gray55"))
abline(h = 0.05, col = "red", lty = 2)
legend("topright", legend = rownames(mat_fdr),
       fill = c("black", "steelblue", "gray55"), bty = "n")

small_alpha <- alpha_mean[alpha_mean$alpha_bin == "<1", ]
mat_alpha <- xtabs(mae_pp0 ~ method + scenario, data = small_alpha)
barplot(mat_alpha, beside = TRUE, las = 2,
        ylab = "mean abs error in pp0, alpha < 1",
        main = "Small-alpha error vs exact",
        col = c("steelblue", "gray55"))
legend("topright", legend = rownames(mat_alpha),
       fill = c("steelblue", "gray55"), bty = "n")
dev.off()

cat("Combined approximation summary:\n")
print(approx_mean)
cat("\nCombined decision summary at FDR 0.05:\n")
print(decision_mean[abs(decision_mean$fdr_target - 0.05) < 1e-12, ])
cat("\nManuscript table saved to:\n")
cat(file.path(result_dir,
              paste0(combined_prefix, "_manuscript_table_fdr005.csv")),
    "\n")
