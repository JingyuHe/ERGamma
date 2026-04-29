#!/usr/bin/env Rscript
# merge_sim62.R - aggregate per-rep RDS into Section 6.2 Table 2 / ROC summary.
args <- commandArgs(FALSE)
this_file <- sub("^--file=", "", args[grepl("^--file=", args)][1])
src_dir <- dirname(normalizePath(this_file, mustWork = TRUE))
source(file.path(src_dir, "benchmark_lib.R"))
bench <- setup_benchmark()

data_choice <- Sys.getenv("ARM62_DATA", unset = "schliep_filtered")
transform_choice <- Sys.getenv(
  "ARM62_TRANSFORM",
  unset = if (data_choice == "schliep_filtered") "log" else "log2_floor20")
n_reps <- env_int("ARM62_REPS", 200)
equalcv <- Sys.getenv("ARM62_EQUALCV", unset = "1") != "0"

suffix <- paste(data_choice, transform_choice, paste0("reps", n_reps),
                if (equalcv) "equalcv" else "freecv", sep = "_")
chunk_dir <- result_file(file.path("sim62_chunks", suffix))

rows_all <- list(); roc_all <- list()
for (st in c("parametric", "nonparametric")) {
  for (f in list.files(file.path(chunk_dir, st),
                       pattern = "^rep_[0-9]+\\.rds$", full.names = TRUE)) {
    obj <- readRDS(f)
    rows_all[[length(rows_all) + 1]] <- obj$rows
    if (!is.null(obj$roc))
      roc_all[[length(roc_all) + 1]] <- obj$roc
  }
}
rows <- do.call(rbind, rows_all)
roc  <- if (length(roc_all)) do.call(rbind, roc_all) else NULL
cat(sprintf("Merged %d rep snapshots, %d row records.\n",
            length(rows_all), nrow(rows)))

rows$has_error <- !is.na(rows$error)
sm <- aggregate(cbind(n_de, true_de, fdr, power, elapsed) ~
                  simulation + method + n_per_group,
                data = rows, FUN = function(v) mean(v, na.rm = TRUE))
errs <- aggregate(has_error ~ simulation + method + n_per_group,
                  data = rows, FUN = sum)
names(errs)[4] <- "n_errors"
sm <- merge(sm, errs, by = c("simulation", "method", "n_per_group"),
            all.x = TRUE)

paper_table2 <- data.frame(
  simulation = rep(c("parametric", "nonparametric"), each = 16),
  method = rep(rep(c("GaGa", "MiGaGa2", "Ga", "limma_BH"), each = 4), 2),
  n_per_group = rep(c(5, 10, 15, 20), 8),
  paper_fdr = c(
    0.011, 0.000, 0.007, 0.002,
    0.011, 0.000, 0.008, 0.002,
    0.159, 0.133, 0.117, 0.105,
    0.012, 0.035, 0.034, 0.036,
    0.043, 0.066, 0.067, 0.065,
    0.047, 0.066, 0.068, 0.068,
    0.342, 0.289, 0.254, 0.239,
    0.047, 0.021, 0.019, 0.024),
  paper_power = c(
    0.066, 0.322, 0.512, 0.608,
    0.066, 0.328, 0.520, 0.615,
    0.434, 0.587, 0.667, 0.712,
    0.063, 0.288, 0.487, 0.580,
    0.054, 0.319, 0.529, 0.660,
    0.057, 0.327, 0.541, 0.671,
    0.397, 0.567, 0.666, 0.740,
    0.029, 0.168, 0.321, 0.449))
sm <- merge(sm, paper_table2,
            by = c("simulation", "method", "n_per_group"), all.x = TRUE)
sm <- sm[order(sm$simulation, sm$method, sm$n_per_group), ]
sm$fdr <- round(sm$fdr, 4); sm$power <- round(sm$power, 4)

write.csv(rows,
  result_file(paste0("armstrong62_simulation_rows_", suffix, ".csv")),
  row.names = FALSE)
write.csv(sm,
  result_file(paste0("armstrong62_table2_summary_", suffix, ".csv")),
  row.names = FALSE)

if (!is.null(roc)) {
  roc_avg <- aggregate(power ~ method + fdr_grid, data = roc, FUN = mean)
  roc_avg <- roc_avg[order(roc_avg$method, roc_avg$fdr_grid), ]
  auc_grid <- function(x, y) sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
  roc_summary <- do.call(rbind, lapply(split(roc_avg, roc_avg$method),
    function(z) data.frame(method = z$method[1],
                           auc = auc_grid(z$fdr_grid, z$power),
                           stringsAsFactors = FALSE)))
  rownames(roc_summary) <- NULL
  write.csv(roc_avg,
    result_file(paste0("armstrong62_roc_curves_", suffix, ".csv")),
    row.names = FALSE)
  write.csv(roc_summary,
    result_file(paste0("armstrong62_roc_summary_", suffix, ".csv")),
    row.names = FALSE)

  plot_roc <- function(file, dev = c("pdf", "png")) {
    dev <- match.arg(dev)
    if (dev == "pdf") grDevices::pdf(file, width = 6.2, height = 5.2)
    else grDevices::png(file, width = 1400, height = 1100, res = 200)
    cols <- c(GaGa = "#1B9E77", MiGaGa2 = "#D95F02",
              Ga = "#666666", limma_BH = "#7570B3")
    plot(NA, xlim = c(0, 1), ylim = c(0, 1),
         xlab = "Average FDR", ylab = "Average power",
         main = "Armstrong nonparametric simulations")
    for (m in names(cols)) {
      z <- roc_avg[roc_avg$method == m, ]
      if (nrow(z) > 0) lines(z$fdr_grid, z$power,
                             col = cols[[m]], lwd = 2)
    }
    leg <- merge(data.frame(method = names(cols)), roc_summary,
                 by = "method", all.x = TRUE, sort = FALSE)
    legend("bottomright",
           legend = paste0(leg$method, " AUC=",
                           sprintf("%.3f", leg$auc)),
           col = cols[leg$method], lwd = 2, bty = "n")
    grDevices::dev.off()
  }
  plot_roc(result_file(paste0("armstrong62_roc_curves_", suffix, ".pdf")), "pdf")
  plot_roc(result_file(paste0("armstrong62_roc_curves_", suffix, ".png")), "png")

  cat("\n=== ROC AUC summary ===\n")
  print(roc_summary, row.names = FALSE)
}

cat("\n=== Section 6.2 Table 2 reproduction (paper merged) ===\n")
print(sm, row.names = FALSE)
cat("\nOutputs suffix:", suffix, "\n")
