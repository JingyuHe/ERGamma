#!/usr/bin/env Rscript
# =============================================================================
# merge_table1.R - aggregate per-rep RDS files into the final summary CSV
# that we compare against Rossell 2009 Table 1.
# =============================================================================

args <- commandArgs(FALSE)
this_file <- sub("^--file=", "", args[grepl("^--file=", args)][1])
src_dir <- dirname(normalizePath(this_file, mustWork = TRUE))
source(file.path(src_dir, "benchmark_lib.R"))
bench <- setup_benchmark()

data_choice <- Sys.getenv("ARMSTRONG_DATA", unset = "schliep_filtered")
transform_choice <- Sys.getenv(
  "ARMSTRONG_TRANSFORM",
  unset = if (data_choice == "schliep_filtered") "log" else "log2_floor20"
)
equalcv <- Sys.getenv("ARMSTRONG_EQUALCV", unset = "1") != "0"
exclude_u95av2 <- Sys.getenv("ARMSTRONG_EXCLUDE_U95AV2", unset = "1") != "0"

suffix <- paste(data_choice, transform_choice,
                if (equalcv) "equalcv" else "freecv",
                if (exclude_u95av2) "exclu95av2" else "allmll",
                sep = "_")
chunk_dir <- result_file(file.path("table1_chunks", suffix))

rep_files <- list.files(chunk_dir, pattern = "^rep_[0-9]+\\.rds$",
                        full.names = TRUE)
if (length(rep_files) == 0) stop("No rep files in ", chunk_dir)
reps <- do.call(rbind, lapply(rep_files, readRDS))
full <- readRDS(file.path(chunk_dir, "full.rds"))

cat(sprintf("Merged %d rep files (%d rows). Reps observed: %s\n",
            length(rep_files), nrow(reps),
            paste(sort(unique(reps$rep)), collapse = ",")))

summary_rows <- aggregate(
  cbind(n_de, reproducibility, elapsed) ~ method + n_per_group,
  data = reps, FUN = function(v) mean(v, na.rm = TRUE))
reps$has_error <- !is.na(reps$error)
errors <- aggregate(has_error ~ method + n_per_group, data = reps, FUN = sum)
names(errors)[3] <- "n_errors"
summary_rows <- merge(summary_rows, errors,
                      by = c("method", "n_per_group"), all.x = TRUE)
summary_rows$n_success <- max(reps$rep) - summary_rows$n_errors

full_summary <- do.call(rbind, lapply(names(full$lists), function(m) {
  data.frame(method = m, n_per_group = "All data",
             n_de = length(full$lists[[m]]),
             reproducibility = NA_real_,
             elapsed = full$elapsed[[m]],
             n_errors = as.integer(!is.na(full$error[[m]])),
             n_success = as.integer(is.na(full$error[[m]])),
             stringsAsFactors = FALSE)
}))
summary_rows <- rbind(summary_rows, full_summary[, names(summary_rows)])
summary_rows$n_per_group <- factor(summary_rows$n_per_group,
                                   levels = c("5", "10", "15", "All data"))
summary_rows <- summary_rows[order(summary_rows$method,
                                   summary_rows$n_per_group), ]

ross <- data.frame(
  method      = rep(c("GaGa", "MiGaGa2", "limma_BH"), each = 4),
  n_per_group = factor(rep(c("5", "10", "15", "All data"), 3),
                       levels = c("5", "10", "15", "All data")),
  paper_n_de  = c(58.5, 431, 784, 991,
                  61.5, 445, 815, 1040,
                  21.5, 181.5, 543, 972),
  paper_repro = c(0.856, 0.893, 0.889, NA,
                  0.860, 0.893, 0.890, NA,
                  0.947, 0.957, 0.946, NA)
)
out <- merge(summary_rows, ross, by = c("method", "n_per_group"), all.x = TRUE)
out$ratio_n_de <- round(out$n_de / out$paper_n_de, 3)
out <- out[order(out$method, out$n_per_group), ]

write.csv(reps,
  result_file(paste0("armstrong_table1_replicates_", suffix, ".csv")),
  row.names = FALSE)
write.csv(out,
  result_file(paste0("armstrong_table1_summary_", suffix, ".csv")),
  row.names = FALSE)
saveRDS(full$lists,
  result_file(paste0("armstrong_table1_full_lists_", suffix, ".rds")))

cat("\n=== Table 1 reproduction (Rossell 2009 reference merged) ===\n")
print(out, row.names = FALSE)
cat("\nOutputs saved with suffix:", suffix, "\n")
