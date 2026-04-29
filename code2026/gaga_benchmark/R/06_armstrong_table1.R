args <- commandArgs(FALSE)
this_file <- sub("^--file=", "", args[grepl("^--file=", args)][1])
source(file.path(dirname(normalizePath(this_file, mustWork = TRUE)), "benchmark_lib.R"))

bench <- setup_benchmark()
require_pkgs(c("gaga", "limma"))

load_armstrong_orange <- function(data_file, exclude_last_mll = TRUE) {
  raw <- read.delim(data_file, check.names = FALSE, stringsAsFactors = FALSE)
  groups <- raw$class[-c(1, 2)]
  xdf <- raw[-c(1, 2), -1, drop = FALSE]
  x <- t(apply(xdf, 2, as.numeric))
  rownames(x) <- colnames(raw)[-1]
  sample_ids <- paste0(groups, "_", ave(seq_along(groups), groups, FUN = seq_along))
  colnames(x) <- sample_ids
  keep <- groups %in% c("ALL", "MLL")
  if (exclude_last_mll) {
    mll_idx <- which(groups == "MLL")
    keep[tail(mll_idx, 2)] <- FALSE
  }
  list(x = x[, keep, drop = FALSE], groups = factor(groups[keep], levels = c("ALL", "MLL")))
}

load_armstrong_schliep <- function(data_file, exclude_last_mll = TRUE) {
  raw <- read.delim(data_file, check.names = FALSE, stringsAsFactors = FALSE)
  groups <- as.character(unlist(raw[1, -1]))
  x <- as.matrix(data.frame(lapply(raw[-1, -1, drop = FALSE], as.numeric), check.names = FALSE))
  rownames(x) <- raw[-1, 1]
  sample_ids <- paste0(groups, "_", ave(seq_along(groups), groups, FUN = seq_along))
  colnames(x) <- sample_ids
  keep <- groups %in% c("ALL", "MLL")
  if (exclude_last_mll) {
    mll_idx <- which(groups == "MLL")
    keep[tail(mll_idx, 2)] <- FALSE
  }
  list(x = x[, keep, drop = FALSE], groups = factor(groups[keep], levels = c("ALL", "MLL")))
}

transform_expr <- function(x, transform) {
  if (transform == "log2_clip") {
    return(log2(pmax(x, 1) + 1))
  }
  if (transform == "clip") {
    return(pmax(x, 1))
  }
  if (transform == "shift") {
    return(x - min(x, na.rm = TRUE) + 1)
  }
  if (transform == "none") {
    if (min(x, na.rm = TRUE) <= 0) {
      stop("Expression matrix has non-positive values; choose a positive transform.")
    }
    return(x)
  }
  stop("Unknown ARMSTRONG_TRANSFORM: ", transform)
}

select_genes <- function(x, n_genes) {
  if (n_genes <= 0 || n_genes >= nrow(x)) {
    return(x)
  }
  vars <- matrixStats::rowVars(x)
  keep <- order(vars, decreasing = TRUE)[seq_len(n_genes)]
  x[keep, , drop = FALSE]
}

fit_gaga_calls <- function(x, groups, nclust) {
  patterns <- matrix(c(0, 0, 0, 1), 2, 2)
  colnames(patterns) <- levels(factor(groups))
  preferred <- Sys.getenv("ARMSTRONG_GAGA_METHOD", unset = "quickEM")
  fit_methods <- unique(c(preferred, "EM"))
  errors <- character()
  for (fit_method in fit_methods) {
    stage <- "fitGG"
    out <- tryCatch({
      fit <- gaga::fitGG(
        x,
        groups,
        patterns = patterns,
        equalcv = TRUE,
        nclust = nclust,
        method = fit_method,
        trace = FALSE
      )
      stage <- "parest"
      fit <- gaga::parest(fit, x = x, groups = groups, alpha = 0.05)
      stage <- "findgenes"
      genes <- gaga::findgenes(fit, x, groups, fdrmax = 0.05, parametric = TRUE)
      rownames(x)[as.logical(genes$d != 0)]
    }, error = function(e) {
      errors <<- c(errors, paste0("method=", fit_method, ", stage=", stage, ": ", conditionMessage(e)))
      NULL
    })
    if (!is.null(out)) {
      return(out)
    }
  }
  stop("GaGa fit failed after quickEM/EM fallback: ", paste(errors, collapse = " | "))
}

fit_limma_calls <- function(x, groups) {
  lim <- limma_or_ttest(x, groups, fdr = 0.05)
  rownames(x)[lim$calls]
}

fit_method_calls <- function(method, x, groups) {
  if (method == "GaGa") {
    return(fit_gaga_calls(x, groups, nclust = 1))
  }
  if (method == "MiGaGa2") {
    return(fit_gaga_calls(x, groups, nclust = 2))
  }
  if (method == "limma_BH") {
    return(fit_limma_calls(x, groups))
  }
  stop("Unknown method: ", method)
}

safe_fit_method_calls <- function(method, x, groups) {
  out <- tryCatch({
    list(calls = fit_method_calls(method, x, groups), error = NA_character_)
  }, error = function(e) {
    list(calls = character(0), error = conditionMessage(e))
  })
  out
}

run_table1 <- function(x, groups, methods, reps, seed) {
  set.seed(seed)
  all_idx <- which(groups == "ALL")
  mll_idx <- which(groups == "MLL")
  if (length(all_idx) < 15 || length(mll_idx) < 15) {
    stop("Need at least 15 ALL and 15 MLL samples for Table 1 subsets.")
  }
  sample_sizes <- c(5, 10, 15)
  full_lists <- list()
  full_rows <- list()
  for (method in methods) {
    elapsed <- system.time({
      fit_out <- safe_fit_method_calls(method, x, groups)
      full_lists[[method]] <- fit_out$calls
    })[["elapsed"]]
    full_rows[[method]] <- data.frame(
      method = method,
      n_per_group = "All data",
      rep = NA_integer_,
      n_de = if (is.na(fit_out$error)) length(full_lists[[method]]) else NA_integer_,
      reproducibility = NA_real_,
      elapsed = elapsed,
      error = fit_out$error
    )
  }

  rep_rows <- list()
  row_id <- 1
  for (rep_id in seq_len(reps)) {
    all_perm <- sample(all_idx)
    mll_perm <- sample(mll_idx)
    for (ss in sample_sizes) {
      cols <- c(all_perm[seq_len(ss)], mll_perm[seq_len(ss)])
      for (method in methods) {
        elapsed <- system.time({
          fit_out <- safe_fit_method_calls(method, x[, cols, drop = FALSE], groups[cols])
          calls <- fit_out$calls
        })[["elapsed"]]
        rep_rows[[row_id]] <- data.frame(
          method = method,
          n_per_group = as.character(ss),
          rep = rep_id,
          n_de = if (is.na(fit_out$error)) length(calls) else NA_integer_,
          reproducibility = if (!is.na(fit_out$error) || length(calls) == 0) {
            NA_real_
          } else {
            mean(calls %in% full_lists[[method]])
          },
          elapsed = elapsed,
          error = fit_out$error
        )
        row_id <- row_id + 1
      }
    }
  }
  list(full = do.call(rbind, full_rows), reps = do.call(rbind, rep_rows),
       full_lists = full_lists)
}

data_choice <- Sys.getenv("ARMSTRONG_DATA", unset = "orange_full")
transform_choice <- Sys.getenv("ARMSTRONG_TRANSFORM", unset = "log2_clip")
exclude_last_mll <- Sys.getenv("ARMSTRONG_EXCLUDE_LAST_MLL", unset = "1") != "0"
reps <- env_int("ARMSTRONG_REPS", 20)
n_genes <- env_int("ARMSTRONG_N_GENES", 0)
seed <- env_int("ARMSTRONG_SEED", 101)
methods <- strsplit(Sys.getenv("ARMSTRONG_METHODS", unset = "GaGa,MiGaGa2,limma_BH"), ",")[[1]]
methods <- trimws(methods)

if (data_choice == "orange_full") {
  data_file <- file.path(bench, "data", "MLL_orange_12533.tab")
  dat <- load_armstrong_orange(data_file, exclude_last_mll = exclude_last_mll)
} else if (data_choice == "schliep_filtered") {
  data_file <- file.path(bench, "data", "armstrong-2002-v2_database.txt")
  dat <- load_armstrong_schliep(data_file, exclude_last_mll = exclude_last_mll)
} else {
  stop("Unknown ARMSTRONG_DATA: ", data_choice)
}

x_raw <- dat$x
x <- transform_expr(x_raw, transform_choice)
if (!requireNamespace("matrixStats", quietly = TRUE)) {
  if (n_genes > 0 && n_genes < nrow(x)) {
    vars <- apply(x, 1, stats::var)
    keep <- order(vars, decreasing = TRUE)[seq_len(n_genes)]
    x <- x[keep, , drop = FALSE]
  }
} else {
  x <- select_genes(x, n_genes)
}
groups <- dat$groups

data_summary <- data.frame(
  data = data_choice,
  source_file = data_file,
  transform = transform_choice,
  exclude_last_mll = exclude_last_mll,
  n_genes = nrow(x),
  n_samples = ncol(x),
  n_all = sum(groups == "ALL"),
  n_mll = sum(groups == "MLL"),
  min_raw = min(x_raw),
  max_raw = max(x_raw),
  min_used = min(x),
  max_used = max(x),
  reps = reps,
  methods = paste(methods, collapse = ",")
)
write.csv(data_summary, result_file("armstrong_table1_data_summary.csv"), row.names = FALSE)
print(data_summary)

res <- run_table1(x, groups, methods, reps, seed)
rep_rows <- res$reps
full_rows <- res$full

summary_rows <- aggregate(
  cbind(n_de, reproducibility, elapsed) ~ method + n_per_group,
  data = rep_rows,
  FUN = function(v) mean(v, na.rm = TRUE)
)
rep_rows$has_error <- !is.na(rep_rows$error)
errors <- aggregate(
  has_error ~ method + n_per_group,
  data = rep_rows,
  FUN = sum
)
names(errors)[3] <- "n_errors"
summary_rows <- merge(summary_rows, errors, by = c("method", "n_per_group"), all.x = TRUE)
summary_rows$n_success <- reps - summary_rows$n_errors
full_summary <- full_rows[, c("method", "n_per_group", "n_de", "reproducibility", "elapsed")]
full_summary$n_errors <- as.integer(!is.na(full_rows$error))
full_summary$n_success <- 1L - full_summary$n_errors
summary_rows <- rbind(summary_rows, full_summary[, names(summary_rows)])
summary_rows$n_per_group <- factor(summary_rows$n_per_group, levels = c("5", "10", "15", "All data"))
summary_rows <- summary_rows[order(summary_rows$method, summary_rows$n_per_group), ]

suffix <- paste(data_choice, transform_choice, if (exclude_last_mll) "exclude2mll" else "allmll", sep = "_")
write.csv(rep_rows, result_file("armstrong_table1_replicates.csv"), row.names = FALSE)
write.csv(summary_rows, result_file("armstrong_table1_summary.csv"), row.names = FALSE)
saveRDS(res$full_lists, result_file("armstrong_table1_full_lists.rds"))
write.csv(data_summary, result_file(paste0("armstrong_table1_data_summary_", suffix, ".csv")), row.names = FALSE)
write.csv(rep_rows, result_file(paste0("armstrong_table1_replicates_", suffix, ".csv")), row.names = FALSE)
write.csv(summary_rows, result_file(paste0("armstrong_table1_summary_", suffix, ".csv")), row.names = FALSE)
saveRDS(res$full_lists, result_file(paste0("armstrong_table1_full_lists_", suffix, ".rds")))
print(summary_rows)
