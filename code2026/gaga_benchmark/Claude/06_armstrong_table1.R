#!/usr/bin/env Rscript
# =============================================================================
# 06_armstrong_table1.R  (Claude/ debugged version)
#
# Reproduces Rossell (2009) Table 1: Armstrong ALL vs MLL, 24 ALL + 18 MLL,
# 20 random subsets at n = 5, 10, 15 per group.
#
# Fixes vs gaga_benchmark/R/06_armstrong_table1.R:
#   * limma is given log-scale data ONCE. The legacy script log2-transformed
#     the matrix in transform_expr() and then limma_or_ttest applied log2()
#     inside, double-logging.
#   * U95Av2 MLL arrays are dropped using CDF info when the data file carries
#     it; otherwise the script falls back to "tail-2" but warns explicitly so
#     the user can verify the file ordering.
#   * equalcv = TRUE is used (paper setting) and exposed via
#     ARMSTRONG_EQUALCV (defaults to TRUE).
#   * Adds Ga (EBarrays GG) as an opt-in baseline. Ga is fit on a positive-
#     scale version of the matrix (2^x for log2 transforms, exp(x) for
#     natural-log transforms), so the GG model sees data on its native scale.
# =============================================================================

args <- commandArgs(FALSE)
this_file <- sub("^--file=", "", args[grepl("^--file=", args)][1])
source(file.path(dirname(normalizePath(this_file, mustWork = TRUE)),
                 "benchmark_lib.R"))

bench <- setup_benchmark()
require_pkgs(c("gaga", "limma"))

# ---------------------------------------------------------------------------
# Loaders
# ---------------------------------------------------------------------------
load_armstrong_orange <- function(data_file) {
  raw <- read.delim(data_file, check.names = FALSE, stringsAsFactors = FALSE)
  groups <- raw$class[-c(1, 2)]
  xdf <- raw[-c(1, 2), -1, drop = FALSE]
  x <- t(apply(xdf, 2, as.numeric))
  rownames(x) <- colnames(raw)[-1]
  sample_ids <- paste0(groups, "_",
                       ave(seq_along(groups), groups, FUN = seq_along))
  colnames(x) <- sample_ids
  list(x = x, groups = groups, cdf = NULL)
}

load_armstrong_schliep <- function(data_file) {
  raw <- read.delim(data_file, check.names = FALSE, stringsAsFactors = FALSE)
  groups <- as.character(unlist(raw[1, -1]))
  x <- as.matrix(data.frame(lapply(raw[-1, -1, drop = FALSE], as.numeric),
                            check.names = FALSE))
  rownames(x) <- raw[-1, 1]
  sample_ids <- paste0(groups, "_",
                       ave(seq_along(groups), groups, FUN = seq_along))
  colnames(x) <- sample_ids
  list(x = x, groups = groups, cdf = NULL)
}

load_armstrong_gcrma <- function(rds_file) {
  obj <- readRDS(rds_file)
  list(x = obj$expr, groups = as.character(obj$groups),
       cdf = obj$sample_info$cdf_name)
}

restrict_to_all_mll <- function(dat, exclude_u95av2 = TRUE) {
  keep <- dat$groups %in% c("ALL", "MLL")
  if (exclude_u95av2) {
    drop_idx <- cdf_safe_drop_u95av2(dat$groups, dat$cdf)
    if (length(drop_idx) > 0) {
      keep[drop_idx] <- FALSE
    }
  }
  list(
    x = dat$x[, keep, drop = FALSE],
    groups = factor(dat$groups[keep], levels = c("ALL", "MLL"))
  )
}

# ---------------------------------------------------------------------------
# Transforms - all return a "modeling-scale" (log) matrix and a "raw-scale"
# version for Ga/GG to use.
# ---------------------------------------------------------------------------
apply_transform <- function(x_raw, transform) {
  if (transform == "log2_clip") {
    x_log <- log2(pmax(x_raw, 1) + 1)
    x_pos <- pmax(x_raw, 1)              # for Ga/GG
    base <- 2
  } else if (transform == "log2_floor20") {
    x_log <- log2(pmax(x_raw, 20))
    x_pos <- pmax(x_raw, 20)
    base <- 2
  } else if (transform == "log") {
    x_log <- log(pmax(x_raw, 1))
    x_pos <- pmax(x_raw, 1)
    base <- exp(1)
  } else if (transform == "none") {
    if (min(x_raw, na.rm = TRUE) <= 0) {
      stop("Transform=none requires positive matrix.")
    }
    x_log <- log2(x_raw)
    x_pos <- x_raw
    base <- 2
  } else {
    stop("Unknown transform: ", transform)
  }
  list(x_log = x_log, x_pos = x_pos, base = base)
}

# ---------------------------------------------------------------------------
# Methods - all consume a "fit_data" list with $x_log and $x_pos.
# ---------------------------------------------------------------------------
make_two_group_pattern <- function(groups) {
  patterns <- matrix(c(0, 0, 0, 1), 2, 2)
  colnames(patterns) <- levels(factor(groups))
  patterns
}

fit_gaga_calls <- function(fit_data, groups, nclust, equalcv,
                           preferred = "quickEM") {
  patterns <- make_two_group_pattern(groups)
  fit_methods <- unique(c(preferred, "EM"))
  errors <- character()
  for (fit_method in fit_methods) {
    out <- tryCatch({
      fit <- gaga::fitGG(fit_data$x_log, groups, patterns = patterns,
                         equalcv = equalcv, nclust = nclust,
                         method = fit_method, trace = FALSE)
      fit <- gaga::parest(fit, x = fit_data$x_log, groups = groups,
                          alpha = 0.05)
      genes <- gaga::findgenes(fit, fit_data$x_log, groups, fdrmax = 0.05,
                               parametric = TRUE)
      rownames(fit_data$x_log)[as.logical(genes$d != 0)]
    }, error = function(e) {
      errors <<- c(errors,
                   paste0("method=", fit_method, ": ", conditionMessage(e)))
      NULL
    })
    if (!is.null(out)) return(out)
  }
  stop("GaGa fit failed: ", paste(errors, collapse = " | "))
}

fit_limma_calls <- function(fit_data, groups) {
  # FIX: pass log-scale data, do NOT log again inside
  res <- limma_BH(fit_data$x_log, groups, fdr = 0.05)
  rownames(fit_data$x_log)[res$calls]
}

fit_ga_calls <- function(fit_data, groups) {
  if (!requireNamespace("EBarrays", quietly = TRUE)) {
    stop("EBarrays not available; install via BiocManager.")
  }
  res <- ebarrays_GG_raw(fit_data$x_pos, groups, fdr = 0.05)
  rownames(fit_data$x_log)[res$calls]
}

fit_method_calls <- function(method, fit_data, groups, equalcv) {
  if (method == "GaGa")    return(fit_gaga_calls(fit_data, groups, 1, equalcv))
  if (method == "MiGaGa2") return(fit_gaga_calls(fit_data, groups, 2, equalcv))
  if (method == "Ga")      return(fit_ga_calls(fit_data, groups))
  if (method == "limma_BH") return(fit_limma_calls(fit_data, groups))
  stop("Unknown method: ", method)
}

safe_fit <- function(method, fit_data, groups, equalcv) {
  tryCatch({
    list(calls = fit_method_calls(method, fit_data, groups, equalcv),
         error = NA_character_)
  }, error = function(e) {
    list(calls = character(0), error = conditionMessage(e))
  })
}

# ---------------------------------------------------------------------------
# Main loop
# ---------------------------------------------------------------------------
run_table1 <- function(x_raw, groups, transform, methods, reps, seed,
                       equalcv) {
  set.seed(seed)
  full_data <- apply_transform(x_raw, transform)
  full_data$x_log <- {
    z <- full_data$x_log
    rownames(z) <- rownames(x_raw)
    z
  }
  full_data$x_pos <- {
    z <- full_data$x_pos
    rownames(z) <- rownames(x_raw)
    z
  }

  full_lists <- list()
  full_rows <- list()
  for (method in methods) {
    elapsed <- system.time({
      out <- safe_fit(method, full_data, groups, equalcv)
      full_lists[[method]] <- out$calls
    })[["elapsed"]]
    full_rows[[method]] <- data.frame(
      method = method, n_per_group = "All data", rep = NA_integer_,
      n_de = if (is.na(out$error)) length(out$calls) else NA_integer_,
      reproducibility = NA_real_, elapsed = elapsed,
      error = out$error, stringsAsFactors = FALSE
    )
  }

  all_idx <- which(groups == "ALL")
  mll_idx <- which(groups == "MLL")
  if (length(all_idx) < 15 || length(mll_idx) < 15) {
    stop("Need at least 15 ALL and 15 MLL samples for Table 1 subsets.")
  }
  rep_rows <- list(); ix <- 1
  for (rep_id in seq_len(reps)) {
    ap <- sample(all_idx); mp <- sample(mll_idx)
    for (ss in c(5, 10, 15)) {
      cols <- c(ap[seq_len(ss)], mp[seq_len(ss)])
      sub_log <- full_data$x_log[, cols, drop = FALSE]
      sub_pos <- full_data$x_pos[, cols, drop = FALSE]
      sub_grp <- factor(rep(c("ALL", "MLL"), each = ss),
                        levels = c("ALL", "MLL"))
      sub_data <- list(x_log = sub_log, x_pos = sub_pos)
      for (method in methods) {
        elapsed <- system.time({
          out <- safe_fit(method, sub_data, sub_grp, equalcv)
        })[["elapsed"]]
        repro <- if (is.na(out$error) && length(out$calls) > 0) {
          mean(out$calls %in% full_lists[[method]])
        } else {
          NA_real_
        }
        rep_rows[[ix]] <- data.frame(
          method = method, n_per_group = as.character(ss), rep = rep_id,
          n_de = if (is.na(out$error)) length(out$calls) else NA_integer_,
          reproducibility = repro, elapsed = elapsed, error = out$error,
          stringsAsFactors = FALSE
        )
        ix <- ix + 1
      }
    }
    cat(sprintf("  rep %d/%d done\n", rep_id, reps))
  }

  list(full = do.call(rbind, full_rows),
       reps = do.call(rbind, rep_rows),
       full_lists = full_lists)
}

# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------
data_choice <- Sys.getenv("ARMSTRONG_DATA", unset = "orange_full")
transform_choice <- Sys.getenv("ARMSTRONG_TRANSFORM", unset = "log2_floor20")
exclude_u95av2 <- Sys.getenv("ARMSTRONG_EXCLUDE_U95AV2", unset = "1") != "0"
reps <- env_int("ARMSTRONG_REPS", 20)
seed <- env_int("ARMSTRONG_SEED", 101)
equalcv <- Sys.getenv("ARMSTRONG_EQUALCV", unset = "1") != "0"
methods <- strsplit(Sys.getenv("ARMSTRONG_METHODS",
                               unset = "GaGa,MiGaGa2,limma_BH"), ",")[[1]]
methods <- trimws(methods)

if (data_choice == "orange_full") {
  data_file <- file.path(bench, "data", "MLL_orange_12533.tab")
  dat <- load_armstrong_orange(data_file)
} else if (data_choice == "schliep_filtered") {
  data_file <- file.path(bench, "data", "armstrong-2002-v2_database.txt")
  dat <- load_armstrong_schliep(data_file)
} else if (data_choice == "gcrma_rds") {
  data_file <- Sys.getenv("ARMSTRONG_GCRMA_RDS",
                          unset = file.path(bench, "results",
                                            "armstrong_gcrma_expression.rds"))
  dat <- load_armstrong_gcrma(data_file)
} else {
  stop("Unknown ARMSTRONG_DATA: ", data_choice)
}

dat <- restrict_to_all_mll(dat, exclude_u95av2 = exclude_u95av2)
x_raw <- dat$x
groups <- dat$groups

data_summary <- data.frame(
  data = data_choice,
  source_file = data_file,
  transform = transform_choice,
  exclude_u95av2 = exclude_u95av2,
  equalcv = equalcv,
  n_genes = nrow(x_raw),
  n_samples = ncol(x_raw),
  n_all = sum(groups == "ALL"),
  n_mll = sum(groups == "MLL"),
  reps = reps,
  methods = paste(methods, collapse = ",")
)
print(data_summary)

res <- run_table1(x_raw, groups, transform_choice, methods, reps, seed, equalcv)

summary_rows <- aggregate(
  cbind(n_de, reproducibility, elapsed) ~ method + n_per_group,
  data = res$reps, FUN = function(v) mean(v, na.rm = TRUE)
)
res$reps$has_error <- !is.na(res$reps$error)
errors <- aggregate(has_error ~ method + n_per_group, data = res$reps,
                    FUN = sum)
names(errors)[3] <- "n_errors"
summary_rows <- merge(summary_rows, errors, by = c("method", "n_per_group"),
                      all.x = TRUE)
summary_rows$n_success <- reps - summary_rows$n_errors
full_summary <- res$full[, c("method", "n_per_group", "n_de",
                             "reproducibility", "elapsed")]
full_summary$n_errors <- as.integer(!is.na(res$full$error))
full_summary$n_success <- 1L - full_summary$n_errors
summary_rows <- rbind(summary_rows, full_summary[, names(summary_rows)])
summary_rows$n_per_group <- factor(summary_rows$n_per_group,
                                   levels = c("5", "10", "15", "All data"))
summary_rows <- summary_rows[order(summary_rows$method,
                                   summary_rows$n_per_group), ]

# Compare to Rossell 2009 Table 1
ross <- data.frame(
  method      = rep(c("GaGa", "MiGaGa2", "limma_BH"), each = 4),
  n_per_group = factor(rep(c("5", "10", "15", "All data"), 3),
                       levels = c("5", "10", "15", "All data")),
  paper_n_de  = c( 58.5, 431,   784,  991,
                   61.5, 445,   815, 1040,
                   21.5, 181.5, 543,  972),
  paper_repro = c(0.856, 0.893, 0.889, NA,
                  0.860, 0.893, 0.890, NA,
                  0.947, 0.957, 0.946, NA)
)
summary_rows <- merge(summary_rows, ross,
                      by = c("method", "n_per_group"), all.x = TRUE)
summary_rows$ratio_n_de <- round(summary_rows$n_de / summary_rows$paper_n_de,
                                 3)
summary_rows <- summary_rows[order(summary_rows$method,
                                   summary_rows$n_per_group), ]

suffix <- paste(data_choice, transform_choice,
                if (equalcv) "equalcv" else "freecv",
                if (exclude_u95av2) "exclu95av2" else "allmll",
                sep = "_")
write.csv(data_summary,
          result_file(paste0("armstrong_table1_data_summary_", suffix, ".csv")),
          row.names = FALSE)
write.csv(res$reps,
          result_file(paste0("armstrong_table1_replicates_", suffix, ".csv")),
          row.names = FALSE)
write.csv(summary_rows,
          result_file(paste0("armstrong_table1_summary_", suffix, ".csv")),
          row.names = FALSE)
saveRDS(res$full_lists,
        result_file(paste0("armstrong_table1_full_lists_", suffix, ".rds")))

cat("\n=== Table 1 reproduction (vs Rossell 2009 Table 1) ===\n")
print(summary_rows, row.names = FALSE)
cat("\nOutputs saved with suffix:", suffix, "\n")
