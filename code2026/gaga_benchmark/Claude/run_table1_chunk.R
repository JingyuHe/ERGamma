#!/usr/bin/env Rscript
# =============================================================================
# run_table1_chunk.R - resumable chunked runner for Table 1.
#
# Each invocation processes pending reps until a soft time budget is hit, then
# exits. Each completed rep is saved to its own RDS file under
# Claude/results/table1_chunks/<suffix>/rep_<id>.rds. The full-data fits are
# done once per (data, transform, equalcv) and saved alongside.
#
# Re-running the script picks up from the next missing rep. After all reps
# are done, the merge script merge_table1.R aggregates them into the
# usual armstrong_table1_summary_*.csv output that we compare against
# Rossell 2009 Table 1.
# =============================================================================

args <- commandArgs(FALSE)
this_file <- sub("^--file=", "", args[grepl("^--file=", args)][1])
src_dir <- dirname(normalizePath(this_file, mustWork = TRUE))
source(file.path(src_dir, "benchmark_lib.R"))

bench <- setup_benchmark()
require_pkgs(c("gaga", "limma"))
if ("EBarrays" %in% strsplit(Sys.getenv("ARMSTRONG_METHODS",
       unset = "GaGa,MiGaGa2,limma_BH"), ",")[[1]]) {
  require_pkgs("EBarrays")
}

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
data_choice <- Sys.getenv("ARMSTRONG_DATA", unset = "schliep_filtered")
transform_choice <- Sys.getenv(
  "ARMSTRONG_TRANSFORM",
  unset = if (data_choice == "schliep_filtered") "log" else "log2_floor20"
)
equalcv <- Sys.getenv("ARMSTRONG_EQUALCV", unset = "1") != "0"
total_reps <- env_int("ARMSTRONG_REPS", 20)
base_seed <- env_int("ARMSTRONG_SEED", 101)
exclude_u95av2 <- Sys.getenv("ARMSTRONG_EXCLUDE_U95AV2", unset = "1") != "0"
methods <- trimws(strsplit(
  Sys.getenv("ARMSTRONG_METHODS", unset = "GaGa,MiGaGa2,limma_BH"), ",")[[1]])
budget_sec <- as.numeric(Sys.getenv("CHUNK_BUDGET_SEC", unset = "30"))

# ---------------------------------------------------------------------------
# Loaders / transforms (subset of 06)
# ---------------------------------------------------------------------------
load_armstrong_orange <- function(data_file) {
  raw <- read.delim(data_file, check.names = FALSE, stringsAsFactors = FALSE)
  groups <- raw$class[-c(1, 2)]
  xdf <- raw[-c(1, 2), -1, drop = FALSE]
  x <- t(apply(xdf, 2, as.numeric))
  rownames(x) <- colnames(raw)[-1]
  list(x = x, groups = groups, cdf = NULL)
}
load_armstrong_schliep <- function(data_file) {
  raw <- read.delim(data_file, check.names = FALSE, stringsAsFactors = FALSE)
  groups <- as.character(unlist(raw[1, -1]))
  x <- as.matrix(data.frame(lapply(raw[-1, -1, drop = FALSE], as.numeric),
                            check.names = FALSE))
  rownames(x) <- raw[-1, 1]
  list(x = x, groups = groups, cdf = NULL)
}
restrict_to_all_mll <- function(dat) {
  keep <- dat$groups %in% c("ALL", "MLL")
  if (exclude_u95av2) {
    drop_idx <- cdf_safe_drop_u95av2(dat$groups, dat$cdf)
    if (length(drop_idx)) keep[drop_idx] <- FALSE
  }
  list(x = dat$x[, keep, drop = FALSE],
       groups = factor(dat$groups[keep], levels = c("ALL", "MLL")))
}
apply_transform <- function(x_raw, transform) {
  if (transform == "log2_floor20")
    list(x_log = log2(pmax(x_raw, 20)), x_pos = pmax(x_raw, 20))
  else if (transform == "log2_clip")
    list(x_log = log2(pmax(x_raw, 1) + 1), x_pos = pmax(x_raw, 1) + 1)
  else if (transform == "log")
    list(x_log = log(pmax(x_raw, 1)), x_pos = pmax(x_raw, 1))
  else stop("Unknown transform: ", transform)
}

make_two_group_pattern <- function(groups) {
  patterns <- matrix(c(0, 0, 0, 1), 2, 2)
  colnames(patterns) <- levels(factor(groups))
  patterns
}
fit_gaga_calls <- function(fit_data, groups, nclust) {
  patterns <- make_two_group_pattern(groups)
  fit <- gaga::fitGG(fit_data$x_log, groups, patterns = patterns,
                     equalcv = equalcv, nclust = nclust,
                     method = "quickEM", trace = FALSE)
  fit <- gaga::parest(fit, x = fit_data$x_log, groups = groups, alpha = 0.05)
  g <- gaga::findgenes(fit, fit_data$x_log, groups, fdrmax = 0.05,
                       parametric = TRUE)
  rownames(fit_data$x_log)[as.logical(g$d != 0)]
}
fit_method_calls <- function(method, fit_data, groups) {
  if (method == "GaGa")    return(fit_gaga_calls(fit_data, groups, 1))
  if (method == "MiGaGa2") return(fit_gaga_calls(fit_data, groups, 2))
  if (method == "limma_BH") {
    res <- limma_BH(fit_data$x_log, groups, fdr = 0.05)
    return(rownames(fit_data$x_log)[res$calls])
  }
  if (method == "Ga") {
    res <- ebarrays_GG_raw(fit_data$x_pos, groups, fdr = 0.05)
    return(rownames(fit_data$x_log)[res$calls])
  }
  stop("Unknown method: ", method)
}
safe_fit <- function(method, fit_data, groups) {
  tryCatch({
    list(calls = fit_method_calls(method, fit_data, groups),
         error = NA_character_)
  }, error = function(e) {
    list(calls = character(0), error = conditionMessage(e))
  })
}

# ---------------------------------------------------------------------------
# Setup paths
# ---------------------------------------------------------------------------
suffix <- paste(data_choice, transform_choice,
                if (equalcv) "equalcv" else "freecv",
                if (exclude_u95av2) "exclu95av2" else "allmll",
                sep = "_")
chunk_dir <- result_file(file.path("table1_chunks", suffix))
dir.create(chunk_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Load + transform
# ---------------------------------------------------------------------------
if (data_choice == "orange_full") {
  data_file <- file.path(bench, "data", "MLL_orange_12533.tab")
  dat <- load_armstrong_orange(data_file)
} else if (data_choice == "schliep_filtered") {
  data_file <- file.path(bench, "data", "armstrong-2002-v2_database.txt")
  dat <- load_armstrong_schliep(data_file)
} else {
  stop("Unknown ARMSTRONG_DATA: ", data_choice)
}
dat <- restrict_to_all_mll(dat)
groups <- dat$groups
tx <- apply_transform(dat$x, transform_choice)
rownames(tx$x_log) <- rownames(dat$x); rownames(tx$x_pos) <- rownames(dat$x)

# ---------------------------------------------------------------------------
# Full-data fits (cached)
# ---------------------------------------------------------------------------
full_file <- file.path(chunk_dir, "full.rds")
if (file.exists(full_file)) {
  full <- readRDS(full_file)
  cat("Loaded cached full-data fits from", full_file, "\n")
} else {
  cat("Computing full-data fits ...\n")
  full <- list(lists = list(), elapsed = list(), error = list())
  for (m in methods) {
    el <- system.time({ out <- safe_fit(m, tx, groups) })[["elapsed"]]
    full$lists[[m]] <- out$calls
    full$elapsed[[m]] <- el
    full$error[[m]] <- out$error
    cat(sprintf("  %s: n_de=%d  elapsed=%.2fs  error=%s\n",
                m, length(out$calls), el,
                ifelse(is.na(out$error), "-", out$error)))
  }
  saveRDS(full, full_file)
}

# ---------------------------------------------------------------------------
# Per-rep loop with time budget
# ---------------------------------------------------------------------------
sample_sizes <- c(5, 10, 15)
all_idx <- which(groups == "ALL"); mll_idx <- which(groups == "MLL")

t_start <- Sys.time()
done_reps <- vapply(seq_len(total_reps), function(r) {
  file.exists(file.path(chunk_dir, sprintf("rep_%03d.rds", r)))
}, logical(1))
todo <- which(!done_reps)
cat(sprintf("Reps total=%d  done=%d  todo=%d  budget=%.0fs\n",
            total_reps, sum(done_reps), length(todo), budget_sec))

for (rep_id in todo) {
  if (as.numeric(Sys.time() - t_start, units = "secs") > budget_sec) {
    cat("Time budget exhausted; exiting (resume by re-running).\n"); break
  }
  set.seed(base_seed + rep_id)
  ap <- sample(all_idx); mp <- sample(mll_idx)
  rep_rows <- list(); ix <- 1
  for (ss in sample_sizes) {
    cols <- c(ap[seq_len(ss)], mp[seq_len(ss)])
    sub_data <- list(x_log = tx$x_log[, cols, drop = FALSE],
                     x_pos = tx$x_pos[, cols, drop = FALSE])
    sub_grp <- factor(rep(c("ALL", "MLL"), each = ss),
                      levels = c("ALL", "MLL"))
    for (m in methods) {
      el <- system.time({ out <- safe_fit(m, sub_data, sub_grp) })[["elapsed"]]
      repro <- if (is.na(out$error) && length(out$calls) > 0)
        mean(out$calls %in% full$lists[[m]]) else NA_real_
      rep_rows[[ix]] <- data.frame(
        method = m, n_per_group = as.character(ss), rep = rep_id,
        n_de = if (is.na(out$error)) length(out$calls) else NA_integer_,
        reproducibility = repro, elapsed = el, error = out$error,
        stringsAsFactors = FALSE
      )
      ix <- ix + 1
    }
  }
  saveRDS(do.call(rbind, rep_rows),
          file.path(chunk_dir, sprintf("rep_%03d.rds", rep_id)))
  cat(sprintf("rep %d/%d done  cum=%.1fs\n", rep_id, total_reps,
              as.numeric(Sys.time() - t_start, units = "secs")))
}

cat(sprintf("Chunk run finished. Reps complete now=%d/%d\n",
            sum(vapply(seq_len(total_reps), function(r)
              file.exists(file.path(chunk_dir, sprintf("rep_%03d.rds", r))),
              logical(1))),
            total_reps))
