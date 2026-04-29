#!/usr/bin/env Rscript
# =============================================================================
# run_sim62_chunk.R - resumable chunked runner for Section 6.2 simulation.
# Each call processes pending reps until the time budget hits, then exits.
# Per-rep state is written under
# Claude/results/sim62_chunks/<suffix>/{parametric,nonparametric}/rep_<id>.rds
# A separate merge step aggregates them into the Table 2 / ROC summaries.
# =============================================================================

args <- commandArgs(FALSE)
this_file <- sub("^--file=", "", args[grepl("^--file=", args)][1])
src_dir <- dirname(normalizePath(this_file, mustWork = TRUE))
source(file.path(src_dir, "benchmark_lib.R"))
bench <- setup_benchmark()
require_pkgs(c("gaga", "EBarrays", "Biobase", "limma"))

data_choice <- Sys.getenv("ARM62_DATA", unset = "schliep_filtered")
transform_choice <- Sys.getenv(
  "ARM62_TRANSFORM",
  unset = if (data_choice == "schliep_filtered") "log" else "log2_floor20"
)
fit_method <- Sys.getenv("ARM62_GAGA_METHOD", unset = "quickEM")
n_reps <- env_int("ARM62_REPS", 200)
seed <- env_int("ARM62_SEED", 620)
fdr_target <- as.numeric(Sys.getenv("ARM62_FDR", unset = "0.05"))
sample_sizes <- as.integer(strsplit(
  Sys.getenv("ARM62_SAMPLE_SIZES", unset = "5,10,15,20"), ",")[[1]])
sample_sizes <- sort(unique(sample_sizes))
budget_sec <- as.numeric(Sys.getenv("CHUNK_BUDGET_SEC", unset = "30"))
equalcv <- Sys.getenv("ARM62_EQUALCV", unset = "1") != "0"
sim_filter <- Sys.getenv("ARM62_SIM_TYPES",
                         unset = "parametric,nonparametric")
sim_types <- trimws(strsplit(sim_filter, ",")[[1]])

# ---------------------------------------------------------------------------
# Loaders
# ---------------------------------------------------------------------------
load_armstrong_schliep <- function(f) {
  raw <- read.delim(f, check.names = FALSE, stringsAsFactors = FALSE)
  groups <- as.character(unlist(raw[1, -1]))
  x <- as.matrix(data.frame(lapply(raw[-1, -1, drop = FALSE], as.numeric),
                            check.names = FALSE))
  rownames(x) <- raw[-1, 1]
  list(x = x, groups = groups, cdf = NULL)
}
load_armstrong_orange <- function(f) {
  raw <- read.delim(f, check.names = FALSE, stringsAsFactors = FALSE)
  groups <- raw$class[-c(1, 2)]
  xdf <- raw[-c(1, 2), -1, drop = FALSE]
  x <- t(apply(xdf, 2, as.numeric))
  rownames(x) <- colnames(raw)[-1]
  list(x = x, groups = groups, cdf = NULL)
}
restrict <- function(d) {
  keep <- d$groups %in% c("ALL", "MLL")
  drop_idx <- cdf_safe_drop_u95av2(d$groups, d$cdf)
  if (length(drop_idx)) keep[drop_idx] <- FALSE
  list(x = d$x[, keep, drop = FALSE],
       groups = factor(d$groups[keep], levels = c("ALL", "MLL")))
}
apply_t <- function(x_raw, t) {
  if (t == "log2_floor20")
    list(x_log = log2(pmax(x_raw, 20)), x_pos = pmax(x_raw, 20), undo = function(z) 2^z)
  else if (t == "log2_clip")
    list(x_log = log2(pmax(x_raw, 1) + 1), x_pos = pmax(x_raw, 1) + 1, undo = function(z) 2^z)
  else if (t == "log")
    list(x_log = log(pmax(x_raw, 1)), x_pos = pmax(x_raw, 1), undo = function(z) exp(z))
  else stop("Unknown transform: ", t)
}
two_grp <- function(g) {
  patterns <- matrix(c(0, 0, 0, 1), 2, 2)
  colnames(patterns) <- levels(factor(g)); patterns
}
fit_gaga <- function(x_log, groups, nclust) {
  fit <- gaga::fitGG(x_log, groups, patterns = two_grp(groups),
                     equalcv = equalcv, nclust = nclust,
                     method = fit_method, trace = FALSE)
  gaga::parest(fit, x = x_log, groups = groups, alpha = fdr_target)
}
score_calls <- function(method, fit_data, groups) {
  out <- tryCatch({
    if (method == "GaGa") {
      f <- fit_gaga(fit_data$x_log, groups, 1)
      g <- gaga::findgenes(f, fit_data$x_log, groups,
                           fdrmax = fdr_target, parametric = TRUE)
      list(calls = as.logical(g$d != 0), score = 1 - f$pp[, 1],
           error = NA_character_)
    } else if (method == "MiGaGa2") {
      f <- fit_gaga(fit_data$x_log, groups, 2)
      g <- gaga::findgenes(f, fit_data$x_log, groups,
                           fdrmax = fdr_target, parametric = TRUE)
      list(calls = as.logical(g$d != 0), score = 1 - f$pp[, 1],
           error = NA_character_)
    } else if (method == "Ga") {
      r <- ebarrays_GG_raw(fit_data$x_pos, groups, fdr = fdr_target)
      list(calls = r$calls, score = r$score, error = NA_character_)
    } else if (method == "limma_BH") {
      r <- limma_BH(fit_data$x_log, groups, fdr = fdr_target)
      list(calls = r$calls, score = r$score, error = NA_character_)
    } else stop("Unknown method: ", method)
  }, error = function(e) list(calls = rep(FALSE, nrow(fit_data$x_log)),
                              score = rep(NA_real_, nrow(fit_data$x_log)),
                              error = conditionMessage(e)))
  names(out$calls) <- rownames(fit_data$x_log)
  names(out$score) <- rownames(fit_data$x_log); out
}
roc_curve <- function(score, truth) {
  ok <- is.finite(score); score <- score[ok]; truth <- truth[ok]
  ord <- order(score, decreasing = TRUE); truth <- truth[ord]
  tp <- cumsum(truth); fp <- cumsum(!truth)
  data.frame(rank = seq_along(score),
             fdr = fp / pmax(tp + fp, 1),
             power = tp / max(sum(truth), 1))
}
power_at_grid <- function(curve, grid) {
  vapply(grid, function(g) {
    ok <- curve$fdr <= g
    if (!any(ok)) 0 else max(curve$power[ok], na.rm = TRUE)
  }, numeric(1))
}

# ---------------------------------------------------------------------------
# Load data, do full-data fit
# ---------------------------------------------------------------------------
suffix <- paste(data_choice, transform_choice, paste0("reps", n_reps),
                if (equalcv) "equalcv" else "freecv", sep = "_")
chunk_dir <- result_file(file.path("sim62_chunks", suffix))
dir.create(chunk_dir, recursive = TRUE, showWarnings = FALSE)

if (data_choice == "schliep_filtered") {
  dat <- restrict(load_armstrong_schliep(file.path(bench, "data",
                  "armstrong-2002-v2_database.txt")))
} else {
  dat <- restrict(load_armstrong_orange(file.path(bench, "data",
                  "MLL_orange_12533.tab")))
}
groups <- dat$groups; tx <- apply_t(dat$x, transform_choice)
rownames(tx$x_log) <- rownames(dat$x); rownames(tx$x_pos) <- rownames(dat$x)

full_file <- file.path(chunk_dir, "full.rds")
if (file.exists(full_file)) {
  full <- readRDS(full_file)
  cat("[cache] loaded full-data fit\n")
} else {
  cat("[fit ] full-data GaGa for simulation truth ...\n")
  set.seed(seed)
  ff <- fit_gaga(tx$x_log, groups, 1)
  fg <- gaga::findgenes(ff, tx$x_log, groups,
                        fdrmax = fdr_target, parametric = TRUE)
  full <- list(par = gaga::getpar(ff),
               truth = as.logical(fg$d != 0),
               n_de = sum(fg$d != 0))
  names(full$truth) <- rownames(tx$x_log)
  saveRDS(full, full_file)
}
cat("Full-data DE genes:", full$n_de, "/", length(full$truth), "\n")

# ---------------------------------------------------------------------------
# Simulators
# ---------------------------------------------------------------------------
sim_parametric <- function(par, n_genes, m_per_group) {
  p_de <- unname(par$probpat[min(2, length(par$probpat))])
  p_de <- min(max(p_de, 0), 1)
  xsim <- gaga::simGG(n = n_genes, m = c(m_per_group, m_per_group),
                      p.de = p_de,
                      a0 = par$a0, nu = par$nu,
                      balpha = par$balpha, nualpha = par$nualpha,
                      equalcv = equalcv, probclus = par$probclus)
  fd <- Biobase::fData(xsim)
  list(x = Biobase::exprs(xsim),
       truth = abs(fd$mean.1 - fd$mean.2) > 1e-12)
}
sim_nonparametric <- function(x_log, groups, truth, m_per_group) {
  all_cols <- which(groups == levels(groups)[1])
  mll_cols <- which(groups == levels(groups)[2])
  idx_ee  <- sample(seq_len(ncol(x_log)), 2 * m_per_group, replace = TRUE)
  idx_all <- sample(all_cols, m_per_group, replace = TRUE)
  idx_mll <- sample(mll_cols, m_per_group, replace = TRUE)
  out <- matrix(NA_real_, nrow = nrow(x_log), ncol = 2 * m_per_group)
  out[!truth, ] <- x_log[!truth, idx_ee, drop = FALSE]
  out[truth, seq_len(m_per_group)] <-
    x_log[truth, idx_all, drop = FALSE]
  out[truth, m_per_group + seq_len(m_per_group)] <-
    x_log[truth, idx_mll, drop = FALSE]
  out
}
subset_sim <- function(xsim, ss) {
  max_m <- ncol(xsim) / 2
  xsim[, c(seq_len(ss), max_m + seq_len(ss)), drop = FALSE]
}

# ---------------------------------------------------------------------------
# Main loop
# ---------------------------------------------------------------------------
methods <- c("GaGa", "MiGaGa2", "Ga", "limma_BH")
roc_grid <- seq(0, 1, by = 0.01)
max_m <- max(sample_sizes)

for (st in sim_types) dir.create(file.path(chunk_dir, st), recursive = TRUE,
                                 showWarnings = FALSE)
done_for <- function(st, r)
  file.exists(file.path(chunk_dir, st, sprintf("rep_%03d.rds", r)))

t_start <- Sys.time()
for (rep_id in seq_len(n_reps)) {
  if (as.numeric(Sys.time() - t_start, units = "secs") > budget_sec) {
    cat("Time budget exhausted; resume by re-running.\n"); break
  }
  for (st in sim_types) {
    if (done_for(st, rep_id)) next
    set.seed(seed + 31 * which(sim_types == st)[1] + rep_id)
    if (st == "parametric") {
      sim <- sim_parametric(full$par, length(full$truth), max_m)
      truth_rep <- sim$truth
      xsim_log <- sim$x
    } else {
      xsim_log <- sim_nonparametric(tx$x_log, groups, full$truth, max_m)
      truth_rep <- full$truth
    }
    xsim_pos <- tx$undo(xsim_log)
    rownames(xsim_log) <- rownames(tx$x_log)
    rownames(xsim_pos) <- rownames(tx$x_log)
    names(truth_rep)   <- rownames(tx$x_log)
    rep_rows <- list(); roc_rows <- list(); ix <- 1; rix <- 1
    for (ss in sample_sizes) {
      sub_log <- subset_sim(xsim_log, ss)
      sub_pos <- subset_sim(xsim_pos, ss)
      sub_grp <- factor(rep(levels(groups), each = ss),
                        levels = levels(groups))
      fit_data <- list(x_log = sub_log, x_pos = sub_pos)
      for (m in methods) {
        elapsed <- system.time({
          sc <- score_calls(m, fit_data, sub_grp)
        })[["elapsed"]]
        ev <- eval_calls(sc$calls, truth_rep)
        rep_rows[[ix]] <- data.frame(
          simulation = st, rep = rep_id, n_per_group = ss,
          method = m, ev, elapsed = elapsed, error = sc$error,
          stringsAsFactors = FALSE)
        ix <- ix + 1
        if (st == "nonparametric" && ss == max_m && is.na(sc$error)) {
          curve <- roc_curve(sc$score, truth_rep)
          pow <- power_at_grid(curve, roc_grid)
          roc_rows[[rix]] <- data.frame(
            rep = rep_id, method = m, fdr_grid = roc_grid, power = pow,
            stringsAsFactors = FALSE)
          rix <- rix + 1
        }
      }
    }
    saveRDS(list(rows = do.call(rbind, rep_rows),
                 roc = if (length(roc_rows)) do.call(rbind, roc_rows)
                       else NULL),
            file.path(chunk_dir, st, sprintf("rep_%03d.rds", rep_id)))
  }
  cat(sprintf("rep %d/%d done  cum=%.1fs\n", rep_id, n_reps,
              as.numeric(Sys.time() - t_start, units = "secs")))
}

done_all <- vapply(seq_len(n_reps), function(r) {
  all(vapply(sim_types, done_for, logical(1), r = r))
}, logical(1))
cat(sprintf("Reps complete now=%d/%d\n", sum(done_all), n_reps))
