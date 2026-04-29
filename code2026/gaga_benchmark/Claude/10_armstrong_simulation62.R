#!/usr/bin/env Rscript
# =============================================================================
# 10_armstrong_simulation62.R  (Claude/ debugged version)
#
# Reproduces Rossell (2009) Section 6.2 simulations:
#   * Full-data GaGa fit defines the "truth" mask and parametric generator
#   * Parametric simulation: simGG with fitted hyperparameters
#   * Nonparametric simulation: bootstrap from ALL/MLL columns for DE genes,
#     bootstrap from any column for EE genes
#   * For each (rep, sample-size, method) compute FDR / power vs truth
#   * Average ROC curves across reps for the n=20 nonparametric case
#
# Fixes vs gaga_benchmark/R/10_armstrong_simulation62.R:
#   * Ga (EBarrays GG) is fit on the *positive-scale* matrix (2^x_log or
#     exp(x_log)), not on log-scale. The legacy script gave it log-scale data.
#   * limma uses limma_BH on the log-scale matrix directly (no double-log,
#     no fragile 2^x trick that broke for natural-log transforms).
#   * Default sample sizes are 5, 10, 15, 20 (paper) but configurable.
# =============================================================================

args <- commandArgs(FALSE)
this_file <- sub("^--file=", "", args[grepl("^--file=", args)][1])
source(file.path(dirname(normalizePath(this_file, mustWork = TRUE)),
                 "benchmark_lib.R"))

bench <- setup_benchmark()
paths <- .libPaths()
sys_paths <- paths[grepl("^/opt/homebrew", paths)]
loc_paths <- paths[!grepl("^/opt/homebrew", paths)]
.libPaths(unique(c(sys_paths, file.path(bench, "r-lib-gcrma"), loc_paths)))
require_pkgs(c("gaga", "EBarrays", "Biobase", "limma"))

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
data_choice <- Sys.getenv("ARM62_DATA", unset = "schliep_filtered")
transform_choice <- Sys.getenv(
  "ARM62_TRANSFORM",
  unset = if (data_choice == "schliep_filtered") "log" else "log2_floor20"
)
fit_method <- Sys.getenv("ARM62_GAGA_METHOD", unset = "quickEM")
n_reps <- env_int("ARM62_REPS", 20)
seed <- env_int("ARM62_SEED", 620)
fdr_target <- as.numeric(Sys.getenv("ARM62_FDR", unset = "0.05"))
sample_sizes <- as.integer(strsplit(
  Sys.getenv("ARM62_SAMPLE_SIZES", unset = "5,10,15,20"), ",")[[1]])
sample_sizes <- sort(unique(sample_sizes))
roc_grid <- seq(0, 1, by = 0.01)
equalcv <- Sys.getenv("ARM62_EQUALCV", unset = "1") != "0"

# ---------------------------------------------------------------------------
# Loaders / transforms (shared with 06; copied to keep this file standalone)
# ---------------------------------------------------------------------------
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

restrict_to_all_mll <- function(dat, exclude_u95av2 = TRUE) {
  keep <- dat$groups %in% c("ALL", "MLL")
  if (exclude_u95av2) {
    drop_idx <- cdf_safe_drop_u95av2(dat$groups, dat$cdf)
    if (length(drop_idx) > 0) keep[drop_idx] <- FALSE
  }
  list(
    x = dat$x[, keep, drop = FALSE],
    groups = factor(dat$groups[keep], levels = c("ALL", "MLL"))
  )
}

apply_transform <- function(x_raw, transform) {
  if (transform == "log2_floor20") {
    list(x_log = log2(pmax(x_raw, 20)),
         x_pos = pmax(x_raw, 20), base = 2)
  } else if (transform == "log2_clip") {
    list(x_log = log2(pmax(x_raw, 1) + 1),
         x_pos = pmax(x_raw, 1) + 1, base = 2)
  } else if (transform == "log") {
    list(x_log = log(pmax(x_raw, 1)),
         x_pos = pmax(x_raw, 1), base = exp(1))
  } else {
    stop("Unknown transform: ", transform)
  }
}

# ---------------------------------------------------------------------------
# Method wrappers
# ---------------------------------------------------------------------------
gaga_two_group <- function(groups) {
  patterns <- matrix(c(0, 0, 0, 1), 2, 2)
  colnames(patterns) <- levels(factor(groups))
  patterns
}

fit_gaga <- function(x_log, groups, nclust) {
  fit <- gaga::fitGG(x_log, groups, patterns = gaga_two_group(groups),
                     equalcv = equalcv, nclust = nclust,
                     method = fit_method, trace = FALSE)
  gaga::parest(fit, x = x_log, groups = groups, alpha = fdr_target)
}

score_calls <- function(method, fit_data, groups) {
  out <- tryCatch({
    if (method == "GaGa") {
      fit <- fit_gaga(fit_data$x_log, groups, nclust = 1)
      genes <- gaga::findgenes(fit, fit_data$x_log, groups,
                               fdrmax = fdr_target, parametric = TRUE)
      list(calls = as.logical(genes$d != 0),
           score = 1 - fit$pp[, 1], error = NA_character_)
    } else if (method == "MiGaGa2") {
      fit <- fit_gaga(fit_data$x_log, groups, nclust = 2)
      genes <- gaga::findgenes(fit, fit_data$x_log, groups,
                               fdrmax = fdr_target, parametric = TRUE)
      list(calls = as.logical(genes$d != 0),
           score = 1 - fit$pp[, 1], error = NA_character_)
    } else if (method == "Ga") {
      ga <- ebarrays_GG_raw(fit_data$x_pos, groups, fdr = fdr_target)
      list(calls = ga$calls, score = ga$score, error = NA_character_)
    } else if (method == "limma_BH") {
      lim <- limma_BH(fit_data$x_log, groups, fdr = fdr_target)
      list(calls = lim$calls, score = lim$score, error = NA_character_)
    } else {
      stop("Unknown method: ", method)
    }
  }, error = function(e) {
    list(calls = rep(FALSE, nrow(fit_data$x_log)),
         score = rep(NA_real_, nrow(fit_data$x_log)),
         error = conditionMessage(e))
  })
  names(out$calls) <- rownames(fit_data$x_log)
  names(out$score) <- rownames(fit_data$x_log)
  out
}

# ---------------------------------------------------------------------------
# ROC helpers
# ---------------------------------------------------------------------------
roc_curve <- function(score, truth) {
  ok <- is.finite(score)
  score <- score[ok]; truth <- truth[ok]
  ord <- order(score, decreasing = TRUE)
  truth <- truth[ord]
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

auc_grid <- function(x, y) sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)

# ---------------------------------------------------------------------------
# Simulators
# ---------------------------------------------------------------------------
simulate_parametric <- function(fit, n_genes, m_per_group) {
  par <- gaga::getpar(fit)
  p_de <- unname(par$probpat[min(2, length(par$probpat))])
  p_de <- min(max(p_de, 0), 1)
  xsim <- gaga::simGG(
    n = n_genes, m = c(m_per_group, m_per_group),
    p.de = p_de,
    a0 = par$a0, nu = par$nu,
    balpha = par$balpha, nualpha = par$nualpha,
    equalcv = equalcv, probclus = par$probclus
  )
  fd <- Biobase::fData(xsim)
  truth <- abs(fd$mean.1 - fd$mean.2) > 1e-12
  list(x = Biobase::exprs(xsim), truth = truth)
}

simulate_nonparametric <- function(x_log, groups, truth, m_per_group) {
  all_cols <- which(groups == levels(groups)[1])
  mll_cols <- which(groups == levels(groups)[2])
  idx_ee <- sample(seq_len(ncol(x_log)), 2 * m_per_group, replace = TRUE)
  idx_all <- sample(all_cols, m_per_group, replace = TRUE)
  idx_mll <- sample(mll_cols, m_per_group, replace = TRUE)
  out <- matrix(NA_real_, nrow = nrow(x_log), ncol = 2 * m_per_group)
  out[!truth, ] <- x_log[!truth, idx_ee, drop = FALSE]
  out[truth,  seq_len(m_per_group)] <-
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
# Load + transform
# ---------------------------------------------------------------------------
if (data_choice == "schliep_filtered") {
  dat <- load_armstrong_schliep(file.path(bench, "data",
                                          "armstrong-2002-v2_database.txt"))
} else if (data_choice == "orange_full") {
  dat <- load_armstrong_orange(file.path(bench, "data",
                                         "MLL_orange_12533.tab"))
} else {
  stop("Unknown ARM62_DATA: ", data_choice)
}
dat <- restrict_to_all_mll(dat)

set.seed(seed)
tx <- apply_transform(dat$x, transform_choice)
x_log <- tx$x_log; rownames(x_log) <- rownames(dat$x)
x_pos <- tx$x_pos; rownames(x_pos) <- rownames(dat$x)
groups <- dat$groups

cat("Armstrong 6.2 data:\n")
cat("  data:", data_choice, "\n")
cat("  transform:", transform_choice, "\n")
cat("  dimensions:", nrow(x_log), "genes x", ncol(x_log), "samples\n")
cat("  groups:", paste(names(table(groups)),
                       as.integer(table(groups)), collapse = ", "), "\n")
cat("  reps =", n_reps, " sample_sizes =",
    paste(sample_sizes, collapse = ","), "\n")

cat("Fitting full-data GaGa to define simulation truth...\n")
full_fit <- fit_gaga(x_log, groups, nclust = 1)
full_genes <- gaga::findgenes(full_fit, x_log, groups,
                              fdrmax = fdr_target, parametric = TRUE)
full_truth <- as.logical(full_genes$d != 0)
names(full_truth) <- rownames(x_log)
full_par <- gaga::getpar(full_fit)
cat("  full-data DE genes:", sum(full_truth), "of", length(full_truth), "\n")
cat("  fitted prior DE prob:",
    signif(unname(full_par$probpat[min(2, length(full_par$probpat))]), 4), "\n")

# ---------------------------------------------------------------------------
# Simulation loop
# ---------------------------------------------------------------------------
max_m <- max(sample_sizes)
methods <- c("GaGa", "MiGaGa2", "Ga", "limma_BH")
rows <- list(); roc_rows <- list()
ix <- 1; rix <- 1

for (sim_type in c("parametric", "nonparametric")) {
  for (rep_id in seq_len(n_reps)) {
    if (sim_type == "parametric") {
      sim <- simulate_parametric(full_fit, nrow(x_log), max_m)
      truth_rep <- sim$truth
      xsim_log <- sim$x
    } else {
      xsim_log <- simulate_nonparametric(x_log, groups, full_truth, max_m)
      truth_rep <- full_truth
    }
    # Restore positive-scale companion for Ga/GG
    xsim_pos <- if (transform_choice == "log") {
      exp(xsim_log)
    } else {
      2 ^ xsim_log
    }
    rownames(xsim_log) <- rownames(x_log)
    rownames(xsim_pos) <- rownames(x_log)
    names(truth_rep)   <- rownames(x_log)
    colnames(xsim_log) <- paste0(rep(levels(groups), each = max_m), "_",
                                 sequence(rep(max_m, 2)))
    colnames(xsim_pos) <- colnames(xsim_log)

    for (ss in sample_sizes) {
      sub_log <- subset_sim(xsim_log, ss)
      sub_pos <- subset_sim(xsim_pos, ss)
      sub_grp <- factor(rep(levels(groups), each = ss),
                        levels = levels(groups))
      fit_data <- list(x_log = sub_log, x_pos = sub_pos)
      for (method in methods) {
        elapsed <- system.time({
          sc <- score_calls(method, fit_data, sub_grp)
        })[["elapsed"]]
        ev <- eval_calls(sc$calls, truth_rep)
        rows[[ix]] <- data.frame(
          simulation = sim_type, rep = rep_id, n_per_group = ss,
          method = method, ev,
          elapsed = elapsed, error = sc$error,
          stringsAsFactors = FALSE
        )
        ix <- ix + 1
        if (sim_type == "nonparametric" && ss == max_m && is.na(sc$error)) {
          curve <- roc_curve(sc$score, truth_rep)
          pow <- power_at_grid(curve, roc_grid)
          roc_rows[[rix]] <- data.frame(
            rep = rep_id, method = method,
            fdr_grid = roc_grid, power = pow, stringsAsFactors = FALSE
          )
          rix <- rix + 1
        }
      }
    }
    cat(sprintf("  %s rep %d/%d done\n", sim_type, rep_id, n_reps))
  }
}

result_rows <- do.call(rbind, rows)
result_rows$has_error <- !is.na(result_rows$error)

summary_rows <- aggregate(
  cbind(n_de, true_de, fdr, power, elapsed) ~
    simulation + method + n_per_group,
  data = result_rows, FUN = function(v) mean(v, na.rm = TRUE)
)
errors <- aggregate(has_error ~ simulation + method + n_per_group,
                    data = result_rows, FUN = sum)
names(errors)[4] <- "n_errors"
summary_rows <- merge(summary_rows, errors,
                      by = c("simulation", "method", "n_per_group"),
                      all.x = TRUE)

# Paper Table 2 values (Rossell 2009)
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
    0.047, 0.021, 0.019, 0.024
  ),
  paper_power = c(
    0.066, 0.322, 0.512, 0.608,
    0.066, 0.328, 0.520, 0.615,
    0.434, 0.587, 0.667, 0.712,
    0.063, 0.288, 0.487, 0.580,
    0.054, 0.319, 0.529, 0.660,
    0.057, 0.327, 0.541, 0.671,
    0.397, 0.567, 0.666, 0.740,
    0.029, 0.168, 0.321, 0.449
  )
)
summary_rows <- merge(summary_rows, paper_table2,
                      by = c("simulation", "method", "n_per_group"),
                      all.x = TRUE)
summary_rows <- summary_rows[order(summary_rows$simulation,
                                   summary_rows$method,
                                   summary_rows$n_per_group), ]

roc_all <- do.call(rbind, roc_rows)
roc_avg <- aggregate(power ~ method + fdr_grid, data = roc_all, FUN = mean)
roc_avg <- roc_avg[order(roc_avg$method, roc_avg$fdr_grid), ]
roc_summary <- do.call(rbind, lapply(split(roc_avg, roc_avg$method),
  function(z) data.frame(method = z$method[1],
                         auc = auc_grid(z$fdr_grid, z$power),
                         stringsAsFactors = FALSE)
))
rownames(roc_summary) <- NULL

suffix <- paste(data_choice, transform_choice, paste0("reps", n_reps),
                if (equalcv) "equalcv" else "freecv", sep = "_")
write.csv(result_rows,
          result_file(paste0("armstrong62_simulation_rows_", suffix, ".csv")),
          row.names = FALSE)
write.csv(summary_rows,
          result_file(paste0("armstrong62_table2_summary_", suffix, ".csv")),
          row.names = FALSE)
write.csv(roc_avg,
          result_file(paste0("armstrong62_roc_curves_", suffix, ".csv")),
          row.names = FALSE)
write.csv(roc_summary,
          result_file(paste0("armstrong62_roc_summary_", suffix, ".csv")),
          row.names = FALSE)
saveRDS(list(summary = summary_rows, rows = result_rows,
             roc = roc_avg, full_truth = full_truth, full_fit = full_fit),
        result_file(paste0("armstrong62_simulation_", suffix, ".rds")))

plot_roc <- function(file, device = c("pdf", "png")) {
  device <- match.arg(device)
  if (device == "pdf") grDevices::pdf(file, width = 6.2, height = 5.2)
  else grDevices::png(file, width = 1400, height = 1100, res = 200)
  cols <- c(GaGa = "#1B9E77", MiGaGa2 = "#D95F02",
            Ga = "#666666", limma_BH = "#7570B3")
  plot(NA, xlim = c(0, 1), ylim = c(0, 1),
       xlab = "Average FDR", ylab = "Average power",
       main = "Armstrong nonparametric simulations")
  for (m in methods) {
    z <- roc_avg[roc_avg$method == m, ]
    lines(z$fdr_grid, z$power, col = cols[[m]], lwd = 2)
  }
  leg <- merge(data.frame(method = methods), roc_summary,
               by = "method", all.x = TRUE, sort = FALSE)
  legend("bottomright",
         legend = paste0(leg$method, " AUC=", sprintf("%.3f", leg$auc)),
         col = cols[leg$method], lwd = 2, bty = "n")
  grDevices::dev.off()
}
plot_roc(result_file(paste0("armstrong62_roc_curves_", suffix, ".pdf")), "pdf")
plot_roc(result_file(paste0("armstrong62_roc_curves_", suffix, ".png")), "png")

cat("\n=== Table 2 reproduction (paper values merged) ===\n")
print(summary_rows, row.names = FALSE)
cat("\n=== ROC AUC summary ===\n")
print(roc_summary, row.names = FALSE)
cat("\nOutputs suffix:", suffix, "\n")
