args <- commandArgs(FALSE)
this_file <- sub("^--file=", "", args[grepl("^--file=", args)][1])
source(file.path(dirname(normalizePath(this_file, mustWork = TRUE)), "benchmark_lib.R"))

bench <- setup_benchmark()
paths <- .libPaths()
system_paths <- paths[grepl("^/opt/homebrew", paths)]
local_paths <- paths[!grepl("^/opt/homebrew", paths)]
.libPaths(unique(c(system_paths, file.path(bench, "r-lib-gcrma"), local_paths)))
require_pkgs(c("gaga", "EBarrays", "Biobase", "limma"))

data_choice <- Sys.getenv("ARM62_DATA", unset = "schliep_filtered")
transform_choice <- Sys.getenv("ARM62_TRANSFORM", unset = if (data_choice == "schliep_filtered") "log" else "log2_floor20")
fit_method <- Sys.getenv("ARM62_GAGA_METHOD", unset = "quickEM")
n_reps <- env_int("ARM62_REPS", 20)
seed <- env_int("ARM62_SEED", 620)
fdr_target <- as.numeric(Sys.getenv("ARM62_FDR", unset = "0.05"))
sample_sizes <- as.integer(strsplit(Sys.getenv("ARM62_SAMPLE_SIZES", unset = "5,10,15,20"), ",")[[1]])
sample_sizes <- sort(unique(sample_sizes))
roc_grid <- seq(0, 1, by = 0.01)

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

transform_expr <- function(x, transform) {
  if (transform == "log") {
    return(log(pmax(x, 1)))
  }
  if (transform == "log1p") {
    return(log(pmax(x, 0) + 1))
  }
  if (transform == "log2_floor20") {
    return(log2(pmax(x, 20)))
  }
  if (transform == "log2p1") {
    return(log2(pmax(x, 0) + 1))
  }
  stop("Unknown ARM62_TRANSFORM: ", transform)
}

gaga_patterns <- function(groups) {
  patterns <- matrix(c(0, 0, 0, 1), 2, 2)
  colnames(patterns) <- levels(factor(groups))
  patterns
}

eb_patterns <- function(groups) {
  EBarrays::ebPatterns(c(
    paste(rep(1, length(groups)), collapse = " "),
    paste(ifelse(groups == levels(factor(groups))[1], 1, 2), collapse = " ")
  ))
}

fit_gaga <- function(x, groups, nclust) {
  fit <- gaga::fitGG(
    x,
    groups,
    patterns = gaga_patterns(groups),
    equalcv = TRUE,
    nclust = nclust,
    method = fit_method,
    trace = FALSE
  )
  gaga::parest(fit, x = x, groups = groups, alpha = fdr_target)
}

fit_ga <- function(x, groups) {
  fit <- EBarrays::emfit(
    data = x,
    family = "GG",
    hypotheses = eb_patterns(groups),
    num.iter = 20,
    verbose = FALSE
  )
  post <- EBarrays::postprob(fit, x)$pattern
  threshold <- EBarrays::crit.fun(post[, 1], fdr_target)
  list(fit = fit, post = post, threshold = threshold)
}

score_calls <- function(method, x, groups) {
  out <- tryCatch({
    if (method == "GaGa") {
      fit <- fit_gaga(x, groups, nclust = 1)
      genes <- gaga::findgenes(fit, x, groups, fdrmax = fdr_target, parametric = TRUE)
      score <- 1 - fit$pp[, 1]
      return(list(calls = as.logical(genes$d != 0), score = score, error = NA_character_))
    }
    if (method == "MiGaGa2") {
      fit <- fit_gaga(x, groups, nclust = 2)
      genes <- gaga::findgenes(fit, x, groups, fdrmax = fdr_target, parametric = TRUE)
      score <- 1 - fit$pp[, 1]
      return(list(calls = as.logical(genes$d != 0), score = score, error = NA_character_))
    }
    if (method == "Ga") {
      fit <- fit_ga(x, groups)
      return(list(calls = fit$post[, 2] > fit$threshold, score = fit$post[, 2], error = NA_character_))
    }
    if (method == "limma_BH") {
      lim <- limma_or_ttest(2^x, groups, fdr = fdr_target)
      return(list(calls = lim$calls, score = lim$score, error = NA_character_))
    }
    stop("Unknown method: ", method)
  }, error = function(e) {
    list(calls = rep(FALSE, nrow(x)), score = rep(NA_real_, nrow(x)), error = conditionMessage(e))
  })
  names(out$calls) <- rownames(x)
  names(out$score) <- rownames(x)
  out
}

roc_curve <- function(score, truth) {
  ok <- is.finite(score)
  score <- score[ok]
  truth <- truth[ok]
  ord <- order(score, decreasing = TRUE)
  truth <- truth[ord]
  tp <- cumsum(truth)
  fp <- cumsum(!truth)
  n_true <- sum(truth)
  fdr <- fp / pmax(tp + fp, 1)
  power <- tp / max(n_true, 1)
  data.frame(rank = seq_along(score), fdr = fdr, power = power)
}

power_at_grid <- function(curve, grid) {
  vapply(grid, function(g) {
    ok <- curve$fdr <= g
    if (!any(ok)) {
      return(0)
    }
    max(curve$power[ok], na.rm = TRUE)
  }, numeric(1))
}

auc_grid <- function(x, y) {
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
}

simulate_parametric <- function(fit, n_genes, m) {
  par <- gaga::getpar(fit)
  p_de <- unname(par$probpat[min(2, length(par$probpat))])
  p_de <- min(max(p_de, 0), 1)
  xsim <- gaga::simGG(
    n = n_genes,
    m = c(m, m),
    p.de = p_de,
    a0 = par$a0,
    nu = par$nu,
    balpha = par$balpha,
    nualpha = par$nualpha,
    equalcv = TRUE,
    probclus = par$probclus
  )
  fd <- Biobase::fData(xsim)
  truth <- abs(fd$mean.1 - fd$mean.2) > 1e-12
  list(x = Biobase::exprs(xsim), truth = truth)
}

simulate_nonparametric <- function(x, groups, truth, m) {
  all_cols <- which(groups == levels(groups)[1])
  mll_cols <- which(groups == levels(groups)[2])
  idx_ee <- sample(seq_len(ncol(x)), 2 * m, replace = TRUE)
  idx_all <- sample(all_cols, m, replace = TRUE)
  idx_mll <- sample(mll_cols, m, replace = TRUE)
  out <- matrix(NA_real_, nrow = nrow(x), ncol = 2 * m)
  out[!truth, ] <- x[!truth, idx_ee, drop = FALSE]
  out[truth, seq_len(m)] <- x[truth, idx_all, drop = FALSE]
  out[truth, m + seq_len(m)] <- x[truth, idx_mll, drop = FALSE]
  out
}

subset_sim <- function(xsim, ss) {
  max_m <- ncol(xsim) / 2
  xsim[, c(seq_len(ss), max_m + seq_len(ss)), drop = FALSE]
}

if (data_choice == "schliep_filtered") {
  dat <- load_armstrong_schliep(file.path(bench, "data", "armstrong-2002-v2_database.txt"))
} else if (data_choice == "orange_full") {
  dat <- load_armstrong_orange(file.path(bench, "data", "MLL_orange_12533.tab"))
} else {
  stop("Unknown ARM62_DATA: ", data_choice)
}

set.seed(seed)
x <- transform_expr(dat$x, transform_choice)
groups <- dat$groups
if (any(!is.finite(x)) || min(x) < 0) {
  stop("Transformed data must be finite and non-negative for gamma-family fits.")
}

cat("Armstrong 6.2 data:\n")
cat("  data:", data_choice, "\n")
cat("  transform:", transform_choice, "\n")
cat("  dimensions:", nrow(x), "genes x", ncol(x), "samples\n")
cat("  groups:", paste(names(table(groups)), as.integer(table(groups)), collapse = ", "), "\n")

cat("Fitting full-data GaGa to define simulation truth and parametric generator...\n")
full_fit <- fit_gaga(x, groups, nclust = 1)
full_genes <- gaga::findgenes(full_fit, x, groups, fdrmax = fdr_target, parametric = TRUE)
full_truth <- as.logical(full_genes$d != 0)
names(full_truth) <- rownames(x)
full_par <- gaga::getpar(full_fit)
cat("  full-data DE genes:", sum(full_truth), "of", length(full_truth), "\n")
cat("  fitted prior DE probability:", signif(unname(full_par$probpat[min(2, length(full_par$probpat))]), 4), "\n")

max_m <- max(sample_sizes)
methods <- c("GaGa", "MiGaGa2", "Ga", "limma_BH")
group_max <- factor(rep(levels(groups), each = max_m), levels = levels(groups))

rows <- list()
roc_rows <- list()
row_id <- 1
roc_id <- 1

for (sim_type in c("parametric", "nonparametric")) {
  for (rep_id in seq_len(n_reps)) {
    sim <- if (sim_type == "parametric") {
      simulate_parametric(full_fit, nrow(x), max_m)
    } else {
      list(x = simulate_nonparametric(x, groups, full_truth, max_m), truth = full_truth)
    }
    xsim_all <- sim$x
    truth <- sim$truth
    rownames(xsim_all) <- rownames(x)
    names(truth) <- rownames(x)
    colnames(xsim_all) <- paste0(rep(levels(groups), each = max_m), "_", sequence(rep(max_m, 2)))
    xsim_all <- pmax(xsim_all, .Machine$double.eps)

    for (ss in sample_sizes) {
      xsub <- subset_sim(xsim_all, ss)
      gsub <- factor(rep(levels(groups), each = ss), levels = levels(groups))
      for (method in methods) {
        elapsed <- system.time({
          sc <- score_calls(method, xsub, gsub)
        })[["elapsed"]]
        ev <- eval_calls(sc$calls, truth)
        rows[[row_id]] <- data.frame(
          simulation = sim_type,
          rep = rep_id,
          n_per_group = ss,
          method = method,
          ev,
          elapsed = elapsed,
          error = sc$error,
          stringsAsFactors = FALSE
        )
        row_id <- row_id + 1
        if (sim_type == "nonparametric" && ss == max_m && is.na(sc$error)) {
          curve <- roc_curve(sc$score, truth)
          pow <- power_at_grid(curve, roc_grid)
          roc_rows[[roc_id]] <- data.frame(
            rep = rep_id,
            method = method,
            fdr_grid = roc_grid,
            power = pow,
            stringsAsFactors = FALSE
          )
          roc_id <- roc_id + 1
        }
      }
    }
    cat("  ", sim_type, " rep ", rep_id, "/", n_reps, " done\n", sep = "")
  }
}

result_rows <- do.call(rbind, rows)
result_rows$has_error <- !is.na(result_rows$error)

summary_rows <- stats::aggregate(
  cbind(n_de, true_de, fdr, power, elapsed) ~ simulation + method + n_per_group,
  data = result_rows,
  FUN = function(v) mean(v, na.rm = TRUE)
)
errors <- stats::aggregate(
  has_error ~ simulation + method + n_per_group,
  data = result_rows,
  FUN = sum
)
names(errors)[4] <- "n_errors"
summary_rows <- merge(summary_rows, errors, by = c("simulation", "method", "n_per_group"), all.x = TRUE)

paper_table2 <- data.frame(
  simulation = rep(c("parametric", "nonparametric"), each = 16),
  method = rep(rep(c("GaGa", "MiGaGa2", "Ga", "limma_BH"), each = 4), 2),
  n_per_group = rep(sample_sizes, 8),
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
summary_rows <- merge(summary_rows, paper_table2, by = c("simulation", "method", "n_per_group"), all.x = TRUE)
summary_rows <- summary_rows[order(summary_rows$simulation, summary_rows$method, summary_rows$n_per_group), ]

roc_all <- do.call(rbind, roc_rows)
roc_avg <- stats::aggregate(power ~ method + fdr_grid, data = roc_all, FUN = mean)
roc_avg <- roc_avg[order(roc_avg$method, roc_avg$fdr_grid), ]
roc_summary <- do.call(rbind, lapply(split(roc_avg, roc_avg$method), function(z) {
  data.frame(method = z$method[1], auc = auc_grid(z$fdr_grid, z$power), stringsAsFactors = FALSE)
}))
rownames(roc_summary) <- NULL

suffix <- paste(data_choice, transform_choice, paste0("reps", n_reps), sep = "_")
write.csv(result_rows, result_file(paste0("armstrong62_simulation_rows_", suffix, ".csv")), row.names = FALSE)
write.csv(summary_rows, result_file(paste0("armstrong62_table2_summary_", suffix, ".csv")), row.names = FALSE)
write.csv(roc_avg, result_file(paste0("armstrong62_roc_curves_", suffix, ".csv")), row.names = FALSE)
write.csv(roc_summary, result_file(paste0("armstrong62_roc_summary_", suffix, ".csv")), row.names = FALSE)
write.csv(data.frame(
  data = data_choice,
  transform = transform_choice,
  n_genes = nrow(x),
  n_samples = ncol(x),
  n_all = sum(groups == levels(groups)[1]),
  n_mll = sum(groups == levels(groups)[2]),
  full_data_truth_de = sum(full_truth),
  parametric_prior_de_probability = unname(full_par$probpat[min(2, length(full_par$probpat))]),
  parametric_generator = "gaga::simGG using full-data GaGa posterior-mean hyperparameters",
  reps = n_reps,
  fdr_target = fdr_target,
  fit_method = fit_method
), result_file(paste0("armstrong62_data_summary_", suffix, ".csv")), row.names = FALSE)
saveRDS(list(summary = summary_rows, rows = result_rows, roc = roc_avg, full_truth = full_truth, full_fit = full_fit),
        result_file(paste0("armstrong62_simulation_", suffix, ".rds")))

plot_roc <- function(file, device = c("pdf", "png")) {
  device <- match.arg(device)
  if (device == "pdf") {
    grDevices::pdf(file, width = 6.2, height = 5.2)
  } else {
    grDevices::png(file, width = 1400, height = 1100, res = 200)
  }
  cols <- c(GaGa = "#1B9E77", MiGaGa2 = "#D95F02", Ga = "#666666", limma_BH = "#7570B3")
  plot(NA, xlim = c(0, 1), ylim = c(0, 1),
       xlab = "Average FDR", ylab = "Average power",
       main = "Armstrong bootstrap simulations")
  for (method in methods) {
    z <- roc_avg[roc_avg$method == method, ]
    lines(z$fdr_grid, z$power, col = cols[[method]], lwd = 2)
  }
  leg <- merge(data.frame(method = methods), roc_summary, by = "method", all.x = TRUE, sort = FALSE)
  legend("bottomright",
         legend = paste0(leg$method, " AUC=", sprintf("%.3f", leg$auc)),
         col = cols[leg$method], lwd = 2, bty = "n")
  grDevices::dev.off()
}

plot_roc(result_file(paste0("armstrong62_roc_curves_", suffix, ".pdf")), "pdf")
plot_roc(result_file(paste0("armstrong62_roc_curves_", suffix, ".png")), "png")

cat("Table 2-style summary:\n")
print(summary_rows)
cat("ROC summary:\n")
print(roc_summary)
cat("Outputs suffix:", suffix, "\n")
