args <- commandArgs(FALSE)
this_file <- sub("^--file=", "", args[grepl("^--file=", args)][1])
source(file.path(dirname(normalizePath(this_file, mustWork = TRUE)), "benchmark_lib.R"))

bench <- setup_benchmark()
paths <- .libPaths()
system_paths <- paths[grepl("^/opt/homebrew", paths)]
local_paths <- paths[!grepl("^/opt/homebrew", paths)]
.libPaths(unique(c(system_paths, file.path(bench, "r-lib-gcrma"), local_paths)))
require_pkgs(c("gaga", "EBarrays", "Biobase"))

data_choice <- Sys.getenv("FIG1_DATA", unset = "schliep_filtered")
transform_choice <- Sys.getenv("FIG1_TRANSFORM", unset = if (data_choice == "schliep_filtered") "log" else "log2_floor20")
fit_method <- Sys.getenv("FIG1_GAGA_METHOD", unset = "quickEM")
n_prior <- env_int("FIG1_N_PRIOR", 10000)
seed <- env_int("FIG1_SEED", 109)

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
  stop("Unknown FIG1_TRANSFORM: ", transform)
}

gaga_patterns <- function(groups) {
  patterns <- matrix(c(0, 0, 0, 1), 2, 2)
  colnames(patterns) <- levels(factor(groups))
  patterns
}

eb_patterns <- function(groups) {
  lv <- levels(factor(groups))
  EBarrays::ebPatterns(c(
    paste(rep(1, length(groups)), collapse = " "),
    paste(ifelse(groups == lv[1], 1, 2), collapse = " ")
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
  gaga::parest(fit, x = x, groups = groups, alpha = 0.05)
}

fit_ga_calls <- function(x, groups) {
  fit <- EBarrays::emfit(
    data = x,
    family = "GG",
    hypotheses = eb_patterns(groups),
    num.iter = 20,
    verbose = FALSE
  )
  post <- EBarrays::postprob(fit, x)$pattern
  threshold <- EBarrays::crit.fun(post[, 1], 0.05)
  post[, 2] > threshold
}

prior_values <- function(fit, n, m) {
  par <- gaga::getpar(fit)
  sim <- gaga::simGG(
    n = n,
    m = m,
    p.de = 0,
    a0 = par$a0,
    nu = par$nu,
    balpha = par$balpha,
    nualpha = par$nualpha,
    equalcv = TRUE,
    probclus = par$probclus
  )
  as.vector(Biobase::exprs(sim))
}

if (data_choice == "schliep_filtered") {
  dat <- load_armstrong_schliep(file.path(bench, "data", "armstrong-2002-v2_database.txt"))
} else if (data_choice == "orange_full") {
  dat <- load_armstrong_orange(file.path(bench, "data", "MLL_orange_12533.tab"))
} else {
  stop("Unknown FIG1_DATA: ", data_choice)
}

set.seed(seed)
x <- transform_expr(dat$x, transform_choice)
groups <- dat$groups
if (any(!is.finite(x)) || min(x) < 0) {
  stop("Transformed data must be finite and non-negative for gamma-family fits.")
}

cat("Armstrong Figure 1 data:\n")
cat("  data:", data_choice, "\n")
cat("  transform:", transform_choice, "\n")
cat("  dimensions:", nrow(x), "genes x", ncol(x), "samples\n")
cat("  groups:", paste(names(table(groups)), as.integer(table(groups)), collapse = ", "), "\n")

cat("Fitting Ga, GaGa and MiGaGa2...\n")
ga_calls <- fit_ga_calls(x, groups)
gaga_fit <- fit_gaga(x, groups, nclust = 1)
migaga_fit <- fit_gaga(x, groups, nclust = 2)

idx <- split(seq_along(groups), groups)
gene_mean <- rowMeans(x)
vars <- lapply(idx, function(cols) apply(x[, cols, drop = FALSE], 1, stats::var))
ns <- vapply(idx, length, integer(1))
pooled_sd <- sqrt(Reduce(`+`, Map(function(v, n) (n - 1) * v, vars, ns)) / (sum(ns) - length(ns)))
gene_cv <- pooled_sd / pmax(gene_mean, .Machine$double.eps)

v_obs <- as.vector(x)
v_gaga <- prior_values(gaga_fit, n_prior, ncol(x))
v_migaga <- prior_values(migaga_fit, n_prior, ncol(x))
obs_max <- max(v_obs, na.rm = TRUE) * 1.05
v_gaga <- v_gaga[is.finite(v_gaga) & v_gaga >= 0 & v_gaga <= obs_max]
v_migaga <- v_migaga[is.finite(v_migaga) & v_migaga >= 0 & v_migaga <= obs_max]

plot_figure1 <- function(file, device = c("pdf", "png")) {
  device <- match.arg(device)
  if (device == "pdf") {
    grDevices::pdf(file, width = 10, height = 5)
  } else {
    grDevices::png(file, width = 2000, height = 1000, res = 200)
  }
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(oldpar)
    grDevices::dev.off()
  })
  graphics::par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))

  graphics::plot(
    gene_mean, gene_cv,
    pch = 20,
    cex = 0.25,
    col = grDevices::rgb(0, 0, 0, 0.35),
    xlab = "Mean expression",
    ylab = "Coefficient of variation",
    main = "(a) Mean and CV"
  )
  graphics::points(gene_mean[ga_calls], gene_cv[ga_calls], pch = 8, cex = 0.65)

  xrng <- c(0, max(v_obs, na.rm = TRUE) + 0.5)
  d_obs <- stats::density(v_obs, from = xrng[1], to = xrng[2], n = 512)
  d_gaga <- stats::density(v_gaga, from = xrng[1], to = xrng[2], n = 512)
  d_migaga <- stats::density(v_migaga, from = xrng[1], to = xrng[2], n = 512)
  ymax <- max(d_obs$y, d_gaga$y, d_migaga$y) * 1.05
  graphics::plot(
    d_obs,
    lwd = 2,
    col = "black",
    xlim = xrng,
    ylim = c(0, ymax),
    xlab = "Expression levels",
    ylab = "Density",
    main = "(b) Marginal density"
  )
  graphics::lines(d_gaga$x, d_gaga$y, lwd = 2, lty = 2)
  graphics::lines(d_migaga$x, d_migaga$y, lwd = 2, lty = 3)
  graphics::legend(
    "topright",
    legend = c("Observed data", "GaGa", "MiGaGa2"),
    lty = c(1, 2, 3),
    lwd = 2,
    bty = "n"
  )
}

suffix <- paste(data_choice, transform_choice, sep = "_")
plot_figure1(result_file(paste0("armstrong_figure1_reproduction_", suffix, ".pdf")), "pdf")
plot_figure1(result_file(paste0("armstrong_figure1_reproduction_", suffix, ".png")), "png")

summary <- data.frame(
  data = data_choice,
  transform = transform_choice,
  n_genes = nrow(x),
  n_samples = ncol(x),
  n_all = sum(groups == levels(groups)[1]),
  n_mll = sum(groups == levels(groups)[2]),
  ga_de = sum(ga_calls),
  gaga_probpat_de = unname(gaga::getpar(gaga_fit)$probpat[min(2, length(gaga::getpar(gaga_fit)$probpat))]),
  migaga_probpat_de = unname(gaga::getpar(migaga_fit)$probpat[min(2, length(gaga::getpar(migaga_fit)$probpat))]),
  fit_method = fit_method,
  n_prior = n_prior
)
write.csv(summary, result_file(paste0("armstrong_figure1_summary_", suffix, ".csv")), row.names = FALSE)
saveRDS(
  list(summary = summary, ga_calls = ga_calls, gaga_fit = gaga_fit, migaga_fit = migaga_fit),
  result_file(paste0("armstrong_figure1_data_", suffix, ".rds"))
)

cat("Figure 1 summary:\n")
print(summary)
cat("Outputs suffix:", suffix, "\n")
