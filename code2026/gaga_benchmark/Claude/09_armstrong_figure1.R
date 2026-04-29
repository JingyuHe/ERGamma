#!/usr/bin/env Rscript
# =============================================================================
# 09_armstrong_figure1.R  (Claude/ debugged version)
# Reproduces Rossell (2009) Figure 1: per-gene mean vs CV (panel a) with stars
# on Ga-flagged DE genes; observed marginal density vs GaGa / MiGaGa2 prior
# predictive (panel b). Same logic as the legacy R/09; just uses Claude/
# benchmark_lib (so the Mach-O library issue is gone).
# =============================================================================
args <- commandArgs(FALSE)
this_file <- sub("^--file=", "", args[grepl("^--file=", args)][1])
source(file.path(dirname(normalizePath(this_file, mustWork = TRUE)),
                 "benchmark_lib.R"))
bench <- setup_benchmark()
require_pkgs(c("gaga", "EBarrays", "Biobase"))

data_choice <- Sys.getenv("FIG1_DATA", unset = "schliep_filtered")
transform_choice <- Sys.getenv(
  "FIG1_TRANSFORM",
  unset = if (data_choice == "schliep_filtered") "log" else "log2_floor20"
)
fit_method <- Sys.getenv("FIG1_GAGA_METHOD", unset = "quickEM")
n_prior <- env_int("FIG1_N_PRIOR", 10000)
seed <- env_int("FIG1_SEED", 109)
equalcv <- Sys.getenv("FIG1_EQUALCV", unset = "1") != "0"

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
restrict_to_all_mll <- function(d) {
  keep <- d$groups %in% c("ALL", "MLL")
  drop_idx <- cdf_safe_drop_u95av2(d$groups, d$cdf)
  if (length(drop_idx)) keep[drop_idx] <- FALSE
  list(x = d$x[, keep, drop = FALSE],
       groups = factor(d$groups[keep], levels = c("ALL", "MLL")))
}
apply_transform <- function(x_raw, t) {
  if (t == "log") log(pmax(x_raw, 1))
  else if (t == "log2_floor20") log2(pmax(x_raw, 20))
  else if (t == "log2_clip") log2(pmax(x_raw, 1) + 1)
  else stop("unknown transform: ", t)
}

if (data_choice == "schliep_filtered") {
  dat <- load_armstrong_schliep(file.path(bench, "data",
                                          "armstrong-2002-v2_database.txt"))
} else {
  dat <- load_armstrong_orange(file.path(bench, "data",
                                         "MLL_orange_12533.tab"))
}
dat <- restrict_to_all_mll(dat)

set.seed(seed)
x <- apply_transform(dat$x, transform_choice)
rownames(x) <- rownames(dat$x); groups <- dat$groups
if (any(!is.finite(x)) || min(x) < 0)
  stop("Transformed data must be finite and non-negative.")

cat("Figure 1 data:", data_choice, transform_choice,
    nrow(x), "x", ncol(x), "\n")

eb_pat <- function(g) {
  lv <- levels(factor(g))
  EBarrays::ebPatterns(c(
    paste(rep(1, length(g)), collapse = " "),
    paste(ifelse(g == lv[1], 1, 2), collapse = " ")))
}
ga_calls <- {
  fit <- EBarrays::emfit(data = pmax(dat$x, .Machine$double.eps),
                         family = "GG", hypotheses = eb_pat(groups),
                         num.iter = 20, verbose = FALSE)
  post <- EBarrays::postprob(fit, pmax(dat$x, .Machine$double.eps))$pattern
  thr <- EBarrays::crit.fun(post[, 1], 0.05)
  post[, 2] > thr
}
cat("Ga DE calls:", sum(ga_calls), "\n")

fit_gaga <- function(nclust) {
  patterns <- matrix(c(0, 0, 0, 1), 2, 2)
  colnames(patterns) <- levels(factor(groups))
  fit <- gaga::fitGG(x, groups, patterns = patterns, equalcv = equalcv,
                     nclust = nclust, method = fit_method, trace = FALSE)
  gaga::parest(fit, x = x, groups = groups, alpha = 0.05)
}
gaga_fit <- fit_gaga(1); migaga_fit <- fit_gaga(2)

prior_values <- function(fit, n, m) {
  par <- gaga::getpar(fit)
  sim <- gaga::simGG(n = n, m = m, p.de = 0,
                     a0 = par$a0, nu = par$nu,
                     balpha = par$balpha, nualpha = par$nualpha,
                     equalcv = equalcv, probclus = par$probclus)
  as.vector(Biobase::exprs(sim))
}
gene_mean <- rowMeans(x)
idx <- split(seq_along(groups), groups)
vars <- lapply(idx, function(c) apply(x[, c, drop = FALSE], 1, stats::var))
ns <- vapply(idx, length, integer(1))
pooled_sd <- sqrt(Reduce(`+`, Map(function(v, n) (n - 1) * v, vars, ns)) /
                  (sum(ns) - length(ns)))
gene_cv <- pooled_sd / pmax(gene_mean, .Machine$double.eps)

v_obs <- as.vector(x)
v_gaga <- prior_values(gaga_fit, n_prior, ncol(x))
v_migaga <- prior_values(migaga_fit, n_prior, ncol(x))
obs_max <- max(v_obs) * 1.05
v_gaga <- v_gaga[is.finite(v_gaga) & v_gaga >= 0 & v_gaga <= obs_max]
v_migaga <- v_migaga[is.finite(v_migaga) & v_migaga >= 0 &
                     v_migaga <= obs_max]

plot_figure1 <- function(file, dev = c("pdf", "png")) {
  dev <- match.arg(dev)
  if (dev == "pdf") grDevices::pdf(file, width = 10, height = 5)
  else grDevices::png(file, width = 2000, height = 1000, res = 200)
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit({ graphics::par(oldpar); grDevices::dev.off() })
  graphics::par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
  graphics::plot(gene_mean, gene_cv, pch = 20, cex = 0.25,
                 col = grDevices::rgb(0, 0, 0, 0.35),
                 xlab = "Mean expression", ylab = "Coefficient of variation",
                 main = "(a) Mean and CV")
  graphics::points(gene_mean[ga_calls], gene_cv[ga_calls],
                   pch = 8, cex = 0.65)
  xrng <- c(0, max(v_obs, na.rm = TRUE) + 0.5)
  d_obs <- stats::density(v_obs,    from = xrng[1], to = xrng[2], n = 512)
  d_g   <- stats::density(v_gaga,   from = xrng[1], to = xrng[2], n = 512)
  d_m   <- stats::density(v_migaga, from = xrng[1], to = xrng[2], n = 512)
  ymax <- max(d_obs$y, d_g$y, d_m$y) * 1.05
  graphics::plot(d_obs, lwd = 2, col = "black",
                 xlim = xrng, ylim = c(0, ymax),
                 xlab = "Expression levels", ylab = "Density",
                 main = "(b) Marginal density")
  graphics::lines(d_g$x, d_g$y, lwd = 2, lty = 2)
  graphics::lines(d_m$x, d_m$y, lwd = 2, lty = 3)
  graphics::legend("topright",
                   legend = c("Observed", "GaGa", "MiGaGa2"),
                   lty = c(1, 2, 3), lwd = 2, bty = "n")
}
suffix <- paste(data_choice, transform_choice, sep = "_")
plot_figure1(result_file(paste0("armstrong_figure1_reproduction_", suffix,
                                ".pdf")), "pdf")
plot_figure1(result_file(paste0("armstrong_figure1_reproduction_", suffix,
                                ".png")), "png")

summary <- data.frame(
  data = data_choice, transform = transform_choice,
  n_genes = nrow(x), n_samples = ncol(x),
  n_all = sum(groups == "ALL"), n_mll = sum(groups == "MLL"),
  ga_de = sum(ga_calls),
  gaga_probpat_de = unname(gaga::getpar(gaga_fit)$probpat[
    min(2, length(gaga::getpar(gaga_fit)$probpat))]),
  migaga_probpat_de = unname(gaga::getpar(migaga_fit)$probpat[
    min(2, length(gaga::getpar(migaga_fit)$probpat))]),
  fit_method = fit_method, n_prior = n_prior)
write.csv(summary,
  result_file(paste0("armstrong_figure1_summary_", suffix, ".csv")),
  row.names = FALSE)
print(summary)
cat("Outputs suffix:", suffix, "\n")
