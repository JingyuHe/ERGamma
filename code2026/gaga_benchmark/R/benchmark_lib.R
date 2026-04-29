script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- args[grepl("^--file=", args)]
  if (length(file_arg) > 0) {
    return(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE))
  }
  normalizePath(getwd(), mustWork = TRUE)
}

benchmark_dir <- function() {
  opt <- getOption("gaga.benchmark.dir")
  if (!is.null(opt)) {
    return(normalizePath(opt, mustWork = TRUE))
  }
  normalizePath(file.path(dirname(script_path()), ".."), mustWork = TRUE)
}

setup_benchmark <- function() {
  bench <- benchmark_dir()
  lib <- file.path(bench, "r-lib")
  legacy_lib <- normalizePath(file.path(bench, "..", "r-lib"), mustWork = FALSE)
  dir.create(file.path(bench, "results"), recursive = TRUE, showWarnings = FALSE)
  dir.create(lib, recursive = TRUE, showWarnings = FALSE)
  .libPaths(unique(c(lib, legacy_lib, .libPaths())))
  bench
}

require_pkgs <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      "Missing packages: ", paste(missing, collapse = ", "),
      ". Run Rscript code2026/gaga_benchmark/R/install_deps.R first.",
      call. = FALSE
    )
  }
}

result_file <- function(name) {
  file.path(benchmark_dir(), "results", name)
}

env_int <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) {
    return(default)
  }
  as.integer(value)
}

safe_div <- function(num, den) {
  ifelse(den == 0, 0, num / den)
}

eval_calls <- function(calls, truth) {
  calls <- as.logical(calls)
  truth <- as.logical(truth)
  tp <- sum(calls & truth, na.rm = TRUE)
  fp <- sum(calls & !truth, na.rm = TRUE)
  fn <- sum(!calls & truth, na.rm = TRUE)
  data.frame(
    n_de = sum(calls, na.rm = TRUE),
    true_de = sum(truth, na.rm = TRUE),
    fdr = safe_div(fp, tp + fp),
    power = safe_div(tp, tp + fn)
  )
}

limma_or_ttest <- function(x, groups, fdr = 0.05) {
  groups <- factor(groups)
  if (requireNamespace("limma", quietly = TRUE)) {
    design <- stats::model.matrix(~ groups)
    fit <- limma::lmFit(log2(x), design)
    fit <- limma::eBayes(fit)
    p <- fit$p.value[, 2]
  } else {
    idx <- split(seq_along(groups), groups)
    if (length(idx) != 2) {
      stop("The base t-test fallback only supports two groups.")
    }
    p <- apply(log2(x), 1, function(row) {
      stats::t.test(row[idx[[1]]], row[idx[[2]]])$p.value
    })
  }
  adj <- stats::p.adjust(p, method = "BH")
  list(p = p, adj_p = adj, calls = adj <= fdr, score = -log10(pmax(p, 1e-300)))
}

gas_log_kernel <- function(x, a, b, c, d, r, s) {
  out <- rep(-Inf, length(x))
  ok <- is.finite(x) & x > 0
  xo <- x[ok]
  out[ok] <- lgamma(a * xo + d) - a * lgamma(xo) +
    (a * xo + d) * (log(xo) - log(r + s * xo)) +
    (b - d - 1) * log(xo) - c * xo
  out
}

gas_approx_shape_rate <- function(a, b, c, d, r, s) {
  shape <- a / 2 + b - 1 / 2
  rate <- c + a * log(s / a)
  c(shape = shape, rate = rate)
}

gas_find_upper <- function(a, b, c, d, r, s) {
  ar <- gas_approx_shape_rate(a, b, c, d, r, s)
  if (is.finite(ar["shape"]) && is.finite(ar["rate"]) &&
      ar["shape"] > 0 && ar["rate"] > 0) {
    upper <- max(20, stats::qgamma(0.999999, shape = ar["shape"], rate = ar["rate"]) * 4)
  } else {
    upper <- 50
  }
  lower <- 1e-8
  for (i in seq_len(8)) {
    opt <- stats::optimize(
      function(x) -gas_log_kernel(x, a, b, c, d, r, s),
      interval = c(lower, upper)
    )
    log_mode <- gas_log_kernel(opt$minimum, a, b, c, d, r, s)
    log_upper <- gas_log_kernel(upper, a, b, c, d, r, s)
    if (is.finite(log_mode) && (!is.finite(log_upper) || log_upper < log_mode - 35)) {
      return(upper)
    }
    upper <- upper * 2
  }
  upper
}

gas_exact_moments <- function(a, b, c, d, r, s, rel.tol = 1e-8) {
  upper <- gas_find_upper(a, b, c, d, r, s)
  z_lower <- log(1e-12)
  z_upper <- log(upper)
  log_integrand <- function(z, pow = 0) {
    x <- exp(z)
    gas_log_kernel(x, a, b, c, d, r, s) + (pow + 1) * z
  }
  log_moment <- function(pow) {
    opt <- stats::optimize(
      function(z) -log_integrand(z, pow),
      interval = c(z_lower, z_upper)
    )
    log_scale <- log_integrand(opt$minimum, pow)
    val <- stats::integrate(
      function(z) exp(log_integrand(z, pow) - log_scale),
      z_lower,
      z_upper,
      subdivisions = 1000,
      rel.tol = rel.tol
    )$value
    log(val) + log_scale
  }
  log_z0 <- log_moment(0)
  log_z1 <- log_moment(1)
  log_z2 <- log_moment(2)
  mean <- exp(log_z1 - log_z0)
  second <- exp(log_z2 - log_z0)
  var <- second - mean^2
  if (var < 0 && abs(var) < 1e-8 * max(1, second)) {
    var <- 0
  }
  mode <- exp(stats::optimize(
    function(z) -gas_log_kernel(exp(z), a, b, c, d, r, s),
    interval = c(z_lower, z_upper)
  )$minimum)
  list(mean = mean, var = var, norm = exp(log_z0), mode = mode, upper = upper)
}

rPIG_ERGamma <- function(n, cc, N) {
  require_pkgs("GIGrvg")
  output <- rep(0, n)
  if (N > 0) {
    for (i in seq_len(N)) {
      if (cc == 0) {
        output <- output + 1 / stats::rgamma(n, shape = 3 / 2) / (4 * i^2)
      } else {
        output <- output + GIGrvg::rgig(
          n = n,
          lambda = -3 / 2,
          chi = 1 / (2 * i^2),
          psi = 2 * cc^2
        )
      }
    }
  }
  if (cc == 0) {
    temp1 <- trigamma(1 + N)
    temp2 <- psigamma(1 + N, deriv = 2)
    rate <- -4 * temp1 / temp2
    shape <- 0.5 * temp1 * rate
  } else {
    temp1 <- digamma(1 + N + cc) - digamma(1 + N)
    temp2 <- temp1 - cc * trigamma(1 + N + cc)
    rate <- 2 * cc^2 * temp1 / temp2
    shape <- rate / (2 * cc) * temp1
  }
  output + stats::rgamma(n, shape = shape, rate = rate)
}

rH <- function(m, a, b, max_try = 1000) {
  tau <- if (b > 0) {
    0.5 + sqrt(0.25 + 2 * a * m / b^2)
  } else {
    -0.5 + sqrt(0.25 + 2 * a * m / b^2)
  }
  v1 <- tau * abs(b) - b
  v2 <- tau * abs(b) / (2 * a)
  for (i in seq_len(max_try)) {
    x <- stats::rgamma(1, m, rate = v1)
    if (stats::runif(1) <= exp(-a * (x - v2)^2)) {
      return(x)
    }
  }
  NA_real_
}

sample_pig_shape_special <- function(delta, mu, n, burnin, trunc, init = 5) {
  gamma_const <- -digamma(1)
  b_pre <- delta * gamma_const - delta * log(delta * mu)
  xi <- init
  out <- numeric(n)
  for (iter in seq_len(n + burnin)) {
    tau <- stats::rgamma(1, shape = delta * xi + 1)
    omega <- rPIG_ERGamma(delta, xi, trunc)
    b <- b_pre + delta * log(tau)
    xi <- rH(delta + 1, sum(omega), b)
    if (!is.finite(xi)) {
      stop("rH failed to accept a draw.")
    }
    if (iter > burnin) {
      out[iter - burnin] <- xi
    }
  }
  out
}

special_log_kernel <- function(x, delta, mu) {
  out <- rep(-Inf, length(x))
  ok <- is.finite(x) & x > 0
  xo <- x[ok]
  out[ok] <- lgamma(delta * xo + 1) -
    delta * lgamma(xo) -
    delta * xo * log(delta * mu)
  out
}

special_exact_moments <- function(delta, mu) {
  upper <- if (delta >= 30) 12 else 30
  opt <- stats::optimize(
    function(x) -special_log_kernel(x, delta, mu),
    interval = c(1e-8, upper)
  )
  log_mode <- special_log_kernel(opt$minimum, delta, mu)
  scaled <- function(x, pow = 0) {
    exp(special_log_kernel(x, delta, mu) - log_mode) * x^pow
  }
  z0 <- stats::integrate(function(x) scaled(x, 0), 0, upper,
                         subdivisions = 1000, rel.tol = 1e-8)$value
  z1 <- stats::integrate(function(x) scaled(x, 1), 0, upper,
                         subdivisions = 1000, rel.tol = 1e-8)$value
  z2 <- stats::integrate(function(x) scaled(x, 2), 0, upper,
                         subdivisions = 1000, rel.tol = 1e-8)$value
  mean <- z1 / z0
  var <- z2 / z0 - mean^2
  list(mean = mean, var = var, norm = z0 * exp(log_mode), mode = opt$minimum)
}

calc_AB_miller <- function(delta, mu, init.x = 1, tol = 1e-8, max.T = 500) {
  x <- init.x
  for (i in seq_len(max.T)) {
    A <- 1 - x^2 * delta * (delta * trigamma(delta * x + 1) - trigamma(x))
    B <- (A - 1) / x -
      delta * (digamma(delta * x + 1) - digamma(x) - log(mu * delta))
    x_new <- A / B
    if (is.finite(x_new) && abs(x / x_new - 1) < tol) {
      return(c(shape = A, rate = B))
    }
    x <- x_new
  }
  c(shape = A, rate = B)
}
