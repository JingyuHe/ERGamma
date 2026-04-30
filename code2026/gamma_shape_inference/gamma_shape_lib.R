gamma_shape_find_repo_root <- function(start = getwd()) {
  path <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (dir.exists(file.path(path, "code2026", "gamma_shape_inference"))) {
      return(path)
    }
    parent <- dirname(path)
    if (identical(parent, path)) {
      stop("Could not find repository root from ", start, call. = FALSE)
    }
    path <- parent
  }
}

gamma_shape_init_libpaths <- function(root = gamma_shape_find_repo_root()) {
  local_lib <- file.path(root, "code2026", "r-lib")
  if (dir.exists(local_lib)) {
    .libPaths(unique(c(normalizePath(local_lib, mustWork = TRUE), .libPaths())))
  }
  invisible(.libPaths())
}

require_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Required package not available: ", pkg, call. = FALSE)
  }
}

euler_gamma_const <- function() {
  -digamma(1)
}

is_integerish <- function(x, tol = 1e-10) {
  is.finite(x) && abs(x - round(x)) < tol
}

log_xi2 <- function(alpha, delta, mu_ratio) {
  out <- rep(-Inf, length(alpha))
  ok <- alpha > 0 & is.finite(alpha)
  if (!is.finite(delta) || delta <= 0 || !is.finite(mu_ratio) || mu_ratio <= 1) {
    return(out)
  }
  a <- alpha[ok]
  out[ok] <- lgamma(delta * a + 1) -
    delta * lgamma(a) -
    delta * a * log(delta * mu_ratio)
  out
}

log_xi2_z <- function(z, delta, mu_ratio) {
  alpha <- exp(z)
  log_xi2(alpha, delta, mu_ratio) + z
}

stirling_gamma_params <- function(delta, mu_ratio) {
  rate <- delta * log(mu_ratio)
  shape <- (delta + 3) / 2
  if (!is.finite(rate) || rate <= 0 || !is.finite(shape) || shape <= 0) {
    return(list(ok = FALSE, shape = NA_real_, rate = NA_real_))
  }
  list(ok = TRUE, shape = shape, rate = rate)
}

xi2_mode_start <- function(delta, mu_ratio) {
  st <- stirling_gamma_params(delta, mu_ratio)
  if (isTRUE(st$ok)) {
    return(max(1e-8, (st$shape - 1) / st$rate))
  }
  1
}

xi2_grid_truth <- function(delta,
                           mu_ratio,
                           grid_size = 20001,
                           tail_nats = 50,
                           min_z = -50,
                           max_z = 50) {
  if (delta <= 0 || mu_ratio <= 1) {
    stop("xi2_grid_truth requires delta > 0 and mu_ratio > 1", call. = FALSE)
  }

  z_center <- log(xi2_mode_start(delta, mu_ratio))
  lower <- max(min_z, z_center - 10)
  upper <- min(max_z, z_center + 10)

  opt <- optimize(
    f = function(z) -log_xi2_z(z, delta, mu_ratio),
    interval = c(lower, upper)
  )
  z_mode <- opt$minimum
  log_mode <- -opt$objective

  left <- z_mode
  while (left > min_z && log_xi2_z(left, delta, mu_ratio) > log_mode - tail_nats) {
    left <- left - 1
  }
  right <- z_mode
  while (right < max_z && log_xi2_z(right, delta, mu_ratio) > log_mode - tail_nats) {
    right <- right + 1
  }

  z <- seq(left, right, length.out = grid_size)
  log_kernel_z <- log_xi2_z(z, delta, mu_ratio)
  dz <- z[2] - z[1]
  trap_w <- rep(dz, length(z))
  trap_w[1] <- dz / 2
  trap_w[length(trap_w)] <- dz / 2

  log_shift <- max(log_kernel_z)
  mass <- exp(log_kernel_z - log_shift) * trap_w
  total_mass <- sum(mass)
  prob <- mass / total_mass
  cdf <- cumsum(prob)
  cdf[length(cdf)] <- 1

  alpha <- exp(z)
  mean_alpha <- sum(prob * alpha)
  var_alpha <- sum(prob * (alpha - mean_alpha)^2)
  sd_alpha <- sqrt(var_alpha)
  skew <- sum(prob * (alpha - mean_alpha)^3) / sd_alpha^3
  kurt <- sum(prob * (alpha - mean_alpha)^4) / var_alpha^2
  qs <- approx(cdf, alpha, xout = c(0.025, 0.5, 0.975), ties = "ordered")$y

  log_norm <- log(total_mass) + log_shift
  density_alpha <- exp(log_xi2(alpha, delta, mu_ratio) - log_norm)

  list(
    delta = delta,
    mu_ratio = mu_ratio,
    z = z,
    alpha = alpha,
    density = density_alpha,
    cdf = cdf,
    prob = prob,
    mode = exp(z_mode),
    mean = mean_alpha,
    var = var_alpha,
    skewness = skew,
    kurtosis = kurt,
    q025 = qs[1],
    q50 = qs[2],
    q975 = qs[3],
    log_norm = log_norm
  )
}

truth_cdf_at <- function(truth, x) {
  approx(truth$alpha, truth$cdf, xout = x, yleft = 0, yright = 1, ties = "ordered")$y
}

truth_quantile_at <- function(truth, p) {
  approx(truth$cdf, truth$alpha, xout = p, ties = "ordered")$y
}

trap_integral <- function(x, y) {
  if (length(x) < 2) {
    return(NA_real_)
  }
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
}

ks_sample <- function(samples, truth) {
  x <- sort(samples[is.finite(samples)])
  n <- length(x)
  if (n < 1) {
    return(NA_real_)
  }
  f0 <- truth_cdf_at(truth, x)
  max(abs(f0 - ((seq_len(n) - 0.5) / n)))
}

w1_sample <- function(samples, truth) {
  x <- sort(samples[is.finite(samples)])
  n <- length(x)
  if (n < 1) {
    return(NA_real_)
  }
  q <- truth_quantile_at(truth, (seq_len(n) - 0.5) / n)
  mean(abs(x - q))
}

acf_lag <- function(samples, lag) {
  samples <- samples[is.finite(samples)]
  if (length(samples) <= lag + 1 || stats::sd(samples) == 0) {
    return(NA_real_)
  }
  stats::acf(samples, lag.max = lag, plot = FALSE)$acf[lag + 1]
}

ess_value <- function(samples) {
  samples <- samples[is.finite(samples)]
  if (length(samples) < 3 || stats::sd(samples) == 0) {
    return(NA_real_)
  }
  if (requireNamespace("coda", quietly = TRUE)) {
    return(as.numeric(coda::effectiveSize(samples)))
  }
  n <- length(samples)
  rho <- stats::acf(samples, plot = FALSE, lag.max = min(1000, n - 2))$acf[-1]
  positive <- rho[rho > 0]
  n / max(1, 1 + 2 * sum(positive))
}

sample_summary_row <- function(samples,
                               truth,
                               method,
                               runtime_sec,
                               n_iter,
                               alpha_true = NA_real_,
                               accept_rate = NA_real_,
                               avg_evals = NA_real_,
                               pig_N = NA_real_,
                               status = "ok",
                               message = "") {
  n_save <- sum(is.finite(samples))
  ci <- if (n_save > 0) {
    as.numeric(stats::quantile(samples, c(0.025, 0.975), names = FALSE, type = 8))
  } else {
    c(NA_real_, NA_real_)
  }
  ess <- ess_value(samples)
  data.frame(
    method = method,
    status = status,
    ks = ks_sample(samples, truth),
    w1 = w1_sample(samples, truth),
    mean = if (n_save > 0) mean(samples, na.rm = TRUE) else NA_real_,
    var = if (n_save > 1) stats::var(samples, na.rm = TRUE) else NA_real_,
    ci_low = ci[1],
    ci_high = ci[2],
    cover = if (is.finite(alpha_true)) ci[1] <= alpha_true && alpha_true <= ci[2] else NA,
    ess = ess,
    ess_per_iter = ess / n_iter,
    ess_per_sec = ess / runtime_sec,
    acf1 = acf_lag(samples, 1),
    acf10 = acf_lag(samples, 10),
    runtime_sec = runtime_sec,
    accept_rate = accept_rate,
    avg_evals = avg_evals,
    pig_N = pig_N,
    n_save = n_save,
    message = message,
    stringsAsFactors = FALSE
  )
}

miller_gamma_params <- function(delta,
                                mu_ratio,
                                init = xi2_mode_start(delta, mu_ratio),
                                tol = 1e-10,
                                max_iter = 500) {
  starts <- unique(pmax(1e-6, c(init, xi2_mode_start(delta, mu_ratio), 0.05, 0.2, 1, 5, 20)))
  last <- list(A = NA_real_, B = NA_real_, x = NA_real_)
  for (x0 in starts) {
    x <- x0
    for (i in seq_len(max_iter)) {
      l1 <- delta * (digamma(delta * x + 1) - digamma(x) - log(delta * mu_ratio))
      l2 <- delta^2 * trigamma(delta * x + 1) - delta * trigamma(x)
      A <- 1 - x^2 * l2
      B <- (A - 1) / x - l1
      last <- list(A = A, B = B, x = x)
      if (!is.finite(A) || !is.finite(B) || A <= 0 || B <= 0) {
        break
      }
      x_new <- A / B
      if (!is.finite(x_new) || x_new <= 0) {
        break
      }
      if (abs(x / x_new - 1) < tol) {
        return(list(ok = TRUE, shape = A, rate = B, iter = i, init = x0))
      }
      x <- x_new
    }
  }
  list(ok = FALSE, shape = last$A, rate = last$B, iter = NA_integer_, init = NA_real_)
}

gamma_approx_row <- function(truth,
                             shape,
                             rate,
                             method,
                             runtime_sec,
                             alpha_true = NA_real_,
                             status = "ok",
                             message = "") {
  if (!is.finite(shape) || !is.finite(rate) || shape <= 0 || rate <= 0) {
    return(data.frame(
      method = method,
      status = "failed",
      ks = NA_real_, w1 = NA_real_, mean = NA_real_, var = NA_real_,
      ci_low = NA_real_, ci_high = NA_real_, cover = NA,
      ess = NA_real_, ess_per_iter = NA_real_, ess_per_sec = NA_real_,
      acf1 = NA_real_, acf10 = NA_real_, runtime_sec = runtime_sec,
      accept_rate = NA_real_, avg_evals = NA_real_, pig_N = NA_real_,
      n_save = NA_integer_, message = message,
      stringsAsFactors = FALSE
    ))
  }

  approx_cdf <- stats::pgamma(truth$alpha, shape = shape, rate = rate)
  ks <- max(abs(approx_cdf - truth$cdf), na.rm = TRUE)
  w1 <- trap_integral(truth$alpha, abs(approx_cdf - truth$cdf))
  ci <- stats::qgamma(c(0.025, 0.975), shape = shape, rate = rate)
  data.frame(
    method = method,
    status = status,
    ks = ks,
    w1 = w1,
    mean = shape / rate,
    var = shape / rate^2,
    ci_low = ci[1],
    ci_high = ci[2],
    cover = if (is.finite(alpha_true)) ci[1] <= alpha_true && alpha_true <= ci[2] else NA,
    ess = NA_real_,
    ess_per_iter = NA_real_,
    ess_per_sec = NA_real_,
    acf1 = NA_real_,
    acf10 = NA_real_,
    runtime_sec = runtime_sec,
    accept_rate = NA_real_,
    avg_evals = NA_real_,
    pig_N = NA_real_,
    n_save = NA_integer_,
    message = message,
    stringsAsFactors = FALSE
  )
}

sample_slice_xi2 <- function(delta,
                             mu_ratio,
                             n_iter,
                             burnin,
                             thin = 1,
                             init = xi2_mode_start(delta, mu_ratio),
                             w = 1,
                             m = 100) {
  logf <- function(z) log_xi2_z(z, delta, mu_ratio)
  z <- log(max(init, 1e-12))
  n_save <- floor((n_iter - burnin) / thin)
  samples <- rep(NA_real_, n_save)
  evals <- 0
  saved <- 0

  for (iter in seq_len(n_iter)) {
    fz <- logf(z)
    evals <- evals + 1
    log_y <- fz - stats::rexp(1)
    u <- stats::runif(1)
    L <- z - w * u
    R <- L + w
    J <- floor(m * stats::runif(1))
    K <- (m - 1) - J

    while (J > 0 && logf(L) > log_y) {
      evals <- evals + 1
      L <- L - w
      J <- J - 1
    }
    while (K > 0 && logf(R) > log_y) {
      evals <- evals + 1
      R <- R + w
      K <- K - 1
    }

    repeat {
      z_prop <- stats::runif(1, L, R)
      f_prop <- logf(z_prop)
      evals <- evals + 1
      if (f_prop >= log_y) {
        z <- z_prop
        break
      }
      if (z_prop < z) {
        L <- z_prop
      } else {
        R <- z_prop
      }
    }

    if (iter > burnin && ((iter - burnin) %% thin == 0)) {
      saved <- saved + 1
      samples[saved] <- exp(z)
    }
  }
  list(samples = samples, avg_evals = evals / n_iter)
}

sample_miller_imh <- function(delta,
                              mu_ratio,
                              shape,
                              rate,
                              n_iter,
                              burnin,
                              thin = 1,
                              init = xi2_mode_start(delta, mu_ratio)) {
  if (!is.finite(shape) || !is.finite(rate) || shape <= 0 || rate <= 0) {
    stop("Invalid Miller proposal parameters", call. = FALSE)
  }
  log_target <- function(x) log_xi2(x, delta, mu_ratio)
  log_q <- function(x) stats::dgamma(x, shape = shape, rate = rate, log = TRUE)

  x <- max(init, 1e-12)
  if (!is.finite(log_target(x))) {
    x <- stats::rgamma(1, shape = shape, rate = rate)
  }
  n_save <- floor((n_iter - burnin) / thin)
  samples <- rep(NA_real_, n_save)
  saved <- 0
  accept <- 0

  for (iter in seq_len(n_iter)) {
    prop <- stats::rgamma(1, shape = shape, rate = rate)
    log_ratio <- log_target(prop) - log_target(x) + log_q(x) - log_q(prop)
    if (is.finite(log_ratio) && log(stats::runif(1)) < log_ratio) {
      x <- prop
      accept <- accept + 1
    }
    if (iter > burnin && ((iter - burnin) %% thin == 0)) {
      saved <- saved + 1
      samples[saved] <- x
    }
  }
  list(samples = samples, accept_rate = accept / n_iter)
}

rH_slice_fallback <- function(m, a, b, n_steps = 30) {
  if (a <= 0) {
    return(stats::rgamma(1, shape = m, rate = max(1e-8, -b)))
  }
  x_mode <- if (m > 1) {
    (b + sqrt(b^2 + 8 * a * (m - 1))) / (4 * a)
  } else {
    max(1e-8, b / (2 * a))
  }
  z <- log(max(x_mode, 1e-8))
  logf <- function(z0) m * z0 - a * exp(2 * z0) + b * exp(z0)
  for (i in seq_len(n_steps)) {
    log_y <- logf(z) - stats::rexp(1)
    L <- z - stats::runif(1)
    R <- L + 1
    while (logf(L) > log_y) L <- L - 1
    while (logf(R) > log_y) R <- R + 1
    repeat {
      z_prop <- stats::runif(1, L, R)
      if (logf(z_prop) >= log_y) {
        z <- z_prop
        break
      }
      if (z_prop < z) L <- z_prop else R <- z_prop
    }
  }
  exp(z)
}

rH <- function(m, a, b, max_try = 1000) {
  if (!is.finite(m) || !is.finite(a) || !is.finite(b) || m <= 0 || a < 0) {
    return(NA_real_)
  }
  if (a == 0) {
    if (b >= 0) return(NA_real_)
    return(stats::rgamma(1, shape = m, rate = -b))
  }
  if (abs(b) < 1e-10) {
    return(sqrt(stats::rgamma(1, shape = m / 2, rate = a)))
  }

  tau <- if (b > 0) {
    0.5 + sqrt(0.25 + 2 * a * m / b^2)
  } else {
    -0.5 + sqrt(0.25 + 2 * a * m / b^2)
  }
  v1 <- tau * abs(b) - b
  v2 <- tau * abs(b) / (2 * a)
  if (!is.finite(v1) || !is.finite(v2) || v1 <= 0) {
    return(rH_slice_fallback(m, a, b))
  }

  for (i in seq_len(max_try)) {
    x <- stats::rgamma(1, m, rate = v1)
    if (stats::runif(1) <= exp(-a * (x - v2)^2)) {
      return(x)
    }
  }
  rH_slice_fallback(m, a, b)
}

rpig_sum <- function(n, cc, N) {
  require_pkg("GIGrvg")
  n <- as.integer(n)
  if (n <= 0) {
    return(0)
  }
  cc <- abs(cc)

  if (N > 0) {
    k <- seq_len(N)
    if (cc < 1e-12) {
      finite_sum <- sum(1 / stats::rgamma(n * N, shape = 3 / 2) /
                          rep(4 * k^2, each = n))
    } else {
      finite_sum <- 0
      # GIGrvg::rgig() uses scalar chi/psi parameters in this setting; keep
      # the truncation terms explicit so each k gets the intended GIG law.
      for (kk in k) {
        finite_sum <- finite_sum + sum(GIGrvg::rgig(
          n = n,
          lambda = -3 / 2,
          chi = 1 / (2 * kk^2),
          psi = 2 * cc^2
        ))
      }
    }
  } else {
    finite_sum <- 0
  }

  if (cc < 1e-12) {
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
  finite_sum + stats::rgamma(1, shape = n * shape, rate = rate)
}

sample_pig_xi2 <- function(delta,
                           mu_ratio,
                           n_iter,
                           burnin,
                           thin = 1,
                           pig_N = 200,
                           init = xi2_mode_start(delta, mu_ratio),
                           max_delta = Inf) {
  if (!is_integerish(delta)) {
    stop("P-IG Gibbs requires integer posterior delta in this implementation", call. = FALSE)
  }
  delta_int <- as.integer(round(delta))
  if (delta_int > max_delta) {
    stop("P-IG skipped because posterior delta exceeds max_delta", call. = FALSE)
  }

  gamma_const <- euler_gamma_const()
  x <- max(init, 1e-8)
  n_save <- floor((n_iter - burnin) / thin)
  samples <- rep(NA_real_, n_save)
  saved <- 0
  b_pre <- delta * (gamma_const - log(delta * mu_ratio))

  for (iter in seq_len(n_iter)) {
    tau <- stats::rgamma(1, shape = delta * x + 1, rate = 1)
    omega_sum <- rpig_sum(delta_int, x, pig_N)
    b <- b_pre + delta * log(tau)
    x_new <- rH(delta + 1, omega_sum, b)
    if (is.finite(x_new) && x_new > 0) {
      x <- x_new
    }
    if (iter > burnin && ((iter - burnin) %% thin == 0)) {
      saved <- saved + 1
      samples[saved] <- x
    }
  }
  list(samples = samples)
}

prior_from_center <- function(delta, center) {
  if (delta <= 0) {
    return(list(delta = 0, eta = NA_real_, mu = NA_real_))
  }
  list(delta = delta, eta = center, mu = exp(digamma(center)))
}

posterior_params_from_stats <- function(n,
                                        x_a,
                                        x_g,
                                        prior_delta = 0,
                                        prior_center = NA_real_) {
  if (prior_delta > 0) {
    if (!is.finite(prior_center) || prior_center <= 0) {
      stop("prior_center must be positive when prior_delta > 0", call. = FALSE)
    }
    prior <- prior_from_center(prior_delta, prior_center)
    delta_post <- prior_delta + n
    eta_post <- (prior_delta * prior$eta + n * x_a) / delta_post
    mu_post <- prior$mu^(prior_delta / delta_post) * x_g^(n / delta_post)
  } else {
    delta_post <- n
    eta_post <- x_a
    mu_post <- x_g
  }
  list(
    delta_post = delta_post,
    eta_post = eta_post,
    mu_post = mu_post,
    mu_ratio = eta_post / mu_post
  )
}

simulate_gamma_stats <- function(alpha_true, n) {
  x <- stats::rgamma(n, shape = alpha_true, rate = 1)
  list(
    x = x,
    x_a = mean(x),
    x_g = exp(mean(log(x)))
  )
}

skip_row <- function(method, message, pig_N = NA_real_) {
  data.frame(
    method = method,
    status = "skipped",
    ks = NA_real_, w1 = NA_real_, mean = NA_real_, var = NA_real_,
    ci_low = NA_real_, ci_high = NA_real_, cover = NA,
    ess = NA_real_, ess_per_iter = NA_real_, ess_per_sec = NA_real_,
    acf1 = NA_real_, acf10 = NA_real_, runtime_sec = NA_real_,
    accept_rate = NA_real_, avg_evals = NA_real_, pig_N = pig_N,
    n_save = NA_integer_, message = message,
    stringsAsFactors = FALSE
  )
}

run_gamma_shape_task <- function(task, settings) {
  set.seed(settings$seed + task$task_id * 1009)

  if (isTRUE(task$fixed_stats)) {
    x_a <- task$x_a
    x_g <- task$x_g
  } else {
    stats <- simulate_gamma_stats(task$alpha_true, task$n)
    x_a <- stats$x_a
    x_g <- stats$x_g
  }

  pp <- posterior_params_from_stats(
    n = task$n,
    x_a = x_a,
    x_g = x_g,
    prior_delta = task$prior_delta,
    prior_center = task$prior_center
  )

  truth <- xi2_grid_truth(
    delta = pp$delta_post,
    mu_ratio = pp$mu_ratio,
    grid_size = settings$grid_size,
    tail_nats = settings$tail_nats
  )

  rows <- list()

  rows[["truth"]] <- data.frame(
    method = "Numerical grid",
    status = "ok",
    ks = 0,
    w1 = 0,
    mean = truth$mean,
    var = truth$var,
    ci_low = truth$q025,
    ci_high = truth$q975,
    cover = truth$q025 <= task$alpha_true && task$alpha_true <= truth$q975,
    ess = NA_real_, ess_per_iter = NA_real_, ess_per_sec = NA_real_,
    acf1 = NA_real_, acf10 = NA_real_, runtime_sec = NA_real_,
    accept_rate = NA_real_, avg_evals = NA_real_, pig_N = NA_real_,
    n_save = NA_integer_, message = "",
    stringsAsFactors = FALSE
  )

  if (settings$run_approximations) {
    t0 <- proc.time()[["elapsed"]]
    mp <- miller_gamma_params(pp$delta_post, pp$mu_ratio, init = truth$mode)
    rows[["miller"]] <- gamma_approx_row(
      truth = truth,
      shape = mp$shape,
      rate = mp$rate,
      method = "Miller-Gamma",
      runtime_sec = proc.time()[["elapsed"]] - t0,
      alpha_true = task$alpha_true,
      status = if (isTRUE(mp$ok)) "ok" else "failed",
      message = if (isTRUE(mp$ok)) "" else "Miller derivative matching did not converge to valid gamma parameters"
    )

    if (settings$run_miller_imh) {
      if (isTRUE(mp$ok)) {
        t0 <- proc.time()[["elapsed"]]
        mih <- sample_miller_imh(
          delta = pp$delta_post,
          mu_ratio = pp$mu_ratio,
          shape = mp$shape,
          rate = mp$rate,
          n_iter = settings$n_iter,
          burnin = settings$burnin,
          thin = settings$thin,
          init = truth$mode
        )
        rows[["miller_imh"]] <- sample_summary_row(
          samples = mih$samples,
          truth = truth,
          method = "Miller-IMH",
          runtime_sec = proc.time()[["elapsed"]] - t0,
          n_iter = settings$n_iter,
          alpha_true = task$alpha_true,
          accept_rate = mih$accept_rate
        )
      } else {
        rows[["miller_imh"]] <- skip_row(
          "Miller-IMH",
          "Miller-IMH skipped: derivative matching failed"
        )
      }
    }

    t0 <- proc.time()[["elapsed"]]
    sp <- stirling_gamma_params(pp$delta_post, pp$mu_ratio)
    rows[["stirling"]] <- gamma_approx_row(
      truth = truth,
      shape = sp$shape,
      rate = sp$rate,
      method = "Stirling-Gamma",
      runtime_sec = proc.time()[["elapsed"]] - t0,
      alpha_true = task$alpha_true,
      status = if (isTRUE(sp$ok)) "ok" else "failed",
      message = if (isTRUE(sp$ok)) "" else "Stirling approximation produced invalid gamma parameters"
    )
  }

  if (settings$run_slice) {
    t0 <- proc.time()[["elapsed"]]
    sl <- sample_slice_xi2(
      delta = pp$delta_post,
      mu_ratio = pp$mu_ratio,
      n_iter = settings$n_iter,
      burnin = settings$burnin,
      thin = settings$thin,
      init = truth$mode,
      w = settings$slice_w,
      m = settings$slice_m
    )
    rows[["slice"]] <- sample_summary_row(
      samples = sl$samples,
      truth = truth,
      method = "Slice-log-alpha",
      runtime_sec = proc.time()[["elapsed"]] - t0,
      n_iter = settings$n_iter,
      alpha_true = task$alpha_true,
      avg_evals = sl$avg_evals
    )
  }

  if (settings$run_pig) {
    pig_N_values <- if (identical(task$suite, "truncation")) settings$truncation_N else settings$pig_N
    for (pig_N in pig_N_values) {
      method_name <- if (identical(task$suite, "truncation")) {
        paste0("P-IG Gibbs N=", pig_N)
      } else {
        "P-IG Gibbs"
      }
      if (!is_integerish(pp$delta_post)) {
        rows[[paste0("pig_", pig_N)]] <- skip_row(
          method_name,
          "P-IG Gibbs skipped: posterior delta is non-integer",
          pig_N = pig_N
        )
        next
      }
      if (pp$delta_post > settings$pig_max_delta) {
        rows[[paste0("pig_", pig_N)]] <- skip_row(
          method_name,
          paste0("P-IG Gibbs skipped: posterior delta ", pp$delta_post,
                 " exceeds pig_max_delta ", settings$pig_max_delta),
          pig_N = pig_N
        )
        next
      }
      t0 <- proc.time()[["elapsed"]]
      pig <- tryCatch(
        sample_pig_xi2(
          delta = pp$delta_post,
          mu_ratio = pp$mu_ratio,
          n_iter = settings$n_iter,
          burnin = settings$burnin,
          thin = settings$thin,
          pig_N = pig_N,
          init = truth$mode,
          max_delta = settings$pig_max_delta
        ),
        error = function(e) e
      )
      runtime <- proc.time()[["elapsed"]] - t0
      if (inherits(pig, "error")) {
        rows[[paste0("pig_", pig_N)]] <- skip_row(method_name, pig$message, pig_N = pig_N)
      } else {
        rows[[paste0("pig_", pig_N)]] <- sample_summary_row(
          samples = pig$samples,
          truth = truth,
          method = method_name,
          runtime_sec = runtime,
          n_iter = settings$n_iter,
          alpha_true = task$alpha_true,
          pig_N = pig_N
        )
      }
    }
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  cbind(
    task[rep(1, nrow(out)), c(
      "task_id", "suite", "scenario", "rep", "alpha_true", "n",
      "prior_delta", "prior_center", "fixed_stats"
    )],
    x_a = x_a,
    x_g = x_g,
    delta_post = pp$delta_post,
    mu_ratio = pp$mu_ratio,
    truth_mode = truth$mode,
    out,
    stringsAsFactors = FALSE
  )
}

make_damsleth_tasks <- function(reps = 100) {
  base <- data.frame(
    suite = "damsleth",
    scenario = c("Damsleth n=5", "Damsleth n=10", "Damsleth n=30"),
    alpha_true = 5,
    n = c(5, 10, 30),
    x_a = c(7.19, 5.57, 5.09),
    x_g = c(6.05, 5.01, 4.26),
    prior_delta = 0,
    prior_center = NA_real_,
    fixed_stats = TRUE,
    stringsAsFactors = FALSE
  )
  tasks <- base[rep(seq_len(nrow(base)), each = reps), ]
  tasks$rep <- rep(seq_len(reps), times = nrow(base))
  tasks
}

make_main_tasks <- function(reps = 200,
                            alpha_values = c(0.1, 0.5, 1, 2, 5, 20),
                            n_values = c(10, 50, 500)) {
  grid <- expand.grid(
    alpha_true = alpha_values,
    n = n_values,
    rep = seq_len(reps),
    KEEP.OUT.ATTRS = FALSE
  )
  data.frame(
    suite = "main",
    scenario = paste0("alpha=", grid$alpha_true, ", n=", grid$n),
    alpha_true = grid$alpha_true,
    n = grid$n,
    x_a = NA_real_,
    x_g = NA_real_,
    prior_delta = 0,
    prior_center = NA_real_,
    fixed_stats = FALSE,
    rep = grid$rep,
    stringsAsFactors = FALSE
  )
}

make_prior_tasks <- function(reps = 100,
                             alpha_true = 2,
                             n = 30,
                             prior_deltas = c(0.5, 2.5, 10, 50),
                             prior_centers = c(2, 0.5, 10)) {
  grid <- expand.grid(
    prior_delta = prior_deltas,
    prior_center = prior_centers,
    rep = seq_len(reps),
    KEEP.OUT.ATTRS = FALSE
  )
  data.frame(
    suite = "prior",
    scenario = paste0("delta=", grid$prior_delta, ", center=", grid$prior_center),
    alpha_true = alpha_true,
    n = n,
    x_a = NA_real_,
    x_g = NA_real_,
    prior_delta = grid$prior_delta,
    prior_center = grid$prior_center,
    fixed_stats = FALSE,
    rep = grid$rep,
    stringsAsFactors = FALSE
  )
}

make_truncation_tasks <- function(reps = 50,
                                  alpha_values = c(0.1, 2, 20),
                                  n_values = c(10, 50)) {
  grid <- expand.grid(
    alpha_true = alpha_values,
    n = n_values,
    rep = seq_len(reps),
    KEEP.OUT.ATTRS = FALSE
  )
  data.frame(
    suite = "truncation",
    scenario = paste0("alpha=", grid$alpha_true, ", n=", grid$n),
    alpha_true = grid$alpha_true,
    n = grid$n,
    x_a = NA_real_,
    x_g = NA_real_,
    prior_delta = 0,
    prior_center = NA_real_,
    fixed_stats = FALSE,
    rep = grid$rep,
    stringsAsFactors = FALSE
  )
}

make_extreme_tasks <- function(reps = 20,
                               alpha_values = c(0.05, 0.1),
                               n = 100) {
  grid <- expand.grid(
    alpha_true = alpha_values,
    n = n,
    rep = seq_len(reps),
    KEEP.OUT.ATTRS = FALSE
  )
  data.frame(
    suite = "extreme",
    scenario = paste0("alpha=", grid$alpha_true, ", n=", grid$n),
    alpha_true = grid$alpha_true,
    n = grid$n,
    x_a = NA_real_,
    x_g = NA_real_,
    prior_delta = 0,
    prior_center = NA_real_,
    fixed_stats = FALSE,
    rep = grid$rep,
    stringsAsFactors = FALSE
  )
}

make_smoke_tasks <- function() {
  tasks <- rbind(
    make_damsleth_tasks(reps = 1)[1, ],
    make_main_tasks(reps = 1, alpha_values = c(0.5), n_values = c(10)),
    make_prior_tasks(reps = 1, prior_deltas = c(0.5), prior_centers = c(10)),
    make_truncation_tasks(reps = 1, alpha_values = c(0.1), n_values = c(10)),
    make_extreme_tasks(reps = 1, alpha_values = c(0.1), n = 100)
  )
  tasks$suite[seq_len(nrow(tasks))] <- c("damsleth", "main", "prior", "truncation", "extreme")
  tasks
}

make_tasks <- function(suite,
                       profile = "smoke",
                       reps = NULL) {
  if (identical(suite, "smoke")) {
    tasks <- make_smoke_tasks()
  } else if (identical(suite, "damsleth")) {
    tasks <- make_damsleth_tasks(reps %||% if (profile == "full") 100 else 2)
  } else if (identical(suite, "main")) {
    if (profile == "full") {
      tasks <- make_main_tasks(reps %||% 200)
    } else {
      tasks <- make_main_tasks(reps %||% 2, alpha_values = c(0.1, 2), n_values = c(10, 50))
    }
  } else if (identical(suite, "prior")) {
    if (profile == "full") {
      tasks <- make_prior_tasks(reps %||% 100)
    } else {
      tasks <- make_prior_tasks(reps %||% 2, prior_deltas = c(0.5, 10), prior_centers = c(2, 10))
    }
  } else if (identical(suite, "truncation")) {
    if (profile == "full") {
      tasks <- make_truncation_tasks(reps %||% 50)
    } else {
      tasks <- make_truncation_tasks(reps %||% 2, alpha_values = c(0.1, 2), n_values = c(10))
    }
  } else if (identical(suite, "extreme")) {
    tasks <- make_extreme_tasks(reps %||% if (profile == "full") 20 else 2)
  } else if (identical(suite, "all")) {
    tasks <- rbind(
      make_tasks("damsleth", profile, reps),
      make_tasks("main", profile, reps),
      make_tasks("prior", profile, reps),
      make_tasks("truncation", profile, reps),
      make_tasks("extreme", profile, reps)
    )
  } else {
    stop("Unknown suite: ", suite, call. = FALSE)
  }
  tasks$task_id <- seq_len(nrow(tasks))
  tasks[, c("task_id", setdiff(names(tasks), "task_id"))]
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

summarise_gamma_shape_results <- function(raw) {
  ok <- raw[raw$status == "ok", , drop = FALSE]
  if (nrow(ok) == 0) {
    return(data.frame())
  }
  group_cols <- c("suite", "scenario", "alpha_true", "n", "prior_delta",
                  "prior_center", "method", "pig_N")
  key_df <- ok[group_cols]
  for (nm in names(key_df)) {
    key_df[[nm]] <- ifelse(is.na(key_df[[nm]]), "<NA>", as.character(key_df[[nm]]))
  }
  group_key <- do.call(paste, c(key_df, sep = "\r"))
  groups <- split(ok, group_key)
  mc_se <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) <= 1) return(0)
    stats::sd(x) / sqrt(length(x))
  }
  first_nonmissing <- function(x) {
    x <- x[!is.na(x)]
    if (length(x)) x[1] else NA
  }
  do.call(rbind, lapply(groups, function(d) {
    data.frame(
      suite = d$suite[1],
      scenario = d$scenario[1],
      alpha_true = d$alpha_true[1],
      n = d$n[1],
      prior_delta = d$prior_delta[1],
      prior_center = d$prior_center[1],
      method = d$method[1],
      pig_N = first_nonmissing(d$pig_N),
      reps = nrow(d),
      ks_mean = mean(d$ks, na.rm = TRUE),
      ks_se = mc_se(d$ks),
      w1_mean = mean(d$w1, na.rm = TRUE),
      w1_se = mc_se(d$w1),
      ess_iter_median = stats::median(d$ess_per_iter, na.rm = TRUE),
      ess_sec_median = stats::median(d$ess_per_sec, na.rm = TRUE),
      acf1_mean = mean(d$acf1, na.rm = TRUE),
      acf10_mean = mean(d$acf10, na.rm = TRUE),
      coverage = mean(d$cover, na.rm = TRUE),
      runtime_median = stats::median(d$runtime_sec, na.rm = TRUE),
      runtime_p90 = as.numeric(stats::quantile(d$runtime_sec, 0.9, na.rm = TRUE, names = FALSE)),
      mean_error = mean(d$mean - d$alpha_true, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
}
