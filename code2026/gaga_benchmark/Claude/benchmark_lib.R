# =============================================================================
# benchmark_lib.R  (Claude/ debugged version)
#
# Differences from gaga_benchmark/R/benchmark_lib.R:
#   * limma_or_ttest -> renamed to limma_BH; takes data on the "modeling" scale
#     (NOT raw). Does NOT apply log2 internally. The caller decides the scale.
#     Old code applied log2 after the caller already log2'd, double-logging.
#   * Added limma_BH_raw: convenience wrapper that takes raw positive expression
#     and applies log2(pmax(., eps)) before limma.
#   * Added ebarrays_GG_raw: ensure EBarrays GG is fit on positive expression
#     and not on log-scale (the GG model is gamma-on-raw, not gamma-on-log).
#   * Added cdf_safe_drop_u95av2: helper to drop U95Av2 MLL arrays based on CDF
#     name when available, else falls back to "tail-2" with an explicit warning.
# =============================================================================

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
  # script lives in gaga_benchmark/Claude, so benchmark dir is one level up.
  normalizePath(file.path(dirname(script_path()), ".."), mustWork = TRUE)
}

setup_benchmark <- function() {
  bench <- benchmark_dir()
  dir.create(file.path(bench, "Claude", "results"), recursive = TRUE,
             showWarnings = FALSE)
  # Only prepend the project's r-lib if it appears to contain libraries built
  # for the current platform. On macOS the legacy r-lib has Darwin .so files
  # that fail to load on Linux; prepending them shadows the system library
  # and produces "package found but cannot load" errors. The sandbox should
  # rely on the system R library exclusively.
  add_lib <- function(p) {
    if (!dir.exists(p)) return(FALSE)
    pkgs <- list.files(p, full.names = TRUE)
    if (length(pkgs) == 0) return(FALSE)
    # Some r-lib trees contain a mix of Mach-O (macOS) and ELF (Linux) .so
    # files. To be safe we scan every .so/.dylib in the tree and reject the
    # whole library if ANY entry is Mach-O while R is running on Linux.
    so_files <- list.files(file.path(pkgs, "libs"),
                           pattern = "\\.(so|dylib)$",
                           full.names = TRUE, recursive = TRUE)
    if (length(so_files) == 0) return(TRUE)
    is_linux_run <- grepl("linux", R.version$os, ignore.case = TRUE)
    if (!is_linux_run) return(TRUE)
    has_macho <- function(f) {
      tryCatch({
        con <- file(f, "rb"); on.exit(close(con))
        magic <- readBin(con, "raw", n = 4)
        m <- paste(format(magic, width = 2), collapse = "")
        grepl("cefaedfe|cffaedfe|feedface|feedfacf", m)
      }, error = function(e) FALSE)
    }
    if (any(vapply(so_files, has_macho, logical(1)))) {
      message("[setup_benchmark] Skipping ", p,
              " (contains Mach-O libs incompatible with Linux R).")
      return(FALSE)
    }
    TRUE
  }
  user_lib <- file.path(bench, "r-lib")
  legacy_lib <- normalizePath(file.path(bench, "..", "r-lib"),
                              mustWork = FALSE)
  add_paths <- character()
  for (p in c(user_lib, legacy_lib)) {
    if (add_lib(p)) add_paths <- c(add_paths, p)
  }
  if (length(add_paths) > 0) {
    .libPaths(unique(c(add_paths, .libPaths())))
  }
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

# All Claude/ outputs go to Claude/results so we never overwrite the legacy
# R/ output.
result_file <- function(name) {
  file.path(benchmark_dir(), "Claude", "results", name)
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

# ---------------------------------------------------------------------------
# limma helpers - the data must already be on the modeling scale (log).
# These do NOT log-transform. Use limma_BH_raw if your input is raw expression.
# ---------------------------------------------------------------------------
limma_BH <- function(x, groups, fdr = 0.05) {
  groups <- factor(groups)
  if (!requireNamespace("limma", quietly = TRUE)) {
    idx <- split(seq_along(groups), groups)
    if (length(idx) != 2) {
      stop("Base R fallback: only two-group t-test supported.")
    }
    p <- apply(x, 1, function(row) {
      stats::t.test(row[idx[[1]]], row[idx[[2]]])$p.value
    })
  } else {
    design <- stats::model.matrix(~ groups)
    fit <- limma::lmFit(x, design)
    fit <- limma::eBayes(fit)
    if (ncol(design) == 2) {
      p <- fit$p.value[, 2]
    } else {
      tt <- limma::topTable(fit, coef = seq_len(ncol(design))[-1],
                            number = Inf, sort.by = "none")
      p <- tt$P.Value
      names(p) <- rownames(tt)
    }
  }
  adj <- stats::p.adjust(p, method = "BH")
  list(p = p, adj_p = adj, calls = adj <= fdr,
       score = -log10(pmax(p, 1e-300)))
}

limma_BH_raw <- function(x_raw, groups, fdr = 0.05,
                         eps = .Machine$double.eps) {
  limma_BH(log2(pmax(x_raw, eps)), groups, fdr)
}

# Backwards-compatible alias used by the legacy R/ scripts. Kept here only so
# scripts that source this file via `source("../R/benchmark_lib.R")` still
# work, but the alias takes log2 ONCE - matching what the old single-log
# expectation should have been for raw input.
limma_or_ttest <- function(x, groups, fdr = 0.05) {
  warning("limma_or_ttest is deprecated; use limma_BH (log-scale input) or ",
          "limma_BH_raw (raw input). This alias treats x as RAW.",
          call. = FALSE)
  limma_BH_raw(x, groups, fdr)
}

# ---------------------------------------------------------------------------
# EBarrays Ga (Gamma-Gamma) wrapper.
# The GG model assumes positive (raw / linear-scale) expression. Passing
# log-scale data here is a model misspecification and was a bug in the
# legacy 10_armstrong_simulation62.R.
# ---------------------------------------------------------------------------
make_eb_two_group_pattern <- function(groups) {
  groups <- factor(groups)
  lv <- levels(groups)
  EBarrays::ebPatterns(c(
    paste(rep(1, length(groups)), collapse = " "),
    paste(ifelse(groups == lv[1], 1, 2), collapse = " ")
  ))
}

ebarrays_GG_raw <- function(x_raw, groups, fdr = 0.05) {
  if (any(!is.finite(x_raw)) || min(x_raw, na.rm = TRUE) <= 0) {
    stop("EBarrays GG requires strictly positive expression values.")
  }
  fit <- EBarrays::emfit(
    data = x_raw,
    family = "GG",
    hypotheses = make_eb_two_group_pattern(groups),
    num.iter = 20,
    verbose = FALSE
  )
  post <- EBarrays::postprob(fit, x_raw)$pattern
  threshold <- EBarrays::crit.fun(post[, 1], fdr)
  list(fit = fit, post = post, threshold = threshold,
       calls = post[, 2] > threshold,
       score = post[, 2])
}

# ---------------------------------------------------------------------------
# CDF-based U95Av2 dropping (replaces blind tail-2)
# ---------------------------------------------------------------------------
cdf_safe_drop_u95av2 <- function(groups, cdf = NULL, n_drop_fallback = 2) {
  groups <- as.character(groups)
  if (!is.null(cdf)) {
    cdf <- as.character(cdf)
    is_u95av2 <- grepl("u95av2", cdf, ignore.case = TRUE) |
                 grepl("hg.*u95av2", cdf, ignore.case = TRUE)
    drop_idx <- which(groups == "MLL" & is_u95av2)
    if (length(drop_idx) > 0) {
      message(sprintf(
        "[cdf_safe_drop_u95av2] Dropping %d MLL arrays whose CDF matches U95Av2.",
        length(drop_idx)
      ))
      return(drop_idx)
    }
    warning("CDF info supplied but no MLL arrays flagged as U95Av2; ",
            "falling back to tail-", n_drop_fallback, ".",
            call. = FALSE)
  }
  mll_idx <- which(groups == "MLL")
  if (length(mll_idx) <= n_drop_fallback) {
    return(integer(0))
  }
  warning("CDF info not available; dropping the LAST ", n_drop_fallback,
          " MLL arrays as a heuristic. Verify your input file ordering.",
          call. = FALSE)
  utils::tail(mll_idx, n_drop_fallback)
}

# ---------------------------------------------------------------------------
# GaGa shape-step kernel and exact moments  (unchanged math; correct)
# ---------------------------------------------------------------------------
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
    upper <- max(20, stats::qgamma(0.999999, shape = ar["shape"],
                                   rate = ar["rate"]) * 4)
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
    if (is.finite(log_mode) && (!is.finite(log_upper) ||
                                log_upper < log_mode - 35)) {
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
      z_lower, z_upper, subdivisions = 1000, rel.tol = rel.tol
    )$value
    log(val) + log_scale
  }
  log_z0 <- log_moment(0)
  log_z1 <- log_moment(1)
  log_z2 <- log_moment(2)
  mean <- exp(log_z1 - log_z0)
  second <- exp(log_z2 - log_z0)
  var <- second - mean^2
  if (var < 0 && abs(var) < 1e-8 * max(1, second)) var <- 0
  mode <- exp(stats::optimize(
    function(z) -gas_log_kernel(exp(z), a, b, c, d, r, s),
    interval = c(z_lower, z_upper)
  )$minimum)
  list(mean = mean, var = var, norm = exp(log_z0), mode = mode, upper = upper)
}

# ---------------------------------------------------------------------------
# PIG sampler  (unchanged math; verified)
# ---------------------------------------------------------------------------
rPIG_ERGamma <- function(n, cc, N) {
  require_pkgs("GIGrvg")
  output <- rep(0, n)
  if (N > 0) {
    for (i in seq_len(N)) {
      if (cc == 0) {
        output <- output + 1 / stats::rgamma(n, shape = 3 / 2) / (4 * i^2)
      } else {
        output <- output + GIGrvg::rgig(
          n = n, lambda = -3 / 2, chi = 1 / (2 * i^2), psi = 2 * cc^2
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
    if (!is.finite(xi)) stop("rH failed to accept a draw.")
    if (iter > burnin) out[iter - burnin] <- xi
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
