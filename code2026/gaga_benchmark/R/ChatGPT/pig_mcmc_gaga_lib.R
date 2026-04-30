# PIG-augmented per-gene posterior sampler for the GaGa model.
#
# This file intentionally does not compute per-pattern marginal likelihoods by
# numerical integration.  For fixed GaGa hyperparameters it samples
# (z_i, alpha_i) directly:
#   z_i | alpha_i, x_i  by exact collapsed pattern probabilities,
#   alpha_i | z_i, x_i  by PIG augmentation plus an exact log-alpha
#                       slice update for the augmented conditional.

pig_rPIG_original <- function(n, cc, N) {
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

pig_pattern_stats <- function(x, groups, patterns) {
  groups <- factor(groups)
  if (!is.null(colnames(patterns))) {
    group_idx <- match(as.character(groups), colnames(patterns))
  } else {
    group_idx <- as.integer(groups)
  }
  if (anyNA(group_idx)) {
    stop("Some group labels are not present in pattern column names.")
  }

  x <- pmax(as.matrix(x), .Machine$double.eps)
  log_x <- log(x)
  out <- vector("list", nrow(patterns))
  for (h in seq_len(nrow(patterns))) {
    cluster <- as.integer(patterns[h, group_idx])
    ids <- sort(unique(cluster))
    n_c <- integer(length(ids))
    sums <- matrix(0, nrow(x), length(ids))
    log_sums <- matrix(0, nrow(x), length(ids))
    for (j in seq_along(ids)) {
      cols <- which(cluster == ids[j])
      n_c[j] <- length(cols)
      sums[, j] <- rowSums(x[, cols, drop = FALSE])
      log_sums[, j] <- rowSums(log_x[, cols, drop = FALSE])
    }
    out[[h]] <- list(clusters = ids, n_c = n_c, sums = sums,
                     log_sums = log_sums)
  }
  out
}

pig_log_likelihood_alpha <- function(alpha, n_c, S_c, L_c, alpha0, nu) {
  if (!is.finite(alpha) || alpha <= 0) {
    return(-Inf)
  }
  rate0 <- alpha0 / nu
  A <- alpha0 + alpha * n_c
  B <- rate0 + alpha * S_c
  if (any(!is.finite(B)) || any(B <= 0)) {
    return(-Inf)
  }
  sum(
    n_c * alpha * log(alpha) -
      n_c * lgamma(alpha) +
      (alpha - 1) * L_c +
      alpha0 * log(rate0) -
      lgamma(alpha0) +
      lgamma(A) -
      A * log(B)
  )
}

pig_log_alpha_target <- function(alpha, n_c, S_c, L_c, hyper) {
  pig_log_likelihood_alpha(alpha, n_c, S_c, L_c, hyper$alpha0, hyper$nu) +
    (hyper$b_alpha - 1) * log(alpha) -
    (hyper$b_alpha / hyper$nu_alpha) * alpha
}

pig_alpha_remainder <- function(alpha, n_c, S_c, hyper) {
  if (!is.finite(alpha) || alpha <= 0) {
    return(-Inf)
  }
  rate0 <- hyper$alpha0 / hyper$nu
  B <- rate0 + alpha * S_c
  if (any(!is.finite(B)) || any(B <= 0)) {
    return(-Inf)
  }
  sum(n_c) * alpha * log(alpha) -
    sum((hyper$alpha0 + alpha * n_c) * log(B))
}

pig_alpha_remainder_slope <- function(alpha, n_c, S_c, hyper) {
  if (!is.finite(alpha) || alpha <= 0) {
    return(NA_real_)
  }
  rate0 <- hyper$alpha0 / hyper$nu
  B <- rate0 + alpha * S_c
  if (any(!is.finite(B)) || any(B <= 0)) {
    return(NA_real_)
  }
  sum(n_c) * (log(alpha) + 1) -
    sum(n_c * log(B) + (hyper$alpha0 + alpha * n_c) * S_c / B)
}

pig_augmented_log_alpha_y <- function(y, power, omega_sum, lin, n_c, S_c,
                                      hyper) {
  if (!is.finite(y)) {
    return(-Inf)
  }
  alpha <- exp(y)
  if (!is.finite(alpha) || alpha <= 0) {
    return(-Inf)
  }
  power * y - omega_sum * alpha^2 + lin * alpha +
    pig_alpha_remainder(alpha, n_c, S_c, hyper)
}

slice_sample_log_alpha <- function(y, log_density, width = 0.5,
                                   max_steps = 50, max_shrink = 200) {
  cur_log <- log_density(y)
  if (!is.finite(cur_log)) {
    return(list(y = y, moved = FALSE))
  }
  level <- cur_log - stats::rexp(1)
  u <- stats::runif(1)
  left <- y - width * u
  right <- left + width

  n_left <- sample.int(max_steps + 1L, 1L) - 1L
  n_right <- max_steps - n_left
  while (n_left > 0L && log_density(left) > level) {
    left <- left - width
    n_left <- n_left - 1L
  }
  while (n_right > 0L && log_density(right) > level) {
    right <- right + width
    n_right <- n_right - 1L
  }

  for (i in seq_len(max_shrink)) {
    y_new <- stats::runif(1, left, right)
    if (log_density(y_new) >= level) {
      return(list(y = y_new, moved = abs(y_new - y) > 1e-10))
    }
    if (y_new < y) {
      left <- y_new
    } else {
      right <- y_new
    }
  }
  list(y = y, moved = FALSE)
}

pig_alpha_step_slice <- function(alpha, n_c, S_c, L_c, hyper, trunc = 80,
                                 slice_width = 0.5, slice_max_steps = 50) {
  if (!is.finite(alpha) || alpha <= 0) {
    alpha <- max(0.1, min(hyper$nu_alpha, 100))
  }
  n_total <- sum(n_c)
  omega_sum <- sum(pig_rPIG_original(n_total, alpha, trunc))
  tau <- stats::rgamma(length(n_c), shape = hyper$alpha0 + alpha * n_c,
                       rate = 1)
  log_tau <- log(pmax(tau, .Machine$double.xmin))

  power <- n_total + hyper$b_alpha
  lin <- n_total * gamma_const +
    sum(L_c) -
    (hyper$b_alpha / hyper$nu_alpha) +
    sum(n_c * log_tau)

  y <- log(alpha)
  log_density <- function(value) {
    pig_augmented_log_alpha_y(value, power, omega_sum, lin, n_c, S_c, hyper)
  }
  step <- slice_sample_log_alpha(
    y,
    log_density,
    width = slice_width,
    max_steps = slice_max_steps
  )
  alpha_new <- exp(step$y)
  if (!is.finite(alpha_new) || alpha_new <= 0) {
    return(list(alpha = alpha, accepted = FALSE, proposed = TRUE))
  }
  list(alpha = alpha_new, accepted = step$moved, proposed = TRUE)
}

pig_alpha_step_ptn <- function(alpha, n_c, S_c, L_c, hyper, trunc = 80,
                               max_try = 5000) {
  if (!is.finite(alpha) || alpha <= 0) {
    alpha <- max(0.1, min(hyper$nu_alpha, 100))
  }
  n_total <- sum(n_c)
  omega_sum <- sum(pig_rPIG_original(n_total, alpha, trunc))
  tau <- stats::rgamma(length(n_c), shape = hyper$alpha0 + alpha * n_c,
                       rate = 1)
  log_tau <- log(pmax(tau, .Machine$double.xmin))

  power <- n_total + hyper$b_alpha
  lin <- n_total * gamma_const +
    sum(L_c) -
    (hyper$b_alpha / hyper$nu_alpha) +
    sum(n_c * log_tau)
  proposal <- rPTN(power = power, quad = omega_sum, lin = lin,
                   max_try = max_try)
  if (!is.finite(proposal) || proposal <= 0) {
    return(list(alpha = alpha, accepted = FALSE, proposed = FALSE))
  }

  log_acc <- pig_alpha_remainder(proposal, n_c, S_c, hyper) -
    pig_alpha_remainder(alpha, n_c, S_c, hyper)
  if (is.finite(log_acc) && log(stats::runif(1)) < log_acc) {
    return(list(alpha = proposal, accepted = TRUE, proposed = TRUE))
  }
  list(alpha = alpha, accepted = FALSE, proposed = TRUE)
}

pig_alpha_step <- function(alpha, n_c, S_c, L_c, hyper, trunc = 80,
                           kernel = c("slice", "ptn"),
                           slice_width = 0.5, slice_max_steps = 50,
                           max_try = 5000) {
  kernel <- match.arg(kernel)
  if (kernel == "slice") {
    return(pig_alpha_step_slice(
      alpha, n_c, S_c, L_c, hyper,
      trunc = trunc,
      slice_width = slice_width,
      slice_max_steps = slice_max_steps
    ))
  }
  pig_alpha_step_ptn(alpha, n_c, S_c, L_c, hyper, trunc = trunc,
                     max_try = max_try)
}

pig_init_alpha <- function(y, hyper) {
  y <- pmax(as.numeric(y), .Machine$double.eps)
  m <- mean(y)
  v <- stats::var(y)
  mom <- if (is.finite(v) && v > 0) m^2 / v else hyper$nu_alpha
  alpha <- if (is.finite(mom) && mom > 0) mom else hyper$nu_alpha
  max(0.05, min(alpha, 5000))
}

pig_pattern_logprob <- function(alpha, gene_index, stats, probpat, hyper) {
  out <- numeric(length(stats))
  for (h in seq_along(stats)) {
    st <- stats[[h]]
    out[h] <- log(probpat[h]) +
      pig_log_likelihood_alpha(
        alpha,
        st$n_c,
        st$sums[gene_index, ],
        st$log_sums[gene_index, ],
        hyper$alpha0,
        hyper$nu
      )
  }
  out
}

sample_log_prob <- function(logp) {
  probs <- pattern_prob_from_log(logp)
  base::sample.int(length(logp), 1, prob = probs)
}

pattern_prob_from_log <- function(logp) {
  if (all(!is.finite(logp))) {
    return(rep(1 / length(logp), length(logp)))
  }
  softmax_log(logp)
}

pig_mcmc_gaga_pp <- function(x, groups, patterns, fit, max_genes = 0,
                             n_iter = 600, burnin = 200, thin = 1,
                             trunc = 80, alpha_steps = 1,
                             alpha_kernel = c("slice", "ptn"),
                             slice_width = 0.5, slice_max_steps = 50,
                             progress_every = 200,
                             seed = NULL) {
  require_pkgs(c("gaga", "GIGrvg"))
  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (max_genes > 0 && max_genes < nrow(x)) {
    x <- x[seq_len(max_genes), , drop = FALSE]
  }
  if (n_iter <= burnin) {
    stop("n_iter must be larger than burnin.")
  }
  if (thin <= 0) {
    stop("thin must be positive.")
  }
  alpha_kernel <- match.arg(alpha_kernel)

  par <- gaga::getpar(fit)
  if (length(par$a0) != 1 || length(par$nu) != 1) {
    stop("pig_mcmc_gaga_pp currently supports GaGa nclust=1 only.")
  }
  hyper <- list(
    alpha0 = as.numeric(par$a0[1]),
    nu = as.numeric(par$nu[1]),
    b_alpha = as.numeric(par$balpha[1]),
    nu_alpha = as.numeric(par$nualpha[1])
  )
  probpat <- as.numeric(par$probpat)
  if (length(probpat) != nrow(patterns) || any(!is.finite(probpat))) {
    probpat <- rep(1 / nrow(patterns), nrow(patterns))
  }
  probpat <- pmax(probpat, .Machine$double.eps)
  probpat <- probpat / sum(probpat)

  x <- pmax(as.matrix(x), .Machine$double.eps)
  stats <- pig_pattern_stats(x, groups, patterns)
  pp <- matrix(NA_real_, nrow = nrow(x), ncol = nrow(patterns))
  colnames(pp) <- paste0("pattern", seq_len(nrow(patterns)) - 1)
  rownames(pp) <- rownames(x)
  alpha_mean <- numeric(nrow(x))
  accept_rate <- numeric(nrow(x))
  kept_draws <- integer(nrow(x))
  z_counts <- matrix(0L, nrow = nrow(x), ncol = nrow(patterns))
  colnames(z_counts) <- colnames(pp)
  rownames(z_counts) <- rownames(x)

  for (i in seq_len(nrow(x))) {
    if (progress_every > 0 && (i == 1 || i %% progress_every == 0)) {
      cat("  PIG-MCMC scoring gene ", i, "/", nrow(x), "\n", sep = "")
    }
    alpha <- pig_init_alpha(x[i, ], hyper)
    z <- sample_log_prob(pig_pattern_logprob(alpha, i, stats, probpat, hyper))
    counts <- integer(nrow(patterns))
    prob_sum <- numeric(nrow(patterns))
    alpha_sum <- 0
    n_keep <- 0L
    n_accept <- 0L
    n_prop <- 0L

    for (iter in seq_len(n_iter)) {
      z_prob <- pattern_prob_from_log(
        pig_pattern_logprob(alpha, i, stats, probpat, hyper)
      )
      z <- base::sample.int(nrow(patterns), 1, prob = z_prob)
      st <- stats[[z]]
      for (s in seq_len(alpha_steps)) {
        step <- pig_alpha_step(
          alpha,
          st$n_c,
          st$sums[i, ],
          st$log_sums[i, ],
          hyper,
          trunc = trunc,
          kernel = alpha_kernel,
          slice_width = slice_width,
          slice_max_steps = slice_max_steps
        )
        alpha <- step$alpha
        n_accept <- n_accept + as.integer(step$accepted)
        n_prop <- n_prop + as.integer(step$proposed)
      }
      if (iter > burnin && ((iter - burnin) %% thin == 0)) {
        z_prob_keep <- pattern_prob_from_log(
          pig_pattern_logprob(alpha, i, stats, probpat, hyper)
        )
        counts[z] <- counts[z] + 1L
        prob_sum <- prob_sum + z_prob_keep
        alpha_sum <- alpha_sum + alpha
        n_keep <- n_keep + 1L
      }
    }
    if (n_keep == 0L) {
      stop("No post-burnin draws retained for gene ", i)
    }
    z_counts[i, ] <- counts
    pp[i, ] <- prob_sum / n_keep
    alpha_mean[i] <- alpha_sum / n_keep
    accept_rate[i] <- if (n_prop > 0) n_accept / n_prop else NA_real_
    kept_draws[i] <- n_keep
  }

  list(
    pp = pp,
    z_counts = z_counts,
    alpha_mean = alpha_mean,
    accept_rate = accept_rate,
    par = par,
    hyper = hyper,
    config = list(n_iter = n_iter, burnin = burnin, thin = thin,
                  trunc = trunc, alpha_steps = alpha_steps,
                  alpha_kernel = alpha_kernel,
                  slice_width = slice_width,
                  slice_max_steps = slice_max_steps)
  )
}

find_pig_mcmc_calls <- function(pp, fdr = 0.05) {
  null <- pp[, 1]
  ord <- order(null)
  fdr_path <- cumsum(null[ord]) / seq_along(ord)
  keep_ord <- which(fdr_path <= fdr)
  calls <- rep(FALSE, nrow(pp))
  if (length(keep_ord) > 0) {
    calls[ord[seq_len(max(keep_ord))]] <- TRUE
  }
  calls
}
