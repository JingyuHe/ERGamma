# =============================================================================
# pig_gaga_slice_lib.R
#
# Per-gene PIG-augmented MCMC for the GaGa model where every iteration
#   * draws PIG/Damsleth auxiliaries (ω_s ~ PIG(α), τ_c ~ Ga(a_0 + n_c α, 1)),
#   * updates α by EXACT slice sampling on log α of the augmented conditional,
#   * accumulates Rao-Blackwellised pattern probabilities.
#
# This replaces Rossell's Stirling-based gamma-shape approximation
# (gaga::fitGG, gapprox = TRUE in ppGG) with exact PIG-augmented MCMC on the
# full integrated conditional. Hyperparameters (a0, nu, b_alpha, nu_alpha,
# pi) are taken from gaga::fitGG so the only change is Stirling → PIG.
#
# Algorithm follows PIG_GAGA_MCMC_FORMULAS.md and PIG_GAGA_MCMC_ALGORITHM.md.
# =============================================================================

if (!requireNamespace("GIGrvg", quietly = TRUE))
  stop("Package GIGrvg required (rgig).")
library(GIGrvg)

gamma_const <- -digamma(1)

# ---------- PIG draw (truncated GIG sum + gamma tail correction) ---------
rPIG_ERGamma <- function(n, cc, N) {
  output <- rep(0, n)
  if (N > 0) {
    for (i in seq_len(N)) {
      if (cc == 0) {
        output <- output + 1 / stats::rgamma(n, shape = 3 / 2) / (4 * i^2)
      } else {
        output <- output + GIGrvg::rgig(
          n = n, lambda = -3 / 2,
          chi = 1 / (2 * i^2), psi = 2 * cc^2
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

# ---------- Sufficient statistics keyed by pattern --------------------------
gene_cluster_stats_pergene <- function(X, group, patterns) {
  G <- nrow(X); P <- nrow(patterns)
  group_int <- as.integer(group)
  log_obs <- log(pmax(X, .Machine$double.eps))
  res <- vector("list", P)
  for (z_idx in seq_len(P)) {
    z <- as.integer(patterns[z_idx, ])
    cluster <- z[group_int]
    uniq <- sort(unique(cluster))
    n_c <- integer(length(uniq))
    S <- matrix(0, G, length(uniq))
    L <- matrix(0, G, length(uniq))
    for (i in seq_along(uniq)) {
      cols <- which(cluster == uniq[i])
      n_c[i] <- length(cols)
      S[, i] <- rowSums(X[, cols, drop = FALSE])
      L[, i] <- rowSums(log_obs[, cols, drop = FALSE])
    }
    res[[z_idx]] <- list(clusters = uniq, n_c = n_c, S = S, sumlogL = L)
  }
  res
}

# ---------- Per-gene log m_ih(α) (closed form given α) --------------------
log_m_ih_one_gene <- function(alpha, n_c_vec, S_vec, L_vec, alpha0, nu) {
  if (alpha <= 0) return(-Inf)
  rate0 <- alpha0 / nu
  K <- length(n_c_vec); out <- 0
  for (c_idx in seq_len(K)) {
    n <- n_c_vec[c_idx]; S <- S_vec[c_idx]; L <- L_vec[c_idx]
    A <- alpha0 + n * alpha
    B <- rate0 + alpha * S
    out <- out +
      n * alpha * log(alpha) +
      (alpha - 1) * L -
      n * lgamma(alpha) +
      alpha0 * log(rate0) - lgamma(alpha0) +
      lgamma(A) - A * log(B)
  }
  out
}

# Residual R(α) = N α log α − Σ_c (a_0 + n_c α) log(r_0 + α S_c)
R_alpha_one_gene <- function(alpha, n_c_vec, S_vec, alpha0, nu) {
  rate0 <- alpha0 / nu
  N <- sum(n_c_vec)
  out <- N * alpha * log(alpha)
  for (c_idx in seq_along(n_c_vec)) {
    A <- alpha0 + n_c_vec[c_idx] * alpha
    B <- rate0 + alpha * S_vec[c_idx]
    out <- out - A * log(B)
  }
  out
}

# ---------- Augmented log-conditional in y = log α ------------------------
# ell(y) = (N + b_alpha) y − Omega · exp(2y) + B · exp(y) + R(exp(y))
ell_y <- function(y, Omega, B, n_c_vec, S_vec, alpha0, nu, N_i, b_alpha) {
  if (!is.finite(y)) return(-Inf)
  alpha <- exp(y)
  if (!is.finite(alpha) || alpha <= 0 || alpha > 1e12) return(-Inf)
  R <- R_alpha_one_gene(alpha, n_c_vec, S_vec, alpha0, nu)
  (N_i + b_alpha) * y - Omega * exp(2 * y) + B * alpha + R
}

# ---------- Stepping-out + shrinkage slice sampler on a 1D log-density ----
# Standard Neal (2003) slice sampler. Given current y0 and target ell(.),
# draws a new y from p(y) ∝ exp(ell(y)) by:
#   1. log_u = ell(y0) − Exponential(1)
#   2. Stepping-out: expand interval [L, R] = [y0 − w·u', y0 + w·(1−u')] where
#      u' ~ Uniform(0,1), until ell(L) < log_u and ell(R) < log_u (capped by
#      max_steps in each direction).
#   3. Shrinkage: draw y* ~ Uniform(L, R); if ell(y*) >= log_u accept, else
#      shrink to [y*, R] or [L, y*] depending on side and retry.
slice_step <- function(y0, ell_fn, w = 1.0, max_step = 25, max_shrink = 200) {
  ell0 <- ell_fn(y0)
  if (!is.finite(ell0)) return(y0)
  log_u <- ell0 - stats::rexp(1)
  u1 <- stats::runif(1)
  L <- y0 - w * u1
  R <- L + w
  for (k in seq_len(max_step)) {
    if (ell_fn(L) < log_u) break
    L <- L - w
  }
  for (k in seq_len(max_step)) {
    if (ell_fn(R) < log_u) break
    R <- R + w
  }
  for (j in seq_len(max_shrink)) {
    y_star <- stats::runif(1, L, R)
    if (ell_fn(y_star) >= log_u) return(y_star)
    if (y_star < y0) L <- y_star else R <- y_star
  }
  y0  # very rare fallback
}

# ---------- Single MCMC step for one gene --------------------------------
mcmc_step_one_gene_slice <- function(alpha_cur, gene_stats, log_pi, hyper,
                                      trunc, slice_w = 1.0) {
  H <- length(gene_stats)
  log_w <- numeric(H)
  for (h in seq_len(H)) {
    log_w[h] <- log_pi[h] + log_m_ih_one_gene(
      alpha_cur, gene_stats[[h]]$n_c,
      gene_stats[[h]]$S, gene_stats[[h]]$sumlogL,
      hyper$alpha0, hyper$nu)
  }
  M <- max(log_w)
  prob <- exp(log_w - M); prob <- prob / sum(prob)
  z_new <- sample.int(H, 1, prob = prob)

  st <- gene_stats[[z_new]]
  N_i <- sum(st$n_c)
  K_h <- length(st$n_c)

  # PIG and Damsleth auxiliaries given alpha_cur
  tau <- stats::rgamma(K_h, shape = hyper$alpha0 + st$n_c * alpha_cur, rate = 1)
  omega <- rPIG_ERGamma(N_i, alpha_cur, trunc)
  Omega <- sum(omega)
  B <- N_i * gamma_const + sum(st$sumlogL) -
       hyper$b_alpha / hyper$nu_alpha +
       sum(st$n_c * log(pmax(tau, .Machine$double.eps)))

  # Slice update on y = log alpha
  ell_fn <- function(y) {
    ell_y(y, Omega, B, st$n_c, st$S,
          hyper$alpha0, hyper$nu, N_i, hyper$b_alpha)
  }
  y_new <- slice_step(log(alpha_cur), ell_fn, w = slice_w)
  alpha_new <- exp(y_new)
  if (!is.finite(alpha_new) || alpha_new <= 0) alpha_new <- alpha_cur

  # Rao-Blackwellised conditional pattern prob given alpha_new (for logging)
  log_w2 <- numeric(H)
  for (h in seq_len(H)) {
    log_w2[h] <- log_pi[h] + log_m_ih_one_gene(
      alpha_new, gene_stats[[h]]$n_c,
      gene_stats[[h]]$S, gene_stats[[h]]$sumlogL,
      hyper$alpha0, hyper$nu)
  }
  M2 <- max(log_w2)
  pp_cond <- exp(log_w2 - M2); pp_cond <- pp_cond / sum(pp_cond)

  list(alpha = alpha_new, z = z_new, pp_cond = pp_cond)
}

gene_stats_view <- function(stats, i) {
  lapply(stats, function(st) {
    list(n_c = st$n_c, S = st$S[i, ], sumlogL = st$sumlogL[i, ])
  })
}

# ---------- Driver ---------------------------------------------------------
pig_gaga_slice_mcmc <- function(X, group, patterns, hyper, pi_v,
                                n_iter = 200, n_burn = 50, trunc = 30,
                                slice_w = 1.0,
                                init_alpha = NULL,
                                init_em_iter = 12,
                                checkpoint_file = NULL,
                                checkpoint_every = 10,
                                verbose = TRUE) {
  G <- nrow(X); H <- nrow(patterns)
  log_pi <- log(pi_v)
  stats <- gene_cluster_stats_pergene(X, group, patterns)

  # Method-of-moments init for alpha (one per gene, using EE pattern).
  if (is.null(init_alpha)) {
    if (verbose) cat("[init] PIG-EM warm start ...\n")
    init_alpha <- numeric(G)
    init_a <- max(0.5, min(hyper$nu_alpha, 50))
    st_ee <- stats[[1]]
    for (i in seq_len(G)) {
      init_alpha[i] <- pig_em_warmstart(
        st_ee$n_c, st_ee$S[i, ], st_ee$sumlogL[i, ],
        hyper$alpha0, hyper$nu, hyper$b_alpha, hyper$nu_alpha,
        init = init_a, n_iter = init_em_iter)
    }
  }
  alpha_cur <- as.numeric(init_alpha)
  z_cur <- rep(1L, G)

  pp_sum <- matrix(0, G, H)        # Rao-Blackwellised
  z_count <- matrix(0L, G, H)      # hard counts (diagnostic)
  alpha_sum <- numeric(G)
  alpha_sq_sum <- numeric(G)
  it_start <- 1L

  if (!is.null(checkpoint_file) && file.exists(checkpoint_file)) {
    ck <- readRDS(checkpoint_file)
    alpha_cur <- ck$alpha_cur; z_cur <- ck$z_cur
    pp_sum <- ck$pp_sum; z_count <- ck$z_count
    alpha_sum <- ck$alpha_sum; alpha_sq_sum <- ck$alpha_sq_sum
    it_start <- ck$it_done + 1L
    if (verbose) cat(sprintf("[ckpt] resumed from iter %d\n", it_start - 1L))
  }

  t0 <- Sys.time()
  for (it in seq.int(it_start, n_iter)) {
    for (i in seq_len(G)) {
      gv <- gene_stats_view(stats, i)
      out <- mcmc_step_one_gene_slice(alpha_cur[i], gv, log_pi, hyper,
                                      trunc, slice_w)
      alpha_cur[i] <- out$alpha
      z_cur[i] <- out$z
      if (it > n_burn) {
        pp_sum[i, ] <- pp_sum[i, ] + out$pp_cond
        z_count[i, z_cur[i]] <- z_count[i, z_cur[i]] + 1L
        alpha_sum[i] <- alpha_sum[i] + alpha_cur[i]
        alpha_sq_sum[i] <- alpha_sq_sum[i] + alpha_cur[i]^2
      }
    }
    if (verbose) {
      el <- as.numeric(Sys.time() - t0, units = "secs")
      cat(sprintf("[slice] iter %d/%d  median alpha=%.2f  cum=%.1fs\n",
                  it, n_iter, median(alpha_cur), el))
    }
    if (!is.null(checkpoint_file) && it %% checkpoint_every == 0) {
      saveRDS(list(alpha_cur = alpha_cur, z_cur = z_cur,
                   pp_sum = pp_sum, z_count = z_count,
                   alpha_sum = alpha_sum, alpha_sq_sum = alpha_sq_sum,
                   it_done = it),
              checkpoint_file)
    }
  }
  T_eff <- max(1, n_iter - n_burn)
  pp <- pp_sum / T_eff
  pp_hard <- z_count / T_eff
  alpha_mean <- alpha_sum / T_eff
  alpha_var <- pmax(alpha_sq_sum / T_eff - alpha_mean^2, 0)
  list(pp = pp, pp_hard = pp_hard,
       alpha_mean = alpha_mean, alpha_sd = sqrt(alpha_var),
       n_iter = n_iter, n_burn = n_burn,
       hyper = hyper, pi = pi_v)
}

# warm-start EM (PIG paper §3.3 closed-form maximizer with linearisation)
pig_em_warmstart <- function(n_c, S_c, sumlogL_c, alpha0, nu,
                             b_alpha, nu_alpha,
                             init = NULL, n_iter = 12, tol = 1e-4) {
  N <- sum(n_c); K <- length(n_c)
  rate0 <- alpha0 / nu
  alpha <- if (is.null(init)) max(nu_alpha, 0.5) else init
  for (it in seq_len(n_iter)) {
    E_omega_sum <- N * (digamma(1 + alpha) - digamma(1)) / (2 * alpha)
    E_logtau    <- digamma(alpha0 + alpha * n_c)
    awkward_slope <- N * (log(alpha) + 1)
    for (c_idx in seq_len(K)) {
      A <- alpha0 + alpha * n_c[c_idx]
      B <- rate0 + alpha * S_c[c_idx]
      awkward_slope <- awkward_slope +
        (-n_c[c_idx] * log(B) - A * S_c[c_idx] / B)
    }
    kappa_1 <- E_omega_sum
    kappa_2 <- N * (-digamma(1)) + sum(sumlogL_c) +
               sum(n_c * E_logtau) - (b_alpha / nu_alpha) + awkward_slope
    p <- b_alpha; M <- N
    disc <- kappa_2^2 + 8 * kappa_1 * (p + M - 1)
    if (disc < 0 || kappa_1 <= 0) break
    alpha_new <- (kappa_2 + sqrt(disc)) / (4 * kappa_1)
    if (!is.finite(alpha_new) || alpha_new <= 0) break
    if (abs(alpha_new - alpha) < tol * max(1, alpha)) {
      alpha <- alpha_new; break
    }
    alpha <- alpha_new
  }
  alpha
}

# DE call helper (Bayesian FDR)
fdr_calls_pp <- function(pp, fdr = 0.05) {
  null_prob <- pp[, 1]
  ord <- order(null_prob)
  cum_fdr <- cumsum(null_prob[ord]) / seq_along(ord)
  cutoff <- max(c(0, which(cum_fdr <= fdr)))
  calls <- rep(FALSE, nrow(pp))
  if (cutoff > 0) calls[ord[seq_len(cutoff)]] <- TRUE
  calls
}
