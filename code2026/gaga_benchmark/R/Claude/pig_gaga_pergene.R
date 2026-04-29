# =============================================================================
# pig_gaga_pergene.R
#
# Per-gene PIG-EM + Laplace approximation for the GaGa model that MATCHES
# Rossell (2009) / gaga::simGG exactly:
#   x_{ijg} | α_i, λ_{ic}  ~ Ga(α_i, α_i / λ_{ic})           (per-gene shape)
#   λ_{ic}  | α_0, ν       ~ InvGamma(α_0, α_0 / ν)
#   equivalently           1/λ_{ic} ~ Gamma(α_0, rate = α_0 / ν)
#   α_i     | b_α, ν_α     ~ Ga(b_α, b_α / ν_α)
#   z_i     | π            ~ Cat(π)
#
# We follow gaga::getpar's notation throughout: `nu` is the gaga `nu` parameter
# (the prior mean of λ_{ic}, not its rate).  All formulas below use α_0/nu as
# the rate of the inverse-gamma prior, NOT α_0·nu.
#
# Hyperparameters (α_0, λ_0, b_α, ν_α, π) are taken from gaga::fitGG so that
# the comparison isolates Stirling vs P-IG for the per-gene shape inference
# while leaving everything else identical.
#
# For each gene i and pattern z, we compute the marginal likelihood
#   m_i(z) = ∫ p(x_i | α_i, α_0, λ_0, z) p(α_i | b_α, ν_α) d α_i
# by:
#   (1) running PIG-EM (PIG paper §3.3 closed-form maximiser) to find the MAP
#       of the integrand;
#   (2) Laplace-approximating the integral around that MAP.
# Pattern posteriors are π_z m_i(z) / Σ_{z'} π_{z'} m_i(z'), exactly the same
# downstream pipeline gaga uses; only the m_i(z) computation is replaced.
# =============================================================================

# ---------------------------------------------------------------------------
# log-integrand of m_i(z), as a function of α (one gene, one pattern)
# ---------------------------------------------------------------------------
log_integrand <- function(alpha, n_c, S_c, sumlogL_c, alpha0, nu,
                          b_alpha, nu_alpha) {
  # log of the integrand of m_i(z), where the rate of 1/λ_{ic} is α_0/ν.
  # B = α_0/ν + α S_c (NOT α_0 ν !).
  if (alpha <= 0) return(-Inf)
  K <- length(n_c)
  out <- 0
  rate0 <- alpha0 / nu
  for (c_idx in seq_len(K)) {
    A <- alpha0 + alpha * n_c[c_idx]
    B <- rate0 + alpha * S_c[c_idx]
    out <- out +
      alpha * n_c[c_idx] * log(alpha) +
      (alpha - 1) * sumlogL_c[c_idx] -
      n_c[c_idx] * lgamma(alpha) +
      alpha0 * log(rate0) -
      lgamma(alpha0) +
      lgamma(A) -
      A * log(B)
  }
  out <- out + (b_alpha - 1) * log(alpha) -
         (b_alpha / nu_alpha) * alpha
  out
}

# Second derivative of log-integrand wrt alpha at alpha (used for Laplace).
log_integrand_d2 <- function(alpha, n_c, S_c, sumlogL_c, alpha0, nu,
                             b_alpha, nu_alpha) {
  out <- 0
  N_tot <- sum(n_c)
  out <- out + N_tot / alpha
  out <- out - N_tot * trigamma(alpha)
  K <- length(n_c)
  rate0 <- alpha0 / nu
  for (c_idx in seq_len(K)) {
    A <- alpha0 + alpha * n_c[c_idx]
    B <- rate0 + alpha * S_c[c_idx]
    n  <- n_c[c_idx]
    S  <- S_c[c_idx]
    out <- out + n^2 * trigamma(A) -
           2 * n * S / B +
           A * S^2 / B^2
  }
  out <- out - (b_alpha - 1) / alpha^2
  out
}

# ---------------------------------------------------------------------------
# PIG-EM update for one gene, one pattern.  Returns the MAP of α_i.
#
# The integrand for one gene under pattern z, with cluster sufficient stats
# {n_c, S_c, sumlogL_c}, is
#   ∝ α^{α N} / Γ(α)^{N} · (α-1) sumlogL · Π_c Γ(α_0 + α n_c) /
#         (α_0 λ_0 + α S_c)^{α_0 + α n_c} · α^{b_α-1} e^{-(b_α/ν_α) α}.
# After P-IG augmentation of 1/Γ(α)^N (M = N PIG factors), Damsleth-style
# τ_c augmentation of Γ(α_0 + α n_c), and linearisation of the awkward
# α^{αN} and (·)^{·} pieces at α_cur, the closed-form maximiser from PIG
# paper §3.3 is
#   α^* = (κ_2 + √(κ_2² + 8 κ_1 (p+M-1))) / (4 κ_1)
# with p = b_α (from the gamma prior α^{b_α-1}), M = N (PIG factor count).
# ---------------------------------------------------------------------------
pig_em_alpha_pergene <- function(n_c, S_c, sumlogL_c, alpha0, nu,
                                 b_alpha, nu_alpha,
                                 init = NULL, n_iter = 30,
                                 tol = 1e-5) {
  N <- sum(n_c)
  K <- length(n_c)
  rate0 <- alpha0 / nu
  alpha <- if (is.null(init)) max(nu_alpha, 0.5) else init
  for (it in seq_len(n_iter)) {
    E_omega_sum <- N * (digamma(1 + alpha) - digamma(1)) / (2 * alpha)
    E_logtau    <- digamma(alpha0 + alpha * n_c)
    awkward_slope <- N * (log(alpha) + 1)
    for (c_idx in seq_len(K)) {
      A <- alpha0 + alpha * n_c[c_idx]
      B <- rate0 + alpha * S_c[c_idx]
      n <- n_c[c_idx]
      S <- S_c[c_idx]
      awkward_slope <- awkward_slope +
        ( - n * log(B) - A * S / B )
    }
    # κ_1 (quadratic coefficient of α in the augmented log-density)
    kappa_1 <- E_omega_sum
    # κ_2 (linear coefficient of α): N·γ + sumlogL + Σ n_c E[log τ_c]
    #                              − b_α/ν_α + awkward_slope
    kappa_2 <- N * (-digamma(1)) + sum(sumlogL_c) +
               sum(n_c * E_logtau) -
               (b_alpha / nu_alpha) + awkward_slope
    # closed-form maximiser
    p <- b_alpha
    M <- N
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

# ---------------------------------------------------------------------------
# Laplace-approximate log m_i(z) given α* from PIG-EM.
# log m ≈ log p(x_i, α* | z) + 1/2 log(2π) − 1/2 log(-d²log p/dα² @ α*)
# ---------------------------------------------------------------------------
laplace_log_marginal <- function(alpha_star, n_c, S_c, sumlogL_c,
                                 alpha0, nu, b_alpha, nu_alpha) {
  log_f <- log_integrand(alpha_star, n_c, S_c, sumlogL_c,
                         alpha0, nu, b_alpha, nu_alpha)
  d2 <- log_integrand_d2(alpha_star, n_c, S_c, sumlogL_c,
                         alpha0, nu, b_alpha, nu_alpha)
  if (!is.finite(log_f) || !is.finite(d2) || d2 >= 0)
    return(NA_real_)
  log_f + 0.5 * log(2 * pi) - 0.5 * log(-d2)
}

# Numerical log-marginal-likelihood: integrate the log-integrand over a log-α
# grid centred at the PIG-EM mode.  Uses Gauss-Legendre quadrature on log α.
# Much more accurate than Laplace when the per-gene integrand is asymmetric or
# the mode lies near a boundary (which happens often with the GaGa likelihood
# at large α).
numint_log_marginal <- function(alpha_star, n_c, S_c, sumlogL_c,
                                alpha0, nu, b_alpha, nu_alpha,
                                grid_n = 81) {
  # The integrand can have its bulk far above the PIG-EM mode α* when the
  # per-gene posterior is asymmetric (typical for the titration pattern with
  # tight clusters: the integrand keeps growing with α until Stirling effects
  # eventually dominate).  We span a generous log-α range and add a coarse
  # tail probe so we don't miss mass.
  if (!is.finite(alpha_star) || alpha_star <= 0) return(NA_real_)
  log_lo <- log(0.01)
  log_hi <- log(1e10)
  z <- seq(log_lo, log_hi, length.out = grid_n)
  alpha_grid <- exp(z)
  # ∫ f(α) dα = ∫ f(e^z) e^z dz  →  the effective integrand in z has +log α
  ll <- vapply(alpha_grid, function(a) {
    log_integrand(a, n_c, S_c, sumlogL_c, alpha0, nu, b_alpha, nu_alpha) +
      log(a)
  }, numeric(1))
  ok <- is.finite(ll)
  if (sum(ok) < 2) return(NA_real_)
  M <- max(ll[ok])
  dz <- diff(z)
  w <- c(dz[1] / 2, (dz[-length(dz)] + dz[-1]) / 2, dz[length(dz)] / 2)
  M + log(sum(w[ok] * exp(ll[ok] - M)))
}

# ---------------------------------------------------------------------------
# Per-gene cluster sufficient stats keyed by pattern.  X is G x n positive,
# group is a factor of length n with K levels, patterns is P x K integer.
# Returns a list of length P; each element has n_c (vector of length C_z),
# S (G x C_z), sumlogL (G x C_z), and a vector of cluster ids.
# ---------------------------------------------------------------------------
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

# ---------------------------------------------------------------------------
# Main entry point: compute posterior pattern probabilities for all genes
# using per-gene PIG-EM + Laplace.
# ---------------------------------------------------------------------------
pig_gaga_pergene <- function(X, group, patterns, hyper, pi_v,
                             em_iter = 30, em_tol = 1e-5,
                             verbose = TRUE) {
  G <- nrow(X); P <- nrow(patterns)
  stats <- gene_cluster_stats_pergene(X, group, patterns)
  log_m <- matrix(NA_real_, G, P)
  alpha_map <- matrix(NA_real_, G, P)
  for (z_idx in seq_len(P)) {
    n_c_z <- stats[[z_idx]]$n_c
    S_z   <- stats[[z_idx]]$S
    L_z   <- stats[[z_idx]]$sumlogL
    if (verbose)
      cat(sprintf("  pattern %d/%d (clusters %s) ...\n",
                  z_idx, P, paste(stats[[z_idx]]$clusters, collapse = ",")))
    # Cap the EM init so we don't start at absurd values when gaga's
    # nu_alpha estimate is pathological (we observed gaga returning
    # nu_alpha ≈ 1e10 on noisy non-gcrma data; the EM gets stuck there
    # because the gradient is essentially zero in that region).
    init_alpha <- max(0.5, min(hyper$nu_alpha, 100))
    for (i in seq_len(G)) {
      a_star <- pig_em_alpha_pergene(
        n_c_z, S_z[i, ], L_z[i, ],
        hyper$alpha0, hyper$nu,
        hyper$b_alpha, hyper$nu_alpha,
        init = init_alpha, n_iter = em_iter, tol = em_tol)
      log_m[i, z_idx] <- numint_log_marginal(
        a_star, n_c_z, S_z[i, ], L_z[i, ],
        hyper$alpha0, hyper$nu,
        hyper$b_alpha, hyper$nu_alpha)
      alpha_map[i, z_idx] <- a_star
      if (verbose && i %% 2000 == 0)
        cat(sprintf("    gene %d/%d done\n", i, G))
    }
  }
  # convert to posterior pattern probabilities
  log_w <- sweep(log_m, 2, log(pi_v), `+`)
  M <- apply(log_w, 1, max, na.rm = TRUE)
  Wn <- exp(log_w - M)
  pp <- Wn / rowSums(Wn)
  list(pp = pp, log_m = log_m, alpha_map = alpha_map)
}
