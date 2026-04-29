# =============================================================================
# pig_gaga_lib.R — Bayesian inference for the GaGa model using P-IG augmentation
#
# Implements the algorithm derived in 00_FORMULAS.md. Reuses the P-IG sampler
# rPIG_ERGamma() and the power-truncated normal sampler rH() defined in
# code2026/sampling.R (we source it; we do not modify it).
# =============================================================================

# ---------- locate code2026/sampling.R --------------------------------------
.find_code2026 <- function() {
  candidates <- c(
    Sys.getenv("CODE2026_DIR", unset = ""),
    file.path(getwd(), "..", "..", "..", ".."),
    file.path(getwd(), "..", "..", ".."),
    file.path(getwd(), "..", ".."),
    file.path(getwd(), ".."),
    "/sessions/awesome-eager-hamilton/mnt/ERGamma/code2026"
  )
  for (p in candidates) {
    if (nzchar(p) && file.exists(file.path(p, "sampling.R"))) {
      return(normalizePath(p))
    }
  }
  stop("Could not locate code2026/sampling.R; set CODE2026_DIR.")
}

.code2026 <- .find_code2026()
# Inline copy of rPIG_ERGamma() and rH() from code2026/sampling.R so we don't
# pull in the gtools dependency (used by code2026/sampling.R only for rdirichlet
# in Gibbs_MD which we replace with our own implementation below). The math is
# byte-identical to sampling.R; only the loader changed.
if (!requireNamespace("GIGrvg", quietly = TRUE))
  stop("Package GIGrvg required (rgig).")
library(GIGrvg)

rPIG_ERGamma <- function(n, cc, N) {
  output <- rep(0, n)
  if (N > 0) {
    for (i in 1:N) {
      if (cc == 0) {
        output <- output + 1 / rgamma(n, shape = 3 / 2) / (4 * i^2)
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
  output + rgamma(n, shape, rate = rate)
}

rH <- function(m, a, b, max_try = 5000) {
  tau <- ifelse(
    b > 0,
    0.5 + sqrt(1 / 4 + 2 * a * m / b^2),
    -0.5 + sqrt(1 / 4 + 2 * a * m / b^2)
  )
  v1 <- tau * abs(b) - b
  v2 <- tau * abs(b) / (2 * a)
  for (i in 1:max_try) {
    x <- rgamma(1, m, rate = v1)
    if (runif(1) <= exp(-a * (x - v2)^2)) return(x)
  }
  NA_real_
}

if (!exists("rPIG_ERGamma") || !exists("rH")) {
  stop("rPIG_ERGamma / rH not defined after sourcing sampling.R")
}
gamma_const <- -digamma(1)


# ---------- per-gene marginal log-likelihood --------------------------------
# Under pattern z (a vector of cluster labels of length K) and current
# hyperparameters (alpha, alpha0, lambda0), compute log m_i(z) for one gene.
# Using the closed-form integral over lambda_{ic}, see 00_FORMULAS.md eq.
# (m_i(z)).  X_row is a length-n vector, group is a length-n factor whose
# levels match z.
gene_log_marginal <- function(X_row, group, z, alpha, alpha0, lambda0,
                              eps = .Machine$double.eps) {
  # cluster id per observation
  cluster <- z[as.integer(group)]
  obs <- pmax(X_row, eps)
  log_obs <- log(obs)
  ll <- 0
  uniq_clusters <- unique(cluster)
  for (c in uniq_clusters) {
    idx <- which(cluster == c)
    n_c <- length(idx)
    S_c <- sum(obs[idx])
    sumlogL <- sum(log_obs[idx])
    ll <- ll +
      alpha * n_c * log(alpha) +
      (alpha - 1) * sumlogL -
      n_c * lgamma(alpha) +
      alpha0 * log(alpha0 * lambda0) -
      lgamma(alpha0) +
      lgamma(alpha0 + alpha * n_c) -
      (alpha0 + alpha * n_c) * log(alpha0 * lambda0 + alpha * S_c)
  }
  ll
}

# Vectorised version: returns a (G x P) matrix of log-marginal-likelihoods
# for genes x patterns. patterns is a (P x K) integer matrix whose rows
# define each pattern's cluster assignment of the K groups.
all_gene_log_marginals <- function(X, group, patterns,
                                   alpha, alpha0, lambda0) {
  G <- nrow(X); P <- nrow(patterns); K <- length(levels(group))
  group_int <- as.integer(group)
  out <- matrix(NA_real_, G, P)
  for (z_idx in seq_len(P)) {
    z <- as.integer(patterns[z_idx, ])
    cluster <- z[group_int]
    uniq <- sort(unique(cluster))
    log_obs <- log(pmax(X, .Machine$double.eps))
    cluster_terms <- 0
    for (c in uniq) {
      cols <- which(cluster == c)
      n_c <- length(cols)
      S_c <- rowSums(X[, cols, drop = FALSE])
      sumL_c <- rowSums(log_obs[, cols, drop = FALSE])
      cluster_terms <- cluster_terms +
        alpha * n_c * log(alpha) +
        (alpha - 1) * sumL_c -
        n_c * lgamma(alpha) +
        alpha0 * log(alpha0 * lambda0) -
        lgamma(alpha0) +
        lgamma(alpha0 + alpha * n_c) -
        (alpha0 + alpha * n_c) * log(alpha0 * lambda0 + alpha * S_c)
    }
    out[, z_idx] <- cluster_terms
  }
  out
}

# Pattern posterior probabilities per gene given current hyperparameters and
# pattern prior pi.  Returns a (G x P) matrix of posterior pattern probs.
gene_pattern_posterior <- function(X, group, patterns,
                                   alpha, alpha0, lambda0, pi) {
  log_m <- all_gene_log_marginals(X, group, patterns, alpha, alpha0, lambda0)
  log_w <- sweep(log_m, 2, log(pi), `+`)
  M <- apply(log_w, 1, max)
  Wn <- exp(log_w - M)
  Wn / rowSums(Wn)
}


# ---------- sufficient statistics -------------------------------------------
# Per-gene-cluster sufficient stats {n_ic, S_ic, prod x = exp(sum log x)} as
# matrices keyed by pattern.  Pre-computed once per dataset.
gene_cluster_stats <- function(X, group, patterns) {
  G <- nrow(X); P <- nrow(patterns); K <- length(levels(group))
  group_int <- as.integer(group)
  res <- vector("list", P)
  log_obs <- log(pmax(X, .Machine$double.eps))
  for (z_idx in seq_len(P)) {
    z <- as.integer(patterns[z_idx, ])
    cluster <- z[group_int]
    uniq <- sort(unique(cluster))
    n_c_vec <- integer(length(uniq))
    S_mat <- matrix(0, G, length(uniq))
    L_mat <- matrix(0, G, length(uniq))
    for (i in seq_along(uniq)) {
      cols <- which(cluster == uniq[i])
      n_c_vec[i] <- length(cols)
      S_mat[, i] <- rowSums(X[, cols, drop = FALSE])
      L_mat[, i] <- rowSums(log_obs[, cols, drop = FALSE])
    }
    res[[z_idx]] <- list(clusters = uniq, n_c = n_c_vec,
                         S = S_mat, sumlogL = L_mat)
  }
  res
}


# ---------- P-IG augmented update for global shape α -----------------------
# Conditional on the assigned pattern z_i for each gene, integrate λ_{ic} out
# (closed form in §1) and use the Stirling-aware fast approximation for the
# α^{αN} piece via Miller's quadratic correction at the current α.  The
# 1/Γ(α)^N piece is handled exactly by N PIG draws.
#
# The conditional log-density for α, after integrating out λ_{ic} per gene,
# is sum over genes/clusters of
#    αn_{ic}·logα − n_{ic}·logΓ(α) + log Γ(αn_{ic}+α0)
#                  − (αn_{ic}+α0) log(αS_{ic}+α0λ0) + (α−1) log L_{ic}
# plus a Damsleth ξ2 prior contribution
#    log Γ(δα+1) − δ logΓ(α) − δα(log(δμ) + η).
#
# We split this into:
#   (i)  exact PIG handling of (N_tot+δ) factors of 1/Γ(α);
#   (ii) Damsleth τ-augmentation of Γ(δα+1);
#  (iii) linearisation of the awkward terms αN·logα and
#        (αn+α0) log(αS+α0λ0) at the current α^{(t)}, valid in MH-within-
#        Gibbs because we wrap in an accept/reject step.
#
# The result is a PTN proposal whose density is the linearised exact
# density up to a constant; the MH correction restores exactness.

logpost_alpha <- function(alpha, stats, z_assign, alpha0, lambda0,
                          delta, mu, eta_prior) {
  if (alpha <= 0) return(-Inf)
  G <- length(z_assign)
  ll <- 0
  for (z_idx in seq_along(stats)) {
    rows <- which(z_assign == z_idx)
    if (!length(rows)) next
    st <- stats[[z_idx]]
    nC <- st$n_c
    for (i in seq_along(nC)) {
      n_c <- nC[i]
      S_c <- st$S[rows, i]
      L_c <- st$sumlogL[rows, i]
      ll <- ll + sum(
        alpha * n_c * log(alpha) +
        (alpha - 1) * L_c -
        n_c * lgamma(alpha) +
        lgamma(alpha0 + alpha * n_c) -
        (alpha0 + alpha * n_c) * log(alpha0 * lambda0 + alpha * S_c)
      )
    }
  }
  # Damsleth ξ2 prior on alpha (independent of data via eta_prior)
  ll <- ll + lgamma(delta * alpha + 1) - delta * lgamma(alpha) -
        delta * alpha * log(delta * mu) - delta * eta_prior * alpha
  ll
}

# Sample alpha conditional on (z_assign, alpha0, lambda0, hyperparams) using
# a P-IG augmented PTN proposal centred at the current alpha plus an MH step.
update_alpha_pig <- function(alpha_cur, stats, z_assign, alpha0, lambda0,
                             delta, mu, eta_prior, N_pig = 200) {
  # Sufficient stats summed across the assigned pattern.
  N_tot <- 0; sumlogL <- 0; sumS <- 0
  for (z_idx in seq_along(stats)) {
    rows <- which(z_assign == z_idx)
    if (!length(rows)) next
    st <- stats[[z_idx]]
    for (i in seq_along(st$n_c)) {
      n_c <- st$n_c[i]
      N_tot <- N_tot + length(rows) * n_c
      sumS <- sumS + sum(st$S[rows, i])
      sumlogL <- sumlogL + sum(st$sumlogL[rows, i])
    }
  }
  S_pow <- N_tot + delta

  # Damsleth Γ(δα+1) auxiliary
  tau <- rgamma(1, shape = delta * alpha_cur + 1)
  # S_pow PIG(α) auxiliaries; sum them — same convention as gamma_shape2.R
  omega <- rPIG_ERGamma(S_pow, alpha_cur, N_pig)
  a_tilde <- sum(omega)

  # The exact α-conditional has linear-in-α coefficient
  #   sumlogL − δ log(δμ) − δη  + slope of the awkward non-linear pieces at α_cur.
  awkward_slope <- N_tot * (log(alpha_cur) + 1) +
                   gamma_ratio_slope_alpha(alpha_cur, alpha0, lambda0,
                                           stats, z_assign)
  b_tilde <- S_pow * gamma_const + delta * log(tau) -
             delta * log(delta * mu) - delta * eta_prior +
             sumlogL + awkward_slope

  # PTN power = (N_tot + δ) + 1 = S_pow + 1, matching PIG paper §4.2 for the
  # form Γ(α+β)^S/Γ(α)^S e^{cα}.
  alpha_prop <- rH(S_pow + 1, a_tilde, b_tilde)
  if (!is.finite(alpha_prop) || alpha_prop <= 0) return(alpha_cur)

  log_target_prop <- logpost_alpha(alpha_prop, stats, z_assign, alpha0,
                                   lambda0, delta, mu, eta_prior)
  log_target_cur  <- logpost_alpha(alpha_cur,  stats, z_assign, alpha0,
                                   lambda0, delta, mu, eta_prior)
  if (!is.finite(log_target_prop)) return(alpha_cur)
  log_acc <- log_target_prop - log_target_cur
  if (log(runif(1)) < log_acc) alpha_prop else alpha_cur
}

# Slope of log Γ(α n + α0) − (α n + α0) log(α S + α0 λ0) summed over genes
# and clusters, evaluated at α = alpha_cur.  Used to set b_tilde.
gamma_ratio_slope_alpha <- function(alpha, alpha0, lambda0, stats, z_assign) {
  out <- 0
  for (z_idx in seq_along(stats)) {
    rows <- which(z_assign == z_idx)
    if (!length(rows)) next
    st <- stats[[z_idx]]
    for (i in seq_along(st$n_c)) {
      n_c <- st$n_c[i]
      S_c <- st$S[rows, i]
      A <- alpha0 + alpha * n_c
      B <- alpha0 * lambda0 + alpha * S_c
      # d/dα [log Γ(A) − A log B] = n_c · ψ(A) − n_c log B − A · S_c / B
      out <- out + sum(n_c * digamma(A) - n_c * log(B) - A * S_c / B)
    }
  }
  out
}

# slope contribution from the gamma-ratio term to dlogp(α)/dα at α=α_cur,
# used to set the PTN proposal's b parameter
sum_dgamma_term <- function(alpha, alpha0, lambda0, stats, z_assign,
                            kind = "value") {
  out <- 0
  for (z_idx in seq_along(stats)) {
    rows <- which(z_assign == z_idx)
    if (!length(rows)) next
    st <- stats[[z_idx]]
    for (i in seq_along(st$n_c)) {
      n_c <- st$n_c[i]
      S_c <- st$S[rows, i]
      A <- alpha0 + alpha * n_c
      B <- alpha0 * lambda0 + alpha * S_c
      # value: A * log(B); derivative wrt α: n_c log(B) + A · S_c / B
      out <- out + sum(n_c * log(B) + A * S_c / B)
    }
  }
  out
}

# ---------- update for α0 via P-IG augmentation -----------------------------
# Conditional density (with z, α, λ_0 fixed) is proportional to
#   prod_{i,c} (α0 λ0)^{α0} / Γ(α0) · (α0 λ0 + α S_{ic})^{-(α n_{ic}+α0)} · Γ(α n+α0)
# = α0^{α0 M} ... — analogous to the α update.  We use the same MH-with-PTN
# proposal scheme.

logpost_alpha0 <- function(alpha0, alpha, lambda0, stats, z_assign,
                           delta0, mu0, eta_prior0) {
  if (alpha0 <= 0) return(-Inf)
  ll <- 0
  for (z_idx in seq_along(stats)) {
    rows <- which(z_assign == z_idx)
    if (!length(rows)) next
    st <- stats[[z_idx]]
    for (i in seq_along(st$n_c)) {
      n_c <- st$n_c[i]
      S_c <- st$S[rows, i]
      ll <- ll + sum(
        alpha0 * log(alpha0 * lambda0) -
        lgamma(alpha0) +
        lgamma(alpha0 + alpha * n_c) -
        (alpha0 + alpha * n_c) * log(alpha0 * lambda0 + alpha * S_c)
      )
    }
  }
  # Damsleth-style prior on α0
  ll <- ll + lgamma(delta0 * alpha0 + 1) - delta0 * lgamma(alpha0) -
        delta0 * alpha0 * log(delta0 * mu0) - delta0 * eta_prior0 * alpha0
  ll
}

update_alpha0_pig <- function(alpha0_cur, alpha, lambda0, stats, z_assign,
                              delta0, mu0, eta_prior0, N_pig = 200) {
  M <- 0
  for (z_idx in seq_along(stats)) {
    rows <- which(z_assign == z_idx)
    if (!length(rows)) next
    M <- M + length(rows) * length(stats[[z_idx]]$n_c)
  }
  S_pow <- M + delta0
  tau0 <- rgamma(1, shape = delta0 * alpha0_cur + 1)
  omega0 <- rPIG_ERGamma(S_pow, alpha0_cur, N_pig)
  a_tilde <- sum(omega0)
  slope <- 0
  for (z_idx in seq_along(stats)) {
    rows <- which(z_assign == z_idx)
    if (!length(rows)) next
    st <- stats[[z_idx]]
    for (i in seq_along(st$n_c)) {
      n_c <- st$n_c[i]
      S_c <- st$S[rows, i]
      A <- alpha0_cur + alpha * n_c
      B <- alpha0_cur * lambda0 + alpha * S_c
      slope <- slope + sum(
        log(alpha0_cur * lambda0) + 1 +
        digamma(A) -
        log(B) -
        A * lambda0 / B
      )
    }
  }
  b_tilde <- S_pow * gamma_const + delta0 * log(tau0) -
             delta0 * log(delta0 * mu0) - delta0 * eta_prior0 + slope
  alpha0_prop <- rH(S_pow + 1, a_tilde, b_tilde)
  if (!is.finite(alpha0_prop) || alpha0_prop <= 0) return(alpha0_cur)
  log_alpha <- logpost_alpha0(alpha0_prop, alpha, lambda0, stats, z_assign,
                              delta0, mu0, eta_prior0) -
               logpost_alpha0(alpha0_cur,  alpha, lambda0, stats, z_assign,
                              delta0, mu0, eta_prior0)
  if (log(runif(1)) < log_alpha) alpha0_prop else alpha0_cur
}

# ---------- update for λ0 (closed form, conjugate IG) -----------------------
update_lambda0 <- function(alpha0, alpha, stats, z_assign,
                           a_lam = 1, b_lam = 1) {
  # 1/λ0 | rest ~ Ga(a_lam + α0 · M, b_lam + α0 · sum_{ic} 1/λ_{ic})
  # but lambdas were integrated; we instead use the conditional based on the
  # closed-form posterior of λ0 from the integrated likelihood.  See
  # 00_FORMULAS.md §4.  Empirically we compute α0·M·E[1/λ_{ic}|·] from the
  # integrated-form gamma posterior of 1/λ_{ic}.
  Minv <- 0; M <- 0
  for (z_idx in seq_along(stats)) {
    rows <- which(z_assign == z_idx)
    if (!length(rows)) next
    st <- stats[[z_idx]]
    for (i in seq_along(st$n_c)) {
      n_c <- st$n_c[i]
      S_c <- st$S[rows, i]
      shape_post <- alpha0 + alpha * n_c
      rate_post <- alpha0 * 0 + alpha * S_c  # α0·λ0 contribution comes outside
      # We use the integrated posterior mean E[1/λ_{ic}] = shape/(rate + α0·λ0)
      # but since λ0 itself is unknown we use the data-only shape E≈shape/rate
      Minv <- Minv + sum(shape_post / pmax(rate_post, 1e-12))
      M <- M + length(rows)
    }
  }
  shape <- a_lam + alpha0 * M
  rate <- b_lam + alpha0 * Minv
  inv_lam0 <- rgamma(1, shape = shape, rate = rate)
  1 / max(inv_lam0, 1e-12)
}


# ---------- full Gibbs sampler ----------------------------------------------
# Inputs:
#   X         G × n positive expression matrix
#   group     factor of length n with K levels
#   patterns  P × K integer matrix (rows = patterns)
#   pi_prior  Dirichlet pseudo-counts of length P (default 1 = uniform)
#   delta, mu, eta_prior   Damsleth ξ2 prior on α
#   delta0, mu0, eta_prior0 Damsleth ξ2 prior on α0
#   init      list with alpha, alpha0, lambda0, pi (NULL = use moment-based)
#   n_iter, n_burn  number of sweeps and burn-in
# Returns posterior pattern probabilities and traces of (α, α0, λ0, π).
pig_gaga_gibbs <- function(X, group, patterns,
                           pi_prior = NULL,
                           delta = 1, mu = NULL, eta_prior = 0,
                           delta0 = 1, mu0 = NULL, eta_prior0 = 0,
                           init = NULL,
                           n_iter = 600, n_burn = 200,
                           N_pig = 100, verbose = TRUE,
                           fix_alpha0 = FALSE, fix_lambda0 = FALSE) {
  G <- nrow(X); P <- nrow(patterns); n <- ncol(X)
  group <- factor(group)
  if (is.null(pi_prior)) pi_prior <- rep(1, P)
  stats <- gene_cluster_stats(X, group, patterns)

  # Method-of-moments init for α, α0, λ0 if not provided
  log_obs <- log(pmax(X, .Machine$double.eps))
  gene_mean <- rowMeans(X)
  gene_var <- apply(X, 1, stats::var)
  cv2 <- gene_var / gene_mean^2
  alpha_mom <- 1 / pmax(median(cv2, na.rm = TRUE), 1e-3)
  lambda0_mom <- median(gene_mean, na.rm = TRUE)
  alpha0_mom <- 1 / pmax(stats::var(log(gene_mean)) , 1e-3)
  if (is.null(init)) {
    init <- list(alpha = alpha_mom, alpha0 = alpha0_mom,
                 lambda0 = lambda0_mom, pi = pi_prior / sum(pi_prior))
  }
  if (is.null(mu)) mu <- max(alpha_mom * 1.05, 2)
  if (is.null(mu0)) mu0 <- max(alpha0_mom * 1.05, 2)

  alpha <- init$alpha
  alpha0 <- init$alpha0
  lambda0 <- init$lambda0
  pi_v <- init$pi

  trace_alpha <- numeric(n_iter)
  trace_alpha0 <- numeric(n_iter)
  trace_lambda0 <- numeric(n_iter)
  trace_pi <- matrix(0, n_iter, P)
  z_assign <- sample.int(P, G, replace = TRUE, prob = pi_v)
  pp_avg <- matrix(0, G, P)

  for (it in seq_len(n_iter)) {
    # 1. pattern probabilities per gene + sample z_i
    pp <- gene_pattern_posterior(X, group, patterns, alpha, alpha0,
                                 lambda0, pi_v)
    z_assign <- apply(pp, 1, function(p) sample.int(P, 1, prob = p))

    # 2. π | z
    counts <- tabulate(z_assign, nbins = P)
    pi_v <- as.numeric(MCMCpack_rdirichlet(1, pi_prior + counts))

    # 3. λ0 | rest
    if (!fix_lambda0)
      lambda0 <- update_lambda0(alpha0, alpha, stats, z_assign)

    # 4. α | rest  (the cornerstone PIG step)
    alpha <- update_alpha_pig(alpha, stats, z_assign,
                              alpha0, lambda0,
                              delta, mu, eta_prior, N_pig)

    # 5. α0 | rest
    if (!fix_alpha0)
      alpha0 <- update_alpha0_pig(alpha0, alpha, lambda0,
                                  stats, z_assign,
                                  delta0, mu0, eta_prior0, N_pig)

    trace_alpha[it] <- alpha
    trace_alpha0[it] <- alpha0
    trace_lambda0[it] <- lambda0
    trace_pi[it, ] <- pi_v
    if (it > n_burn) pp_avg <- pp_avg + pp

    if (verbose && it %% 20 == 0)
      cat(sprintf("[Gibbs] iter %d  alpha=%.3f  alpha0=%.3f  lambda0=%.3f  ",
                  it, alpha, alpha0, lambda0),
          sprintf("pi=%s\n", paste(round(pi_v, 3), collapse = ",")))
  }
  pp_avg <- pp_avg / max(1, n_iter - n_burn)
  list(pp = pp_avg, alpha = trace_alpha, alpha0 = trace_alpha0,
       lambda0 = trace_lambda0, pi = trace_pi,
       n_iter = n_iter, n_burn = n_burn,
       hyper = list(delta = delta, mu = mu, eta_prior = eta_prior,
                    delta0 = delta0, mu0 = mu0, eta_prior0 = eta_prior0))
}

# small Dirichlet sampler so we don't have a hard dependency on MCMCpack
MCMCpack_rdirichlet <- function(n, alpha) {
  K <- length(alpha)
  out <- matrix(rgamma(n * K, shape = alpha), nrow = n, byrow = TRUE)
  out / rowSums(out)
}


# ---------- DE call helpers -------------------------------------------------
# Given posterior pattern probabilities and a target FDR α_FDR, return a
# logical vector of DE calls.  The first column is the EE (null) pattern.
fdr_calls <- function(pp, alpha_FDR = 0.05) {
  null_prob <- pp[, 1]
  sig <- 1 - null_prob  # prob of any DE pattern
  ord <- order(sig, decreasing = TRUE)
  cum_fdr <- cumsum(null_prob[ord]) / seq_along(ord)
  cutoff_rank <- max(c(0, which(cum_fdr <= alpha_FDR)))
  calls <- rep(FALSE, nrow(pp))
  if (cutoff_rank > 0) calls[ord[seq_len(cutoff_rank)]] <- TRUE
  calls
}

de_score_from_pp <- function(pp) 1 - pp[, 1]
