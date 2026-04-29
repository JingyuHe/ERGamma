# =============================================================================
# pig_gaga_mcmc_lib.R
#
# Per-gene PIG-augmented MCMC for the GaGa model, following the formulas in
# PIG_GAGA_MCMC_FORMULAS.md.  Targets the EXACT collapsed α-conditional, no
# linearisation: PIG augmentation handles 1/Γ(α)^N, gamma augmentation handles
# Γ(a_0 + n_c α), the residual R(α) = N α log α − Σ_c (a_0 + n_c α) log(r_0 + α S_c)
# is corrected by Metropolis–Hastings.
# =============================================================================

if (!requireNamespace("GIGrvg", quietly = TRUE))
  stop("Package GIGrvg required (rgig).")
library(GIGrvg)

gamma_const <- -digamma(1)

# ---------- PIG draw (truncated GIG sum + gamma tail correction) ----------
# Identical to code2026/sampling.R::rPIG_ERGamma, inlined to avoid the
# legacy gtools dependency.  See PIG paper §4.1 (truncated GIG sum).
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

# ---------- PTN draw  rH(power, quad, lin)  -------------------------------
# Density ∝ x^{power-1} exp(-quad x^2 + lin x), x > 0.  Damsleth's gamma
# rejection sampler with proposal Gamma(power, rate = τ |lin| − lin).
rPTN <- function(power, quad, lin, max_try = 5000) {
  if (!is.finite(power) || power <= 0 || !is.finite(quad) || quad <= 0)
    return(NA_real_)
  if (!is.finite(lin)) return(NA_real_)
  if (abs(lin) < 1e-12) lin <- sign(lin) * 1e-12
  if (lin == 0) lin <- 1e-12
  tau <- if (lin > 0) {
    0.5 + sqrt(0.25 + 2 * quad * power / lin^2)
  } else {
    -0.5 + sqrt(0.25 + 2 * quad * power / lin^2)
  }
  rate_g <- tau * abs(lin) - lin
  centre <- tau * abs(lin) / (2 * quad)
  for (i in seq_len(max_try)) {
    x <- stats::rgamma(1, shape = power, rate = rate_g)
    if (stats::runif(1) <= exp(-quad * (x - centre)^2)) return(x)
  }
  NA_real_
}

# ---------- gene-cluster sufficient statistics keyed by pattern ----------
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

# ---------- log m_ih(α) for one gene given its sufficient stats ---------
log_m_ih_one_gene <- function(alpha, n_c_vec, S_vec, L_vec, alpha0, nu) {
  if (alpha <= 0) return(-Inf)
  rate0 <- alpha0 / nu
  K <- length(n_c_vec)
  out <- 0
  for (c_idx in seq_len(K)) {
    n <- n_c_vec[c_idx]
    S <- S_vec[c_idx]
    L <- L_vec[c_idx]
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

# Full augmented log-conditional log p(α | ω, τ, x, z, θ), used by the
# random-walk fallback's MH ratio.  For fixed (ω, τ), this is
# (N+b_α-1) log α − Ω α^2 + B α + R(α).
log_full_aug_cond <- function(alpha, Omega, B, n_c, S_c, alpha0, nu,
                              N_i, b_alpha) {
  if (alpha <= 0) return(-Inf)
  R <- R_alpha_one_gene(alpha, n_c, S_c, alpha0, nu)
  (N_i + b_alpha - 1) * log(alpha) - Omega * alpha^2 + B * alpha + R
}

# ---------- Single MCMC step for one gene  --------------------------------
# Hybrid kernel: try PTN-PIG proposal first (driven by the PIG augmentation),
# fall back to a log-α random walk if PTN is rejected.  Both moves preserve
# the same target distribution; the chain's stationary distribution is the
# exact joint posterior on (α, z, ω, τ).
mcmc_step_one_gene <- function(alpha_cur, gene_stats, log_pi, hyper, trunc,
                               rw_sd = 0.15) {
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

  tau <- stats::rgamma(K_h, shape = hyper$alpha0 + st$n_c * alpha_cur, rate = 1)
  omega <- rPIG_ERGamma(N_i, alpha_cur, trunc)
  Omega <- sum(omega)
  B <- N_i * gamma_const + sum(st$sumlogL) -
       hyper$b_alpha / hyper$nu_alpha +
       sum(st$n_c * log(pmax(tau, .Machine$double.eps)))

  alpha_star <- rPTN(N_i + hyper$b_alpha, Omega, B)
  ptn_ok <- is.finite(alpha_star) && alpha_star > 0
  accepted_ptn <- FALSE
  if (ptn_ok) {
    R_star <- R_alpha_one_gene(alpha_star, st$n_c, st$S,
                               hyper$alpha0, hyper$nu)
    R_cur  <- R_alpha_one_gene(alpha_cur,  st$n_c, st$S,
                               hyper$alpha0, hyper$nu)
    if (is.finite(R_star - R_cur) &&
        log(stats::runif(1)) < (R_star - R_cur)) {
      alpha_cur <- alpha_star
      accepted_ptn <- TRUE
    }
  }

  # Random-walk fallback on log α (always run, gives a backup mix-step that
  # keeps the chain moving even when the PTN proposal is far from the joint
  # mode due to R(α) curvature).  Same conditional given (ω, τ) — Jacobian
  # of α → log α adds  log α* − log α_cur  to the MH ratio.
  eps <- stats::rnorm(1, sd = rw_sd)
  alpha_rw <- alpha_cur * exp(eps)
  if (is.finite(alpha_rw) && alpha_rw > 0) {
    log_target_rw  <- log_full_aug_cond(alpha_rw, Omega, B, st$n_c, st$S,
                                        hyper$alpha0, hyper$nu,
                                        N_i, hyper$b_alpha)
    log_target_cur <- log_full_aug_cond(alpha_cur, Omega, B, st$n_c, st$S,
                                        hyper$alpha0, hyper$nu,
                                        N_i, hyper$b_alpha)
    log_q_jac <- log(alpha_rw) - log(alpha_cur)
    log_acc_rw <- log_target_rw - log_target_cur + log_q_jac
    if (is.finite(log_acc_rw) &&
        log(stats::runif(1)) < log_acc_rw) {
      alpha_cur <- alpha_rw
    }
  }

  list(alpha = alpha_cur, z = z_new, accepted = accepted_ptn, log_w = log_w)
}

# Extract the per-gene cluster sufficient statistics from the dataset-level
# stats (returned by gene_cluster_stats_pergene) for gene index i.
gene_stats_view <- function(stats, i) {
  lapply(stats, function(st) {
    list(n_c = st$n_c, S = st$S[i, ], sumlogL = st$sumlogL[i, ])
  })
}

# Internal PIG-EM (closed-form maximizer, PIG paper §3.3) used for warm-start.
# Same math as pig_em_alpha_pergene() in pig_gaga_pergene.R; included here so
# the MCMC library is self-contained.
pig_em_alpha_pergene_internal <- function(n_c, S_c, sumlogL_c, alpha0, nu,
                                          b_alpha, nu_alpha,
                                          init = NULL, n_iter = 20,
                                          tol = 1e-5) {
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

# ---------- Driver: per-gene PIG-MCMC across the whole dataset ------------
pig_gaga_mcmc <- function(X, group, patterns, hyper, pi_v,
                          n_iter = 200, n_burn = 50, trunc = 30,
                          init_alpha = NULL,
                          init_em_iter = 15,
                          rw_sd = 0.15,
                          checkpoint_file = NULL,
                          checkpoint_every = 25,
                          verbose = TRUE) {
  G <- nrow(X); H <- nrow(patterns)
  log_pi <- log(pi_v)
  stats <- gene_cluster_stats_pergene(X, group, patterns)

  if (is.null(init_alpha) && init_em_iter > 0) {
    # Initialize each gene's α at its PIG-EM mode under the most-likely
    # pattern (the EE pattern under the prior is a reasonable default; EM
    # converges quickly so we just use it as a warm start for the MCMC).
    if (verbose) cat("[init] PIG-EM warm start ...\n")
    init_alpha <- numeric(G)
    init_a <- max(0.5, min(hyper$nu_alpha, 50))
    st_ee <- stats[[1]]
    for (i in seq_len(G)) {
      init_alpha[i] <- pig_em_alpha_pergene_internal(
        st_ee$n_c, st_ee$S[i, ], st_ee$sumlogL[i, ],
        hyper$alpha0, hyper$nu, hyper$b_alpha, hyper$nu_alpha,
        init = init_a, n_iter = init_em_iter)
    }
  }
  if (is.null(init_alpha)) {
    init_alpha <- rep(max(0.5, min(hyper$nu_alpha, 50)), G)
  }
  alpha_cur <- as.numeric(init_alpha)
  z_cur <- rep(1L, G)
  pp_count <- matrix(0, G, H)
  alpha_sum <- numeric(G)
  alpha_sq_sum <- numeric(G)
  acc_count <- integer(G)
  it_start <- 1L

  # resume from checkpoint
  if (!is.null(checkpoint_file) && file.exists(checkpoint_file)) {
    ck <- readRDS(checkpoint_file)
    alpha_cur <- ck$alpha_cur; z_cur <- ck$z_cur
    pp_count <- ck$pp_count; alpha_sum <- ck$alpha_sum
    alpha_sq_sum <- ck$alpha_sq_sum; acc_count <- ck$acc_count
    it_start <- ck$it_done + 1L
    if (verbose)
      cat(sprintf("[ckpt] resumed from iter %d\n", it_start - 1L))
  }

  t0 <- Sys.time()
  for (it in seq.int(it_start, n_iter)) {
    for (i in seq_len(G)) {
      gv <- gene_stats_view(stats, i)
      out <- mcmc_step_one_gene(alpha_cur[i], gv, log_pi, hyper, trunc,
                                rw_sd = rw_sd)
      alpha_cur[i] <- out$alpha
      z_cur[i] <- out$z
      if (out$accepted) acc_count[i] <- acc_count[i] + 1L
      if (it > n_burn) {
        pp_count[i, z_cur[i]] <- pp_count[i, z_cur[i]] + 1
        alpha_sum[i] <- alpha_sum[i] + alpha_cur[i]
        alpha_sq_sum[i] <- alpha_sq_sum[i] + alpha_cur[i]^2
      }
    }
    if (verbose) {
      el <- as.numeric(Sys.time() - t0, units = "secs")
      cat(sprintf("[mcmc] iter %d/%d  acc=%.1f%%  cum=%.1fs\n",
                  it, n_iter, 100 * mean(acc_count) / it, el))
    }
    if (!is.null(checkpoint_file) && it %% checkpoint_every == 0) {
      saveRDS(list(alpha_cur = alpha_cur, z_cur = z_cur,
                   pp_count = pp_count, alpha_sum = alpha_sum,
                   alpha_sq_sum = alpha_sq_sum, acc_count = acc_count,
                   it_done = it),
              checkpoint_file)
    }
  }
  T_eff <- max(1, n_iter - n_burn)
  pp <- pp_count / T_eff
  alpha_mean <- alpha_sum / T_eff
  alpha_var  <- pmax(alpha_sq_sum / T_eff - alpha_mean^2, 0)
  list(pp = pp, alpha_mean = alpha_mean, alpha_sd = sqrt(alpha_var),
       acc_rate = acc_count / n_iter,
       n_iter = n_iter, n_burn = n_burn,
       hyper = hyper, pi = pi_v)
}

# ---------- DE call helper (Bayesian FDR) -------------------------------
fdr_calls_pp <- function(pp, fdr = 0.05) {
  null_prob <- pp[, 1]
  ord <- order(null_prob)
  cum_fdr <- cumsum(null_prob[ord]) / seq_along(ord)
  cutoff <- max(c(0, which(cum_fdr <= fdr)))
  calls <- rep(FALSE, nrow(pp))
  if (cutoff > 0) calls[ord[seq_len(cutoff)]] <- TRUE
  calls
}
