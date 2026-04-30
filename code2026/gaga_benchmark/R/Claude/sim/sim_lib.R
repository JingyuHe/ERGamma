# =============================================================================
# sim_lib.R - shared utilities for the GaGa Stirling-vs-PIG simulation study
# =============================================================================

suppressMessages({
  library(gaga)
  library(Biobase)
  library(GIGrvg)
})

# Source PIG helpers from the parent folder.
.parent <- normalizePath(file.path(dirname(sys.frame(1)$ofile), ".."),
                         mustWork = FALSE)
if (!nzchar(.parent) || !dir.exists(.parent)) {
  .parent <- "/sessions/awesome-eager-hamilton/mnt/ERGamma/code2026/gaga_benchmark/R/Claude"
}
source(file.path(.parent, "pig_gaga_pergene.R"))    # log_integrand, gene_cluster_stats_pergene
source(file.path(.parent, "pig_gaga_slice_lib.R"))  # pig_gaga_slice_mcmc

# ---------- 4 oracle scenarios + alpha sweep grid -------------------------
scenarios <- list(
  A = list(name = "Easy",    G = 200, n = 10,
           b_alpha = 20, nu_alpha = 20, a0 = 20, nu = 5, p_de = 0.20),
  B = list(name = "Medium",  G = 200, n = 5,
           b_alpha = 8,  nu_alpha = 5,  a0 = 10, nu = 5, p_de = 0.20),
  C = list(name = "Hard",    G = 200, n = 4,
           b_alpha = 4,  nu_alpha = 1.5, a0 = 5, nu = 5, p_de = 0.20),
  D = list(name = "Extreme", G = 200, n = 3,
           b_alpha = 2,  nu_alpha = 0.5, a0 = 3, nu = 5, p_de = 0.20)
)

alpha_sweep_grid <- c(0.5, 1, 2, 4, 8, 16, 32, 64)  # ν_α values

# ---------- simulate from gaga::simGG (the exact gaga model) --------------
simulate_data <- function(spec, seed) {
  set.seed(seed)
  xsim <- gaga::simGG(
    n = spec$G, m = c(spec$n, spec$n), p.de = spec$p_de,
    a0 = spec$a0, nu = spec$nu,
    balpha = spec$b_alpha, nualpha = spec$nu_alpha,
    equalcv = TRUE)
  X <- Biobase::exprs(xsim)
  groups <- factor(Biobase::pData(xsim)$group)
  fd <- Biobase::fData(xsim)
  truth <- abs(fd$mean.1 - fd$mean.2) > 1e-12
  list(X = X, groups = groups, truth = truth, sim = xsim)
}

# ---------- gold standard pp via fine 1D log-α quadrature -----------------
gold_pp <- function(X, groups, patterns, hyper, pi_v, n_grid = 600,
                    log_alpha_min = log(0.01), log_alpha_max = log(1e6)) {
  G <- nrow(X); P <- nrow(patterns)
  stats <- gene_cluster_stats_pergene(X, groups, patterns)
  log_m <- matrix(NA_real_, G, P)
  z_grid <- seq(log_alpha_min, log_alpha_max, length.out = n_grid)
  alpha_grid <- exp(z_grid)
  dz <- diff(z_grid)
  w <- c(dz[1] / 2, (dz[-length(dz)] + dz[-1]) / 2, dz[length(dz)] / 2)
  for (h in seq_len(P)) {
    n_c <- stats[[h]]$n_c
    S_mat <- stats[[h]]$S
    L_mat <- stats[[h]]$sumlogL
    for (i in seq_len(G)) {
      ll <- vapply(alpha_grid, function(a) {
        v <- log_integrand(a, n_c, S_mat[i, ], L_mat[i, ],
                           hyper$alpha0, hyper$nu,
                           hyper$b_alpha, hyper$nu_alpha) + log(a)
        if (is.finite(v)) v else -Inf
      }, numeric(1))
      good <- is.finite(ll)
      if (!any(good)) { log_m[i, h] <- NA_real_; next }
      M <- max(ll[good])
      log_m[i, h] <- M + log(sum(w[good] * exp(ll[good] - M)))
    }
  }
  log_w <- sweep(log_m, 2, log(pi_v), `+`)
  Mrow <- apply(log_w, 1, max, na.rm = TRUE)
  Wn <- exp(log_w - Mrow)
  pp <- Wn / rowSums(Wn)
  list(pp = pp, log_m = log_m)
}

# ---------- run gaga + gold + PIG-MCMC on one dataset --------------------
run_methods <- function(X, groups, patterns,
                        n_iter = 200, n_burn = 60, trunc = 30,
                        slice_w = 1.0, n_grid = 600,
                        verbose = FALSE) {
  # Step 1: gaga::fitGG -> hyperparameters from quickEM
  gg_fit <- gaga::fitGG(X, groups, patterns = patterns, equalcv = TRUE,
                        nclust = 1, method = "quickEM", trace = FALSE)
  gg_fit <- gaga::parest(gg_fit, x = X, groups = groups, alpha = 0.05)
  gg_par <- gaga::getpar(gg_fit)
  hyper <- list(alpha0 = gg_par$a0, nu = gg_par$nu,
                b_alpha = gg_par$balpha, nu_alpha = gg_par$nualpha)
  pi_v <- as.numeric(gg_par$probpat)
  pp_gaga <- gg_fit$pp

  # Step 2: gold standard via fine 1D quadrature using the SAME hyperparameters
  gold <- gold_pp(X, groups, patterns, hyper, pi_v, n_grid = n_grid)

  # Step 3: PIG-MCMC slice sampler with the SAME hyperparameters
  out <- pig_gaga_slice_mcmc(X, groups, patterns, hyper, pi_v,
                             n_iter = n_iter, n_burn = n_burn,
                             trunc = trunc, slice_w = slice_w,
                             init_em_iter = 10, verbose = verbose)
  list(hyper = hyper, pi = pi_v,
       pp_gold = gold$pp, pp_gaga = pp_gaga, pp_pig = out$pp,
       gg_fit = gg_fit, pig_alpha_mean = out$alpha_mean)
}

# ---------- metrics: methods vs gold and methods vs truth -----------------
fdr_call_pp <- function(pp, fdr = 0.05) {
  null_prob <- pp[, 1]
  ord <- order(null_prob)
  cum_fdr <- cumsum(null_prob[ord]) / seq_along(ord)
  cutoff <- max(c(0, which(cum_fdr <= fdr)))
  calls <- rep(FALSE, nrow(pp))
  if (cutoff > 0) calls[ord[seq_len(cutoff)]] <- TRUE
  calls
}

compute_auc <- function(score, truth) {
  ord <- order(score, decreasing = TRUE)
  tr <- truth[ord]
  P <- sum(tr); N <- sum(!tr)
  if (P == 0 || N == 0) return(NA_real_)
  tp <- cumsum(tr); fp <- cumsum(!tr)
  fpr <- c(0, fp / N); tpr <- c(0, tp / P)
  sum(diff(fpr) * (head(tpr, -1) + tail(tpr, -1)) / 2)
}

compute_metrics <- function(res, truth, fdr = 0.05) {
  pp_gold <- res$pp_gold; pp_gaga <- res$pp_gaga; pp_pig <- res$pp_pig
  s_gold <- 1 - pp_gold[, 1]; s_gaga <- 1 - pp_gaga[, 1]; s_pig <- 1 - pp_pig[, 1]
  c_gold <- fdr_call_pp(pp_gold, fdr)
  c_gaga <- fdr_call_pp(pp_gaga, fdr)
  c_pig  <- fdr_call_pp(pp_pig, fdr)

  per_method <- function(calls, score, label) {
    data.frame(
      method = label,
      n_de = sum(calls),
      tp = sum(calls & truth),
      fp = sum(calls & !truth),
      empirical_fdr = if (sum(calls) > 0) sum(calls & !truth) / sum(calls) else 0,
      power = sum(calls & truth) / max(sum(truth), 1),
      auc = compute_auc(score, truth),
      stringsAsFactors = FALSE)
  }
  per_method_tbl <- rbind(
    per_method(c_gold, s_gold, "gold"),
    per_method(c_gaga, s_gaga, "gaga"),
    per_method(c_pig,  s_pig,  "pig"))

  vs <- function(pp_a, pp_b, s_a, s_b, c_a, c_b, lab) {
    data.frame(
      pair = lab,
      mean_abs_pp_diff = mean(abs(pp_a[, 1] - pp_b[, 1]), na.rm = TRUE),
      max_abs_pp_diff  = max(abs(pp_a[, 1] - pp_b[, 1]), na.rm = TRUE),
      spearman_score   = stats::cor(s_a, s_b, method = "spearman"),
      jaccard_calls    = sum(c_a & c_b) / max(sum(c_a | c_b), 1),
      stringsAsFactors = FALSE)
  }
  cross <- rbind(
    vs(pp_gaga, pp_gold, s_gaga, s_gold, c_gaga, c_gold, "gaga_vs_gold"),
    vs(pp_pig,  pp_gold, s_pig,  s_gold, c_pig,  c_gold, "pig_vs_gold"),
    vs(pp_pig,  pp_gaga, s_pig,  s_gaga, c_pig,  c_gaga, "pig_vs_gaga"))

  list(per_method = per_method_tbl, cross = cross,
       hyper = res$hyper)
}
