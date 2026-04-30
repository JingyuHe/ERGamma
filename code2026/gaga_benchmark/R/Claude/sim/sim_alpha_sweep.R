#!/usr/bin/env Rscript
# Sweep over nu_alpha (the per-gene shape prior mean) at fixed n, fixed
# remaining hyperparameters. Maps the curve "Stirling error vs alpha scale".

src_dir <- normalizePath(dirname(sub("^--file=", "",
  commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][1])))
source(file.path(src_dir, "sim_lib.R"))
results_dir <- file.path(src_dir, "results"); dir.create(results_dir, showWarnings = FALSE)

# fixed parameters; vary nu_alpha (== prior E[alpha])
sweep_spec_fn <- function(nu_alpha) {
  list(name = sprintf("alpha=%g", nu_alpha),
       G = 200, n = 5,
       b_alpha = max(2, nu_alpha * 0.5),  # weak/proportional concentration
       nu_alpha = nu_alpha,
       a0 = 5, nu = 5, p_de = 0.20)
}

n_reps_per_alpha <- as.integer(Sys.getenv("SWEEP_REPS", unset = "3"))
budget <- as.numeric(Sys.getenv("SIM_BUDGET_SEC", unset = "35"))
base_seed <- as.integer(Sys.getenv("SIM_BASE_SEED", unset = "5000"))
n_iter <- as.integer(Sys.getenv("SIM_N_ITER", unset = "200"))
n_burn <- as.integer(Sys.getenv("SIM_N_BURN", unset = "60"))
trunc  <- as.integer(Sys.getenv("SIM_TRUNC", unset = "30"))

todo <- expand.grid(nu_alpha = alpha_sweep_grid,
                    rep = seq_len(n_reps_per_alpha))
todo <- todo[order(todo$rep, todo$nu_alpha), ]

t_start <- Sys.time()
for (k in seq_len(nrow(todo))) {
  if (as.numeric(Sys.time() - t_start, units = "secs") > budget) {
    cat("[budget] stopping; resume by re-running.\n"); break
  }
  nu_a <- todo$nu_alpha[k]; rep <- todo$rep[k]
  out_file <- file.path(results_dir,
                        sprintf("sweep_nua_%05g_rep_%02d.rds", nu_a * 100, rep))
  if (file.exists(out_file)) next
  spec <- sweep_spec_fn(nu_a)
  cat(sprintf("[sweep] nu_alpha=%g rep=%d ...\n", nu_a, rep))
  dat <- simulate_data(spec, seed = base_seed + 100 * rep + round(nu_a * 10))
  patterns <- matrix(c(0, 0, 0, 1), 2, 2)
  colnames(patterns) <- levels(dat$groups)
  t0 <- Sys.time()
  res <- run_methods(dat$X, dat$groups, patterns,
                     n_iter = n_iter, n_burn = n_burn, trunc = trunc)
  metrics <- compute_metrics(res, dat$truth)
  el <- as.numeric(Sys.time() - t0, units = "secs")
  saveRDS(list(nu_alpha = nu_a, rep = rep, spec = spec,
               truth = dat$truth, hyper_oracle = list(
                 alpha0 = spec$a0, nu = spec$nu,
                 b_alpha = spec$b_alpha, nu_alpha = spec$nu_alpha,
                 p_de = spec$p_de),
               hyper_estimated = res$hyper, pi_v = res$pi,
               pp_gold = res$pp_gold, pp_gaga = res$pp_gaga,
               pp_pig = res$pp_pig,
               pig_alpha_mean = res$pig_alpha_mean,
               metrics = metrics, elapsed_sec = el),
          out_file)
  cat(sprintf("  done %.1fs  gaga|Δpp|=%.4g  pig|Δpp|=%.4g\n",
              el,
              metrics$cross$mean_abs_pp_diff[metrics$cross$pair == "gaga_vs_gold"],
              metrics$cross$mean_abs_pp_diff[metrics$cross$pair == "pig_vs_gold"]))
}
