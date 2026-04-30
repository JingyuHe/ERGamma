#!/usr/bin/env Rscript
# Run one (scenario, replicate) of the simulation and save to RDS.
#
# Env vars:
#   SIM_SCENARIO   = A | B | C | D
#   SIM_REP        = 1, 2, ...
#   SIM_BASE_SEED  = base seed (default 100)
#   SIM_N_ITER, SIM_N_BURN, SIM_TRUNC, SIM_GRID
#
# Output: sim/results/scenario_<id>_rep_<r>.rds

src_dir <- normalizePath(dirname(sub("^--file=", "",
  commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][1])))
source(file.path(src_dir, "sim_lib.R"))

results_dir <- file.path(src_dir, "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

s_id <- Sys.getenv("SIM_SCENARIO", unset = "A")
rep <- as.integer(Sys.getenv("SIM_REP", unset = "1"))
base_seed <- as.integer(Sys.getenv("SIM_BASE_SEED", unset = "100"))
n_iter <- as.integer(Sys.getenv("SIM_N_ITER", unset = "200"))
n_burn <- as.integer(Sys.getenv("SIM_N_BURN", unset = "60"))
trunc  <- as.integer(Sys.getenv("SIM_TRUNC", unset = "30"))
n_grid <- as.integer(Sys.getenv("SIM_GRID", unset = "600"))

stopifnot(s_id %in% names(scenarios))
spec <- scenarios[[s_id]]
out_file <- file.path(results_dir, sprintf("scenario_%s_rep_%02d.rds", s_id, rep))
if (file.exists(out_file)) {
  cat(sprintf("[skip] %s already exists\n", out_file)); quit(status = 0)
}

cat(sprintf("[sim] scenario %s rep %d: G=%d n=%d nu_alpha=%g a0=%g nu=%g\n",
            s_id, rep, spec$G, spec$n, spec$nu_alpha, spec$a0, spec$nu))

dat <- simulate_data(spec, seed = base_seed + 1000 * (match(s_id, names(scenarios)) - 1) + rep)
patterns <- matrix(c(0, 0, 0, 1), 2, 2)
colnames(patterns) <- levels(dat$groups)

t0 <- Sys.time()
res <- run_methods(dat$X, dat$groups, patterns,
                   n_iter = n_iter, n_burn = n_burn,
                   trunc = trunc, n_grid = n_grid, verbose = FALSE)
metrics <- compute_metrics(res, dat$truth)
elapsed <- as.numeric(Sys.time() - t0, units = "secs")

cat(sprintf("[done] %.1fs  hyper=(a0=%.2f, nu=%.2f, b_a=%.2f, nu_a=%.2f)\n",
            elapsed, res$hyper$alpha0, res$hyper$nu,
            res$hyper$b_alpha, res$hyper$nu_alpha))
cat("[methods]\n"); print(metrics$per_method, row.names = FALSE)
cat("[cross]\n");   print(metrics$cross, row.names = FALSE)

saveRDS(list(
  scenario = s_id, rep = rep, spec = spec,
  truth = dat$truth, hyper_oracle = list(
    alpha0 = spec$a0, nu = spec$nu,
    b_alpha = spec$b_alpha, nu_alpha = spec$nu_alpha,
    p_de = spec$p_de),
  hyper_estimated = res$hyper,
  pi_v = res$pi,
  pp_gold = res$pp_gold, pp_gaga = res$pp_gaga, pp_pig = res$pp_pig,
  pig_alpha_mean = res$pig_alpha_mean,
  metrics = metrics, elapsed_sec = elapsed),
  out_file)
cat(sprintf("[saved] %s\n", out_file))
