#!/usr/bin/env Rscript
# Run a batch of (scenario, rep) pairs in one R process to avoid startup
# overhead. Skips already-done reps. Stops gracefully when bash kills it; the
# next call resumes since results are saved per-rep.
#
# Env vars:
#   SIM_TODO  = "A:1,A:2,A:3" or "all" (defaults to remaining of full design)
#   SIM_DESIGN_REPS  = 10 (default)

src_dir <- normalizePath(dirname(sub("^--file=", "",
  commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][1])))
source(file.path(src_dir, "sim_lib.R"))
results_dir <- file.path(src_dir, "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

design_reps <- as.integer(Sys.getenv("SIM_DESIGN_REPS", unset = "10"))
todo_env <- Sys.getenv("SIM_TODO", unset = "all")

if (todo_env == "all") {
  todo <- expand.grid(scenario = names(scenarios),
                      rep = seq_len(design_reps),
                      stringsAsFactors = FALSE)
} else {
  pairs <- strsplit(todo_env, ",")[[1]]
  todo <- do.call(rbind, lapply(pairs, function(p) {
    pp <- strsplit(p, ":")[[1]]
    data.frame(scenario = pp[1], rep = as.integer(pp[2]),
               stringsAsFactors = FALSE)
  }))
}

base_seed <- as.integer(Sys.getenv("SIM_BASE_SEED", unset = "100"))
n_iter <- as.integer(Sys.getenv("SIM_N_ITER", unset = "200"))
n_burn <- as.integer(Sys.getenv("SIM_N_BURN", unset = "60"))
trunc  <- as.integer(Sys.getenv("SIM_TRUNC", unset = "30"))
n_grid <- as.integer(Sys.getenv("SIM_GRID", unset = "600"))

t_start <- Sys.time()
budget_sec <- as.numeric(Sys.getenv("SIM_BUDGET_SEC", unset = "35"))

n_done_call <- 0L
for (k in seq_len(nrow(todo))) {
  if (as.numeric(Sys.time() - t_start, units = "secs") > budget_sec) {
    cat("[budget] stopping early; resume by re-running\n")
    break
  }
  s_id <- todo$scenario[k]; rep <- todo$rep[k]
  out_file <- file.path(results_dir,
                        sprintf("scenario_%s_rep_%02d.rds", s_id, rep))
  if (file.exists(out_file)) next
  spec <- scenarios[[s_id]]
  cat(sprintf("[%s rep %d] G=%d n=%d nu_alpha=%g ...\n",
              s_id, rep, spec$G, spec$n, spec$nu_alpha))
  dat <- simulate_data(spec, seed = base_seed +
                       1000 * (match(s_id, names(scenarios)) - 1) + rep)
  patterns <- matrix(c(0, 0, 0, 1), 2, 2)
  colnames(patterns) <- levels(dat$groups)
  t0 <- Sys.time()
  res <- run_methods(dat$X, dat$groups, patterns,
                     n_iter = n_iter, n_burn = n_burn,
                     trunc = trunc, n_grid = n_grid)
  metrics <- compute_metrics(res, dat$truth)
  el <- as.numeric(Sys.time() - t0, units = "secs")
  saveRDS(list(scenario = s_id, rep = rep, spec = spec,
               truth = dat$truth,
               hyper_oracle = list(alpha0 = spec$a0, nu = spec$nu,
                                   b_alpha = spec$b_alpha,
                                   nu_alpha = spec$nu_alpha,
                                   p_de = spec$p_de),
               hyper_estimated = res$hyper, pi_v = res$pi,
               pp_gold = res$pp_gold, pp_gaga = res$pp_gaga,
               pp_pig = res$pp_pig,
               pig_alpha_mean = res$pig_alpha_mean,
               metrics = metrics, elapsed_sec = el),
          out_file)
  n_done_call <- n_done_call + 1L
  cat(sprintf("  [done] %.1fs  gaga|Δpp|=%.4g  pig|Δpp|=%.4g\n",
              el,
              metrics$cross$mean_abs_pp_diff[metrics$cross$pair == "gaga_vs_gold"],
              metrics$cross$mean_abs_pp_diff[metrics$cross$pair == "pig_vs_gold"]))
}
cat(sprintf("[batch] this call completed %d reps. Wall=%.1fs\n",
            n_done_call,
            as.numeric(Sys.time() - t_start, units = "secs")))
