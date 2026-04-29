args <- commandArgs(FALSE)
this_file <- sub("^--file=", "", args[grepl("^--file=", args)][1])
source(file.path(dirname(normalizePath(this_file, mustWork = TRUE)), "benchmark_lib.R"))

bench <- setup_benchmark()
require_pkgs(c("GIGrvg", "coda"))

set.seed(env_int("PIG_SEED", 13))
n_draws <- env_int("PIG_DRAWS", 3000)
burnin <- env_int("PIG_BURNIN", 1000)
trunc <- env_int("PIG_TRUNC", 300)

cases <- data.frame(
  delta = c(5, 10, 30),
  mu = c(7.19 / 6.05, 5.57 / 5.01, 5.09 / 4.26)
)

rows <- vector("list", nrow(cases))
for (i in seq_len(nrow(cases))) {
  delta <- cases$delta[i]
  mu <- cases$mu[i]
  exact <- special_exact_moments(delta, mu)
  ab <- calc_AB_miller(delta, mu)
  miller_draws <- stats::rgamma(n_draws, shape = ab["shape"], rate = ab["rate"])
  pig_time <- system.time({
    pig_draws <- sample_pig_shape_special(delta, mu, n_draws, burnin, trunc)
  })[["elapsed"]]
  rows[[i]] <- data.frame(
    delta = delta,
    mu = mu,
    exact_mean = exact$mean,
    exact_var = exact$var,
    miller_mean = ab["shape"] / ab["rate"],
    miller_var = ab["shape"] / ab["rate"]^2,
    miller_sample_mean = mean(miller_draws),
    miller_sample_var = stats::var(miller_draws),
    pig_mean = mean(pig_draws),
    pig_var = stats::var(pig_draws),
    pig_ess = as.numeric(coda::effectiveSize(pig_draws)),
    pig_ess_per_draw = as.numeric(coda::effectiveSize(pig_draws)) / length(pig_draws),
    pig_time = pig_time,
    pig_mean_rel_error = (mean(pig_draws) - exact$mean) / exact$mean,
    pig_var_rel_error = (stats::var(pig_draws) - exact$var) / exact$var,
    miller_mean_rel_error = (ab["shape"] / ab["rate"] - exact$mean) / exact$mean,
    miller_var_rel_error = (ab["shape"] / ab["rate"]^2 - exact$var) / exact$var
  )
}

out <- do.call(rbind, rows)
write.csv(out, result_file("special_pig_vs_miller.csv"), row.names = FALSE)
print(out)
