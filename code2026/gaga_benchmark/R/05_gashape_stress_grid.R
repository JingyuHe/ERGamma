args <- commandArgs(FALSE)
this_file <- sub("^--file=", "", args[grepl("^--file=", args)][1])
source(file.path(dirname(normalizePath(this_file, mustWork = TRUE)), "benchmark_lib.R"))

bench <- setup_benchmark()
require_pkgs("gaga")

set.seed(env_int("GAGA_STRESS_SEED", 14))
n_cases <- env_int("GAGA_STRESS_CASES", 120)

draw_params <- function(n) {
  a <- sample(c(1, 2, 3, 5, 10, 20, 50), n, replace = TRUE)
  b <- exp(stats::runif(n, log(0.15), log(12)))
  d <- exp(stats::runif(n, log(0.15), log(50)))
  r <- exp(stats::runif(n, log(0.01), log(500)))
  s <- exp(stats::runif(n, log(0.05), log(100)))
  approx_rate <- exp(stats::runif(n, log(0.03), log(10)))
  c <- approx_rate - a * log(s / a)
  data.frame(a = a, b = b, c = c, d = d, r = r, s = s,
             approx_rate = approx_rate)
}

params <- draw_params(n_cases)
rows <- vector("list", nrow(params))
for (i in seq_len(nrow(params))) {
  p <- params[i, ]
  res <- tryCatch({
    exact <- gas_exact_moments(p$a, p$b, p$c, p$d, p$r, p$s, rel.tol = 1e-7)
    gm <- gaga::mcgamma(p$a, p$b, p$c, p$d, p$r, p$s)
    data.frame(
      case = i,
      a = p$a,
      b = p$b,
      c = p$c,
      d = p$d,
      r = p$r,
      s = p$s,
      exact_mean = exact$mean,
      exact_var = exact$var,
      exact_mode = exact$mode,
      gaga_mean = gm$m,
      gaga_var = gm$v,
      mean_rel_error = (gm$m - exact$mean) / exact$mean,
      var_rel_error = (gm$v - exact$var) / exact$var,
      abs_mean_rel_error = abs((gm$m - exact$mean) / exact$mean),
      abs_var_rel_error = abs((gm$v - exact$var) / exact$var),
      ok = TRUE,
      error = NA_character_
    )
  }, error = function(e) {
    data.frame(
      case = i,
      a = p$a,
      b = p$b,
      c = p$c,
      d = p$d,
      r = p$r,
      s = p$s,
      exact_mean = NA_real_,
      exact_var = NA_real_,
      exact_mode = NA_real_,
      gaga_mean = NA_real_,
      gaga_var = NA_real_,
      mean_rel_error = NA_real_,
      var_rel_error = NA_real_,
      abs_mean_rel_error = NA_real_,
      abs_var_rel_error = NA_real_,
      ok = FALSE,
      error = conditionMessage(e)
    )
  })
  rows[[i]] <- res
}

out <- do.call(rbind, rows)
out <- out[order(-out$abs_mean_rel_error, -out$abs_var_rel_error), ]
write.csv(out, result_file("gas_stress_grid.csv"), row.names = FALSE)
write.csv(utils::head(out, 25), result_file("gas_stress_worst.csv"), row.names = FALSE)
print(utils::head(out, 10))
