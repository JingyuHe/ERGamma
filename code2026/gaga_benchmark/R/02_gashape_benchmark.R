args <- commandArgs(FALSE)
this_file <- sub("^--file=", "", args[grepl("^--file=", args)][1])
source(file.path(dirname(normalizePath(this_file, mustWork = TRUE)), "benchmark_lib.R"))

bench <- setup_benchmark()
require_pkgs(c("gaga", "Biobase"))

set.seed(11)
n_genes <- env_int("GAGA_GAS_N", 300)
n_cases <- env_int("GAGA_GAS_CASES", 12)
n_draws <- env_int("GAGA_GAS_DRAWS", 5000)

a0 <- 25.5
nu <- 0.109
balpha <- 1.183
nualpha <- 1683
xsim <- gaga::simGG(n_genes, m = c(10, 10), p.de = 0.05,
                    a0, nu, balpha, nualpha, equalcv = TRUE)
x <- Biobase::exprs(xsim)
groups <- Biobase::pData(xsim)$group

param_rows <- list()
row_id <- 1
for (gene in seq_len(nrow(x))) {
  for (grp in unique(groups)) {
    vals <- x[gene, groups == grp]
    a <- length(vals)
    b <- balpha
    c <- balpha / nualpha - sum(log(vals))
    d <- a0
    r <- a0 / nu
    s <- sum(vals)
    ar <- gas_approx_shape_rate(a, b, c, d, r, s)
    if (all(is.finite(ar)) && ar["shape"] > 0 && ar["rate"] > 0) {
      param_rows[[row_id]] <- data.frame(
        gene = gene, group = grp, a = a, b = b, c = c, d = d, r = r, s = s,
        approx_mean = ar["shape"] / ar["rate"]
      )
      row_id <- row_id + 1
    }
  }
}
params <- do.call(rbind, param_rows)
ord <- order(params$approx_mean)
pick <- unique(round(seq(1, nrow(params), length.out = min(n_cases, nrow(params)))))
params <- params[ord[pick], ]

rows <- vector("list", nrow(params))
for (i in seq_len(nrow(params))) {
  p <- params[i, ]
  exact_time <- system.time({
    exact <- gas_exact_moments(p$a, p$b, p$c, p$d, p$r, p$s)
  })[["elapsed"]]
  gaga_time <- system.time({
    gm <- gaga::mcgamma(p$a, p$b, p$c, p$d, p$r, p$s)
    draws <- gaga::rcgamma(n_draws, p$a, p$b, p$c, p$d, p$r, p$s)
  })[["elapsed"]]
  rows[[i]] <- data.frame(
    gene = p$gene,
    group = p$group,
    a = p$a,
    b = p$b,
    c = p$c,
    d = p$d,
    r = p$r,
    s = p$s,
    exact_mean = exact$mean,
    exact_var = exact$var,
    exact_mode = exact$mode,
    gaga_mcgamma_mean = gm$m,
    gaga_mcgamma_var = gm$v,
    gaga_sample_mean = mean(draws),
    gaga_sample_var = stats::var(draws),
    mean_rel_error = (gm$m - exact$mean) / exact$mean,
    var_rel_error = (gm$v - exact$var) / exact$var,
    exact_time = exact_time,
    gaga_time = gaga_time
  )
}

out <- do.call(rbind, rows)
write.csv(out, result_file("gas_approx_benchmark.csv"), row.names = FALSE)
print(out)
