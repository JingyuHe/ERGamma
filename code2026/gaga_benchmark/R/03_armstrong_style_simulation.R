args <- commandArgs(FALSE)
this_file <- sub("^--file=", "", args[grepl("^--file=", args)][1])
source(file.path(dirname(normalizePath(this_file, mustWork = TRUE)), "benchmark_lib.R"))

bench <- setup_benchmark()
require_pkgs(c("gaga", "Biobase"))

set.seed(env_int("GAGA_SEED", 12))
n_reps <- env_int("GAGA_REPS", 5)
n_genes <- env_int("GAGA_N_GENES", 500)
sample_sizes <- c(5, 10, 15, 20)
max_m <- max(sample_sizes)

a0 <- 25.5
nu <- 0.109
balpha <- 1.183
nualpha <- 1683
probpat <- c(0.95, 0.05)
patterns <- matrix(c(0, 0, 0, 1), 2, 2)
colnames(patterns) <- c("group 1", "group 2")

run_gaga_calls <- function(x, groups, nclust) {
  fit <- gaga::fitGG(
    x,
    groups,
    patterns = patterns,
    nclust = nclust,
    method = "quickEM",
    trace = FALSE
  )
  fit <- gaga::parest(fit, x = x, groups = groups, alpha = 0.05)
  genes <- gaga::findgenes(fit, x, groups, fdrmax = 0.05, parametric = TRUE)
  genes$d != 0
}

rows <- list()
idx <- 1
for (rep_id in seq_len(n_reps)) {
  xsim <- gaga::simGG(n_genes, m = c(max_m, max_m), p.de = probpat[2],
                      a0, nu, balpha, nualpha, equalcv = TRUE)
  x_all <- Biobase::exprs(xsim)
  groups_all <- Biobase::pData(xsim)$group
  truth <- abs(Biobase::fData(xsim)$mean.1 - Biobase::fData(xsim)$mean.2) > 1e-12
  for (ss in sample_sizes) {
    cols <- c(which(groups_all == "group 1")[seq_len(ss)],
              which(groups_all == "group 2")[seq_len(ss)])
    x <- x_all[, cols]
    groups <- groups_all[cols]
    elapsed <- system.time({
      calls_gaga <- run_gaga_calls(x, groups, nclust = 1)
    })[["elapsed"]]
    rows[[idx]] <- cbind(rep = rep_id, n_per_group = ss, method = "GaGa_quickEM",
                         eval_calls(calls_gaga, truth), elapsed = elapsed)
    idx <- idx + 1

    elapsed <- system.time({
      calls_migaga <- run_gaga_calls(x, groups, nclust = 2)
    })[["elapsed"]]
    rows[[idx]] <- cbind(rep = rep_id, n_per_group = ss, method = "MiGaGa2_quickEM",
                         eval_calls(calls_migaga, truth), elapsed = elapsed)
    idx <- idx + 1

    elapsed <- system.time({
      lim <- limma_or_ttest(x, groups, fdr = 0.05)
    })[["elapsed"]]
    rows[[idx]] <- cbind(rep = rep_id, n_per_group = ss, method = "limma_or_ttest_BH",
                         eval_calls(lim$calls, truth), elapsed = elapsed)
    idx <- idx + 1
  }
}

out <- do.call(rbind, rows)
write.csv(out, result_file("armstrong_style_simulation.csv"), row.names = FALSE)

summary <- stats::aggregate(
  cbind(n_de, fdr, power, elapsed) ~ n_per_group + method,
  data = out,
  FUN = mean
)
write.csv(summary, result_file("armstrong_style_summary.csv"), row.names = FALSE)
print(summary)
