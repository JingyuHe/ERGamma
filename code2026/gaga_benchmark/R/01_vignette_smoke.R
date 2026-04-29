args <- commandArgs(FALSE)
this_file <- sub("^--file=", "", args[grepl("^--file=", args)][1])
source(file.path(dirname(normalizePath(this_file, mustWork = TRUE)), "benchmark_lib.R"))

bench <- setup_benchmark()
require_pkgs(c("gaga", "Biobase"))

set.seed(10)
n <- env_int("GAGA_SMOKE_N", 100)
m <- c(6, 6)
a0 <- 25.5
nu <- 0.109
balpha <- 1.183
nualpha <- 1683
probpat <- c(0.95, 0.05)

xsim <- gaga::simGG(n, m = m, p.de = probpat[2], a0, nu, balpha, nualpha, equalcv = TRUE)
x <- Biobase::exprs(xsim)
groups_full <- Biobase::pData(xsim)$group
train_cols <- setdiff(seq_len(ncol(x)), c(6, 12))
groups <- groups_full[train_cols]
patterns <- matrix(c(0, 0, 0, 1), 2, 2)
colnames(patterns) <- c("group 1", "group 2")
truth <- abs(Biobase::fData(xsim)$mean.1 - Biobase::fData(xsim)$mean.2) > 1e-12

run_gaga <- function(label, nclust) {
  elapsed <- system.time({
    fit <- gaga::fitGG(
      x[, train_cols],
      groups,
      patterns = patterns,
      nclust = nclust,
      method = "quickEM",
      trace = FALSE
    )
    fit <- gaga::parest(fit, x = x[, train_cols], groups = groups, alpha = 0.05)
    genes <- gaga::findgenes(fit, x[, train_cols], groups, fdrmax = 0.05, parametric = TRUE)
  })[["elapsed"]]
  cbind(method = label, eval_calls(genes$d != 0, truth), elapsed = elapsed)
}

rows <- list(
  run_gaga("GaGa_quickEM", 1),
  run_gaga("MiGaGa2_quickEM", 2)
)

lim <- limma_or_ttest(x[, train_cols], groups, fdr = 0.05)
rows[[length(rows) + 1]] <- cbind(method = "limma_or_ttest_BH",
                                  eval_calls(lim$calls, truth),
                                  elapsed = NA_real_)

out <- do.call(rbind, rows)
write.csv(out, result_file("vignette_smoke_summary.csv"), row.names = FALSE)
print(out)
