args <- commandArgs(trailingOnly = FALSE)
file_arg <- args[grepl("^--file=", args)]
if (length(file_arg) == 0) {
  stop("Run this script with Rscript so the benchmark directory can be found.")
}

this_file <- normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)
bench_dir <- normalizePath(file.path(dirname(this_file), ".."), mustWork = TRUE)
lib_dir <- file.path(bench_dir, "r-lib-gcrma")
dir.create(lib_dir, recursive = TRUE, showWarnings = FALSE)
.libPaths(unique(c(lib_dir, .libPaths())))

cran_repo <- "https://cloud.r-project.org"
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = cran_repo, lib = lib_dir)
}

BiocManager::install(
  c(
    "affy",
    "affyio",
    "gcrma",
    "hgu95acdf",
    "hgu95aprobe",
    "hgu95av2cdf",
    "hgu95av2probe"
  ),
  lib = lib_dir,
  ask = FALSE,
  update = FALSE
)

cat("Library paths:\n")
print(.libPaths())
cat("Installed package checks:\n")
for (pkg in c("affy", "affyio", "gcrma", "hgu95acdf", "hgu95aprobe", "hgu95av2cdf", "hgu95av2probe")) {
  cat(pkg, requireNamespace(pkg, quietly = TRUE), "\n")
}
