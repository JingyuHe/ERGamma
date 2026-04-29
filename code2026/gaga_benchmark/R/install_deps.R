args <- commandArgs(trailingOnly = FALSE)
file_arg <- args[grepl("^--file=", args)]
if (length(file_arg) == 0) {
  stop("Run this script with Rscript so the benchmark directory can be found.")
}

this_file <- normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)
bench_dir <- normalizePath(file.path(dirname(this_file), ".."), mustWork = TRUE)
lib_dir <- file.path(bench_dir, "r-lib")
legacy_lib <- normalizePath(file.path(bench_dir, "..", "r-lib"), mustWork = FALSE)
dir.create(lib_dir, recursive = TRUE, showWarnings = FALSE)
.libPaths(unique(c(lib_dir, legacy_lib, .libPaths())))

cran_repo <- "https://cloud.r-project.org"
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = cran_repo, lib = lib_dir)
}

BiocManager::install(
  c("gaga", "limma"),
  lib = lib_dir,
  ask = FALSE,
  update = FALSE
)

cat("Library paths:\n")
print(.libPaths())
cat("Installed package checks:\n")
for (pkg in c("gaga", "limma", "Biobase", "EBarrays", "coda", "GIGrvg")) {
  cat(pkg, requireNamespace(pkg, quietly = TRUE), "\n")
}
