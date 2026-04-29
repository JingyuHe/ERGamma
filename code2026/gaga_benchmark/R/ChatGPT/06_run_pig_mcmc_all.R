args <- commandArgs(FALSE)
this_file <- sub("^--file=", "", args[grepl("^--file=", args)][1])
src_dir <- dirname(normalizePath(this_file, mustWork = TRUE))

run_script <- function(name) {
  cat("\n=== Running ", name, " ===\n", sep = "")
  status <- system2(file.path(R.home("bin"), "Rscript"),
                    file.path(src_dir, name))
  if (!identical(status, 0L)) {
    stop(name, " failed with status ", status)
  }
}

run_script("04_armstrong_pig_mcmc_gaga.R")
run_script("05_maqc_pig_mcmc_gaga.R")
