args <- commandArgs(FALSE)
this_file <- sub("^--file=", "", args[grepl("^--file=", args)][1])
root <- dirname(normalizePath(this_file, mustWork = TRUE))

run_script <- function(name) {
  path <- file.path(root, name)
  cat("\n== ", name, " ==\n", sep = "")
  status <- system2(file.path(R.home("bin"), "Rscript"), path)
  if (!identical(status, 0L)) {
    stop(name, " failed with status ", status, call. = FALSE)
  }
}

run_script("01_armstrong_quad_gaga.R")
run_script("02_maqc_quad_gaga.R")
