#!/usr/bin/env Rscript
# Stage 1: load MAQC GEO matrices, parse annotations, cache to RDS so the
# downstream fitting steps don't re-read 30+ MB of compressed text each call.
args <- commandArgs(FALSE)
this_file <- sub("^--file=", "", args[grepl("^--file=", args)][1])
src_dir <- dirname(normalizePath(this_file, mustWork = TRUE))
source(file.path(src_dir, "benchmark_lib.R"))
bench <- setup_benchmark()
data_dir <- file.path(bench, "data")
cache_file <- result_file("maqc_loaded.rds")
if (file.exists(cache_file)) {
  cat("Cache already exists at:", cache_file, "\n"); quit(status = 0)
}

parse_meta_line <- function(line) {
  x <- scan(text = line, what = character(), sep = "\t", quote = "\"",
            quiet = TRUE, comment.char = ""); x[-1]
}
read_series_metadata <- function(file) {
  con <- gzfile(file, "rt"); on.exit(close(con), add = TRUE)
  out <- list()
  repeat {
    line <- readLines(con, n = 1)
    if (!length(line) || identical(line, "!series_matrix_table_begin")) break
    if (startsWith(line, "!Sample_title"))
      out$title <- parse_meta_line(line)
    else if (startsWith(line, "!Sample_geo_accession"))
      out$geo_accession <- parse_meta_line(line)
  }
  data.frame(gsm = out$geo_accession, title = out$title,
             stringsAsFactors = FALSE)
}
read_series_matrix <- function(file) {
  con <- gzfile(file, "rt"); on.exit(close(con), add = TRUE)
  dat <- read.delim(con, comment.char = "!", check.names = FALSE,
                    stringsAsFactors = FALSE)
  ids <- dat[[1]]
  x <- as.matrix(data.frame(lapply(dat[-1], as.numeric), check.names = FALSE))
  rownames(x) <- ids; x
}
read_platform_table <- function(file, gz = TRUE) {
  con <- if (gz) gzfile(file, "rt") else file(file, "rt")
  on.exit(close(con), add = TRUE)
  lines <- readLines(con, warn = FALSE)
  begin <- match("!platform_table_begin", lines)
  end <- match("!platform_table_end", lines)
  read.delim(text = paste(lines[(begin + 1):(end - 1)], collapse = "\n"),
             check.names = FALSE, stringsAsFactors = FALSE)
}

cat("Loading MAQC affy series matrix...\n")
affy_meta <- read_series_metadata(
  file.path(data_dir, "GSE5350-GPL570_series_matrix.txt.gz"))
affy <- read_series_matrix(
  file.path(data_dir, "GSE5350-GPL570_series_matrix.txt.gz"))
cat("Loading qPCR series matrix...\n")
qpcr_meta <- read_series_metadata(
  file.path(data_dir, "GSE5350-GPL4097_TaqMan_series_matrix.txt.gz"))
qpcr <- read_series_matrix(
  file.path(data_dir, "GSE5350-GPL4097_TaqMan_series_matrix.txt.gz"))
cat("Loading GPL570 annotation...\n")
gpl570 <- read_platform_table(file.path(data_dir, "GPL570.annot.gz"), gz = TRUE)
cat("Loading GPL4097 annotation...\n")
gpl4097 <- read_platform_table(
  file.path(data_dir, "GPL4097_family.soft.gz"), gz = TRUE)

saveRDS(list(affy_meta = affy_meta, affy = affy,
             qpcr_meta = qpcr_meta, qpcr = qpcr,
             gpl570 = gpl570, gpl4097 = gpl4097),
        cache_file)
cat("Saved cache:", cache_file, "\n")
