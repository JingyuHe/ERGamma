args <- commandArgs(FALSE)
this_file <- sub("^--file=", "", args[grepl("^--file=", args)][1])
this_file <- normalizePath(this_file, mustWork = TRUE)

bench <- normalizePath(file.path(dirname(this_file), ".."), mustWork = TRUE)
raw_dir <- file.path(bench, "data", "armstrong_broad_raw")
cel_dir <- file.path(raw_dir, "CEL")
out_dir <- file.path(bench, "results")
dir.create(cel_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

gcrma_lib <- file.path(bench, "r-lib-gcrma")
legacy_lib <- file.path(bench, "r-lib")
.libPaths(unique(c(gcrma_lib, .libPaths(), legacy_lib)))

require_pkgs <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      "Missing packages: ", paste(missing, collapse = ", "),
      ". Run Rscript code2026/gaga_benchmark/R/install_gcrma_deps.R first.",
      call. = FALSE
    )
  }
}

read_sample_key <- function(path) {
  key <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  names(key) <- c("sample_name", "scan_name", "linear_scaling_factor", "paper_figure")
  key$group <- sub("_.*$", "", key$sample_name)
  key
}

extract_archives <- function(raw_dir, cel_dir) {
  archives <- file.path(
    raw_dir,
    c(
      "mll_scans_ALL1.tar.gz",
      "mll_scans_ALL2.tar.gz",
      "mll_scans_MLL1.tar.gz",
      "mll_scans_MLL2.tar.gz",
      "mll_scans_AML1.tar.gz",
      "mll_scans_AML2.tar.gz"
    )
  )
  missing <- archives[!file.exists(archives)]
  if (length(missing) > 0) {
    stop(
      "Raw CEL archives are missing:\n",
      paste(basename(missing), collapse = "\n"),
      "\nPlace the Broad archives in ", raw_dir,
      " and rerun. The archived Broad page lists these files, but Wayback does not contain the tar payloads.",
      call. = FALSE
    )
  }
  for (archive in archives) {
    utils::untar(archive, exdir = cel_dir)
  }
}

cel_files <- function(cel_dir) {
  list.files(cel_dir, pattern = "\\.CEL(\\.gz)?$", recursive = TRUE,
             full.names = TRUE, ignore.case = TRUE)
}

scan_from_cel <- function(path) {
  x <- basename(path)
  x <- sub("\\.gz$", "", x, ignore.case = TRUE)
  sub("\\.CEL$", "", x, ignore.case = TRUE)
}

read_cdf_name <- function(file) {
  hdr <- affyio::read.celfile.header(file)
  cdf <- hdr$cdfName
  if (is.null(cdf)) {
    cdf <- hdr$chiptype
  }
  as.character(cdf)[1]
}

metadata_for_cels <- function(files, key) {
  meta <- data.frame(
    cel_file = files,
    scan_name = scan_from_cel(files),
    stringsAsFactors = FALSE
  )
  meta$cdf_name <- vapply(files, read_cdf_name, character(1))
  merge(meta, key, by = "scan_name", all.x = TRUE, sort = FALSE)
}

run_gcrma <- function(files) {
  gcrma::just.gcrma(filenames = files)
}

scale_gcrma_if_requested <- function(expr, meta) {
  if (Sys.getenv("ARMSTRONG_APPLY_BROAD_SCALING", unset = "0") == "1") {
    factors <- meta$linear_scaling_factor
    names(factors) <- meta$sample_name
    expr <- sweep(expr, 2, factors[colnames(expr)], "*")
  }
  expr
}

require_pkgs(c("affy", "affyio", "gcrma"))

key_file <- file.path(raw_dir, "scaling_factors_and_fig_key.txt")
if (!file.exists(key_file)) {
  stop("Missing sample key: ", key_file, call. = FALSE)
}
sample_key <- read_sample_key(key_file)

if (length(cel_files(cel_dir)) == 0) {
  extract_archives(raw_dir, cel_dir)
}

all_cels <- cel_files(cel_dir)
if (length(all_cels) == 0) {
  stop("No CEL files found in ", cel_dir, call. = FALSE)
}

meta <- metadata_for_cels(all_cels, sample_key)
write.csv(meta, file.path(out_dir, "armstrong_gcrma_cel_metadata_all.csv"), row.names = FALSE)

keep <- meta$group %in% c("ALL", "MLL")
if (Sys.getenv("ARMSTRONG_KEEP_U95AV2", unset = "0") != "1") {
  keep <- keep & !grepl("U95AV2", meta$cdf_name, ignore.case = TRUE)
}
meta_keep <- meta[keep, , drop = FALSE]

if (length(unique(meta_keep$cdf_name)) != 1) {
  stop(
    "Selected CEL files contain multiple CDFs: ",
    paste(unique(meta_keep$cdf_name), collapse = ", "),
    ". Inspect armstrong_gcrma_cel_metadata_all.csv and adjust filtering.",
    call. = FALSE
  )
}

eset <- run_gcrma(meta_keep$cel_file)
expr <- Biobase::exprs(eset)
colnames(expr) <- scan_from_cel(colnames(expr))
idx <- match(meta_keep$scan_name, colnames(expr))
expr <- expr[, idx, drop = FALSE]
colnames(expr) <- meta_keep$sample_name
expr <- scale_gcrma_if_requested(expr, meta_keep)

groups <- factor(meta_keep$group, levels = c("ALL", "MLL"))
names(groups) <- meta_keep$sample_name

saveRDS(
  list(expr = expr, groups = groups, sample_info = meta_keep),
  file.path(out_dir, "armstrong_gcrma_expression.rds")
)
write.csv(expr, file.path(out_dir, "armstrong_gcrma_expression.csv"))
write.csv(meta_keep, file.path(out_dir, "armstrong_gcrma_cel_metadata_used.csv"), row.names = FALSE)

summary <- data.frame(
  n_probes = nrow(expr),
  n_samples = ncol(expr),
  n_all = sum(groups == "ALL"),
  n_mll = sum(groups == "MLL"),
  cdf_name = unique(meta_keep$cdf_name),
  min_expr = min(expr, na.rm = TRUE),
  max_expr = max(expr, na.rm = TRUE),
  applied_broad_scaling = Sys.getenv("ARMSTRONG_APPLY_BROAD_SCALING", unset = "0") == "1",
  stringsAsFactors = FALSE
)
write.csv(summary, file.path(out_dir, "armstrong_gcrma_summary.csv"), row.names = FALSE)
print(summary)
