args <- commandArgs(FALSE)
this_file <- sub("^--file=", "", args[grepl("^--file=", args)][1])
source(file.path(dirname(normalizePath(this_file, mustWork = TRUE)), "benchmark_lib.R"))

bench <- setup_benchmark()
paths <- .libPaths()
system_paths <- paths[grepl("^/opt/homebrew", paths)]
local_paths <- paths[!grepl("^/opt/homebrew", paths)]
.libPaths(unique(c(system_paths, file.path(bench, "r-lib-gcrma"), local_paths)))
require_pkgs(c("gaga", "limma", "EBarrays"))

data_dir <- file.path(bench, "data")

parse_meta_line <- function(line) {
  x <- scan(text = line, what = character(), sep = "\t", quote = "\"",
            quiet = TRUE, comment.char = "")
  x[-1]
}

read_series_metadata <- function(file) {
  con <- gzfile(file, "rt")
  on.exit(close(con), add = TRUE)
  out <- list()
  repeat {
    line <- readLines(con, n = 1)
    if (!length(line) || identical(line, "!series_matrix_table_begin")) {
      break
    }
    if (startsWith(line, "!Sample_title")) {
      out$title <- parse_meta_line(line)
    } else if (startsWith(line, "!Sample_geo_accession")) {
      out$geo_accession <- parse_meta_line(line)
    } else if (startsWith(line, "!Sample_source_name_ch1")) {
      out$source <- parse_meta_line(line)
    } else if (startsWith(line, "!Sample_platform_id")) {
      out$platform <- parse_meta_line(line)
    }
  }
  if (is.null(out$title) || is.null(out$geo_accession)) {
    stop("Could not parse GEO sample metadata from ", file)
  }
  data.frame(
    gsm = out$geo_accession,
    title = out$title,
    source = if (is.null(out$source)) NA_character_ else out$source,
    platform = if (is.null(out$platform)) NA_character_ else out$platform,
    stringsAsFactors = FALSE
  )
}

read_series_matrix <- function(file) {
  con <- gzfile(file, "rt")
  on.exit(close(con), add = TRUE)
  dat <- read.delim(
    con,
    comment.char = "!",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  ids <- dat[[1]]
  x <- as.matrix(data.frame(lapply(dat[-1], as.numeric), check.names = FALSE))
  rownames(x) <- ids
  x
}

read_platform_table <- function(file, gz = TRUE) {
  con <- if (gz) gzfile(file, "rt") else file(file, "rt")
  on.exit(close(con), add = TRUE)
  lines <- readLines(con, warn = FALSE)
  begin <- match("!platform_table_begin", lines)
  end <- match("!platform_table_end", lines)
  if (is.na(begin) || is.na(end)) {
    stop("Could not find platform table in ", file)
  }
  read.delim(
    text = paste(lines[(begin + 1):(end - 1)], collapse = "\n"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

split_symbols <- function(x) {
  x <- toupper(trimws(x))
  x <- strsplit(x, "///", fixed = TRUE)
  lapply(x, function(z) {
    z <- trimws(z)
    z[nzchar(z) & z != "---" & z != "NA"]
  })
}

make_affy_symbol_map <- function(gpl570) {
  sym_col <- "Gene symbol"
  if (!sym_col %in% names(gpl570)) {
    stop("GPL570 annotation missing Gene symbol column.")
  }
  syms <- split_symbols(gpl570[[sym_col]])
  rows <- lapply(seq_along(syms), function(i) {
    if (!length(syms[[i]])) {
      return(NULL)
    }
    data.frame(probe_id = gpl570$ID[i], symbol = syms[[i]], stringsAsFactors = FALSE)
  })
  unique(do.call(rbind, rows))
}

maqc_patterns <- function() {
  pat <- matrix(
    c(
      0, 0, 0, 0,
      0, 1, 1, 1,
      0, 0, 1, 1,
      0, 0, 0, 1,
      0, 1, 2, 3
    ),
    ncol = 4,
    byrow = TRUE
  )
  colnames(pat) <- c("A", "C", "D", "B")
  pat
}

maqc_eb_patterns <- function(groups) {
  pattern_values <- list(
    rep(1, length(groups)),
    ifelse(groups == "A", 1, 2),
    ifelse(groups %in% c("A", "C"), 1, 2),
    ifelse(groups == "B", 2, 1),
    match(groups, c("A", "C", "D", "B"))
  )
  EBarrays::ebPatterns(vapply(pattern_values, paste, collapse = " ", character(1)))
}

fit_gaga_maqc <- function(x, groups, nclust, method = "EM") {
  patterns <- maqc_patterns()
  fit_methods <- unique(c(method, "quickEM", "EM"))
  errors <- character()
  for (fit_method in fit_methods) {
    fit <- tryCatch({
      gaga::fitGG(
        x,
        groups,
        patterns = patterns,
        equalcv = TRUE,
        nclust = nclust,
        method = fit_method,
        trace = FALSE
      )
    }, error = function(e) {
      errors <<- c(errors, paste0("fitGG/", fit_method, ": ", conditionMessage(e)))
      NULL
    })
    if (is.null(fit)) {
      next
    }
    fit <- tryCatch({
      gaga::parest(fit, x = x, groups = groups, alpha = 0.05)
    }, error = function(e) {
      errors <<- c(errors, paste0("parest/", fit_method, ": ", conditionMessage(e)))
      NULL
    })
    if (!is.null(fit)) {
      fit$method_used <- fit_method
      return(fit)
    }
  }
  stop("GaGa/MiGaGa fit failed: ", paste(errors, collapse = " | "))
}

fit_ga_maqc <- function(x, groups) {
  patterns <- maqc_eb_patterns(groups)
  attempts <- list(
    raw = x,
    log2 = log2(pmax(x, .Machine$double.eps))
  )
  errors <- character()
  for (nm in names(attempts)) {
    xx <- attempts[[nm]]
    fit <- tryCatch({
      EBarrays::emfit(
        data = xx,
        family = "GG",
        hypotheses = patterns,
        num.iter = 20,
        verbose = FALSE
      )
    }, error = function(e) {
      errors <<- c(errors, paste0(nm, ": ", conditionMessage(e)))
      NULL
    })
    if (!is.null(fit)) {
      post <- EBarrays::postprob(fit, xx)$pattern
      return(list(fit = fit, pp = post, method_used = paste0("EBarrays_GG_", nm)))
    }
  }
  stop("Ga/Gamma-Gamma fit failed: ", paste(errors, collapse = " | "))
}

limma_f_score <- function(x, groups) {
  groups <- factor(groups)
  design <- stats::model.matrix(~ groups)
  fit <- limma::lmFit(log2(pmax(x, .Machine$double.eps)), design)
  fit <- limma::eBayes(fit)
  tt <- limma::topTable(fit, coef = seq_len(ncol(design))[-1], number = Inf,
                        sort.by = "none")
  p <- tt$P.Value
  names(p) <- rownames(tt)
  -log10(pmax(p, 1e-300))
}

score_by_symbol <- function(scores, affy_map) {
  idx <- match(affy_map$probe_id, names(scores))
  map <- affy_map[!is.na(idx), , drop = FALSE]
  map$score <- scores[idx[!is.na(idx)]]
  stats::aggregate(score ~ symbol, data = map, FUN = max, na.rm = TRUE)
}

limma_qpcr_validation <- function(qpcr, groups) {
  groups <- factor(groups)
  design <- stats::model.matrix(~ groups)
  fit <- limma::lmFit(log2(pmax(qpcr, .Machine$double.eps)), design)
  fit <- limma::eBayes(fit)
  tt <- limma::topTable(fit, coef = seq_len(ncol(design))[-1], number = Inf,
                        sort.by = "none")
  p <- tt$P.Value
  names(p) <- rownames(tt)
  p
}

collapse_probe_symbols <- function(affy_map) {
  stats::aggregate(symbol ~ probe_id, data = affy_map, FUN = function(z) paste(unique(z), collapse = ";"))
}

top_pattern4_lists <- function(fits, affy_fit, groups, probe_symbols, n = 1000) {
  means <- sapply(c("A", "B"), function(g) rowMeans(affy_fit[, groups == g, drop = FALSE]))
  rownames(means) <- rownames(affy_fit)
  rows <- list()
  row_id <- 1
  for (method in names(fits)) {
    pp4 <- fits[[method]]$pp[, 5]
    names(pp4) <- rownames(affy_fit)
    for (direction in c("B_gt_A", "A_gt_B")) {
      keep <- if (direction == "B_gt_A") means[, "B"] > means[, "A"] else means[, "A"] > means[, "B"]
      ord <- order(pp4[keep], decreasing = TRUE)
      ids <- names(pp4[keep])[ord][seq_len(min(n, sum(keep)))]
      out <- data.frame(
        method = method,
        direction = direction,
        rank = seq_along(ids),
        probe_id = ids,
        pattern4_prob = pp4[ids],
        mean_A = means[ids, "A"],
        mean_B = means[ids, "B"],
        log2_B_over_A = log2(pmax(means[ids, "B"], .Machine$double.eps) / pmax(means[ids, "A"], .Machine$double.eps)),
        stringsAsFactors = FALSE
      )
      out <- merge(out, probe_symbols, by = "probe_id", all.x = TRUE, sort = FALSE)
      out <- out[order(out$rank), ]
      rows[[row_id]] <- out
      row_id <- row_id + 1
    }
  }
  do.call(rbind, rows)
}

validation_curve <- function(score, validated) {
  ok <- is.finite(score) & !is.na(validated)
  score <- score[ok]
  validated <- as.logical(validated[ok])
  ord <- order(score, decreasing = TRUE)
  validated <- validated[ord]
  n <- length(validated)
  x <- c(0, cumsum(!validated) / n)
  y <- c(0, cumsum(validated) / n)
  auc <- sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
  std_auc <- if (sum(!validated) == 0 || sum(validated) == 0) {
    NA_real_
  } else {
    auc / ((sum(!validated) / n) * (sum(validated) / n))
  }
  data.frame(
    rank = seq_along(x) - 1,
    prop_nonvalidated = x,
    prop_validated = y,
    auc_scaled = auc,
    auc_standard = std_auc
  )
}

first_site_filter <- function(meta) {
  grepl("^MAQC_AFX_1_[ABCD][1-5]$", meta$title)
}

group_from_title <- function(title) {
  sub("^MAQC_[A-Z]+_[0-9]+_([ABCD]).*$", "\\1", title)
}

n_genes <- env_int("MAQC_N_GENES", 0)
gaga_method <- Sys.getenv("MAQC_GAGA_METHOD", unset = "EM")

affy_file <- file.path(data_dir, "GSE5350-GPL570_series_matrix.txt.gz")
qpcr_file <- file.path(data_dir, "GSE5350-GPL4097_TaqMan_series_matrix.txt.gz")
gpl570_file <- file.path(data_dir, "GPL570.annot.gz")
gpl4097_file <- file.path(data_dir, "GPL4097_family.soft.gz")

affy_meta <- read_series_metadata(affy_file)
affy <- read_series_matrix(affy_file)
qpcr_meta <- read_series_metadata(qpcr_file)
qpcr <- read_series_matrix(qpcr_file)
gpl570 <- read_platform_table(gpl570_file, gz = TRUE)
gpl4097 <- read_platform_table(gpl4097_file, gz = TRUE)

keep <- first_site_filter(affy_meta)
affy_meta_1 <- affy_meta[keep, , drop = FALSE]
affy_1 <- affy[, affy_meta_1$gsm, drop = FALSE]
colnames(affy_1) <- affy_meta_1$title
groups_affy <- group_from_title(affy_meta_1$title)

if (n_genes > 0 && n_genes < nrow(affy_1)) {
  log_affy_1 <- log2(pmax(affy_1, .Machine$double.eps))
  vars <- if (requireNamespace("matrixStats", quietly = TRUE)) {
    matrixStats::rowVars(log_affy_1)
  } else {
    apply(log_affy_1, 1, stats::var)
  }
  keep_genes <- order(vars, decreasing = TRUE)[seq_len(n_genes)]
  affy_fit <- affy_1[keep_genes, , drop = FALSE]
} else {
  affy_fit <- affy_1
}

qpcr_groups <- group_from_title(qpcr_meta$title)
colnames(qpcr) <- qpcr_meta$title
qpcr_p <- limma_qpcr_validation(qpcr, qpcr_groups)
qpcr_annot <- gpl4097[, c("ID", "ORF", "AssayID"), drop = FALSE]
names(qpcr_annot) <- c("qpcr_id", "symbol", "assay_id")
qpcr_annot$symbol <- toupper(trimws(qpcr_annot$symbol))
qpcr_info <- merge(
  data.frame(qpcr_id = rownames(qpcr), qpcr_p = qpcr_p, stringsAsFactors = FALSE),
  qpcr_annot,
  by = "qpcr_id",
  all.x = TRUE,
  sort = FALSE
)
qpcr_info$validated <- qpcr_info$qpcr_p < 0.05

affy_map <- make_affy_symbol_map(gpl570)

message("Fitting GaGa on ", nrow(affy_fit), " probes and ", ncol(affy_fit), " arrays.")
gg_fit <- fit_gaga_maqc(affy_fit, groups_affy, nclust = 1, method = gaga_method)
message("Fitting MiGaGa2 on ", nrow(affy_fit), " probes and ", ncol(affy_fit), " arrays.")
mg_fit <- fit_gaga_maqc(affy_fit, groups_affy, nclust = 2, method = gaga_method)
message("Fitting Ga/Gamma-Gamma on ", nrow(affy_fit), " probes and ", ncol(affy_fit), " arrays.")
ga_error <- NA_character_
ga_fit <- tryCatch({
  fit_ga_maqc(affy_fit, groups_affy)
}, error = function(e) {
  ga_error <<- conditionMessage(e)
  warning("Skipping Ga/Gamma-Gamma: ", ga_error)
  NULL
})

gaga_score <- 1 - gg_fit$pp[, 1]
migaga_score <- 1 - mg_fit$pp[, 1]
names(gaga_score) <- rownames(affy_fit)
names(migaga_score) <- rownames(affy_fit)
limma_score <- limma_f_score(affy_fit, groups_affy)
if (!is.null(ga_fit)) {
  ga_score <- 1 - ga_fit$pp[, 1]
  names(ga_score) <- rownames(affy_fit)
}

paper_counts <- data.frame(
  method = rep(c("GaGa", "MiGaGa2"), each = 5),
  pattern = rep(0:4, 2),
  paper_count = c(20272, 0, 1429, 3935, 29039, 16328, 0, 1323, 3697, 33327)
)

pattern_rows <- list()
if (!is.null(ga_fit)) {
  pattern_rows[["Ga"]] <- data.frame(method = "Ga", pattern = 0:4,
                                     count = as.integer(tabulate(max.col(ga_fit$pp, ties.method = "first"), nbins = 5)))
}
pattern_counts <- do.call(rbind, c(pattern_rows, list(
  data.frame(method = "GaGa", pattern = 0:4,
             count = as.integer(tabulate(max.col(gg_fit$pp, ties.method = "first"), nbins = 5))),
  data.frame(method = "MiGaGa2", pattern = 0:4,
             count = as.integer(tabulate(max.col(mg_fit$pp, ties.method = "first"), nbins = 5)))
)))
pattern_counts <- merge(pattern_counts, paper_counts, by = c("method", "pattern"), all.x = TRUE)
pattern_counts$diff_from_paper <- pattern_counts$count - pattern_counts$paper_count

method_scores <- list(
  GaGa = score_by_symbol(gaga_score, affy_map),
  MiGaGa2 = score_by_symbol(migaga_score, affy_map),
  limma_BH = score_by_symbol(limma_score, affy_map)
)
if (!is.null(ga_fit)) {
  method_scores <- c(list(Ga = score_by_symbol(ga_score, affy_map)), method_scores)
}

probe_symbols <- collapse_probe_symbols(affy_map)
top_pattern4 <- top_pattern4_lists(
  fits = list(GaGa = gg_fit, MiGaGa2 = mg_fit),
  affy_fit = affy_fit,
  groups = groups_affy,
  probe_symbols = probe_symbols,
  n = 1000
)

curve_rows <- list()
score_rows <- list()
summary_rows <- list()
for (method in names(method_scores)) {
  dat <- merge(qpcr_info, method_scores[[method]], by = "symbol", all.x = FALSE, sort = FALSE)
  curve <- validation_curve(dat$score, dat$validated)
  curve$method <- method
  curve_rows[[method]] <- curve
  dat$method <- method
  score_rows[[method]] <- dat
  final <- curve[nrow(curve), ]
  summary_rows[[method]] <- data.frame(
    method = method,
    n_qpcr_assays = nrow(qpcr_info),
    n_mapped_assays = nrow(dat),
    n_validated_mapped = sum(dat$validated),
    prop_nonvalidated_mapped = mean(!dat$validated),
    auc_scaled = final$auc_scaled,
    auc_standard = final$auc_standard,
    stringsAsFactors = FALSE
  )
}

data_summary <- data.frame(
  affy_file = affy_file,
  qpcr_file = qpcr_file,
  n_affy_probes = nrow(affy),
  n_affy_samples = ncol(affy),
  n_first_site_arrays = ncol(affy_1),
  first_site_titles = paste(affy_meta_1$title, collapse = ","),
  n_fit_probes = nrow(affy_fit),
  n_qpcr_assays = nrow(qpcr),
  n_qpcr_samples = ncol(qpcr),
  n_qpcr_validated = sum(qpcr_info$validated),
  prop_qpcr_nonvalidated = mean(!qpcr_info$validated),
  n_affy_symbol_links = nrow(affy_map),
  gaga_method_requested = gaga_method,
  ga_method_used = if (is.null(ga_fit)) NA_character_ else ga_fit$method_used,
  ga_error = ga_error,
  gaga_method_used = gg_fit$method_used,
  migaga_method_used = mg_fit$method_used,
  stringsAsFactors = FALSE
)

curve_all <- do.call(rbind, curve_rows)
score_all <- do.call(rbind, score_rows)
roc_summary <- do.call(rbind, summary_rows)

write.csv(data_summary, result_file("maqc_data_summary.csv"), row.names = FALSE)
write.csv(pattern_counts, result_file("maqc_pattern_counts.csv"), row.names = FALSE)
write.csv(roc_summary, result_file("maqc_roc_summary.csv"), row.names = FALSE)
write.csv(curve_all, result_file("maqc_roc_curves.csv"), row.names = FALSE)
write.csv(score_all, result_file("maqc_qpcr_scores.csv"), row.names = FALSE)
write.csv(top_pattern4, result_file("maqc_top_pattern4_genes.csv"), row.names = FALSE)

plot_maqc <- function(file, device = c("pdf", "png")) {
  device <- match.arg(device)
  if (device == "pdf") {
    grDevices::pdf(file, width = 6.5, height = 5)
  } else {
    grDevices::png(file, width = 1300, height = 1000, res = 200)
  }
  plot(
    NA,
    xlim = c(0, max(curve_all$prop_nonvalidated, na.rm = TRUE)),
    ylim = c(0, max(curve_all$prop_validated, na.rm = TRUE)),
    xlab = "Proportion nonvalidated qPCR assays",
    ylab = "Proportion validated qPCR assays",
    main = "MAQC validation curves"
  )
  cols <- c(Ga = "#666666", GaGa = "#1B9E77", MiGaGa2 = "#D95F02", limma_BH = "#7570B3")
  for (method in names(method_scores)) {
    cc <- curve_all[curve_all$method == method, ]
    lines(cc$prop_nonvalidated, cc$prop_validated, col = cols[[method]], lwd = 2)
  }
  leg <- merge(data.frame(method = names(method_scores)), roc_summary, by = "method", sort = FALSE)
  legend("bottomright",
         legend = paste0(leg$method, " AUC=", sprintf("%.4f", leg$auc_scaled)),
         col = cols[leg$method], lwd = 2, bty = "n")
  grDevices::dev.off()
}

plot_maqc(result_file("maqc_roc_curves.pdf"), "pdf")
plot_maqc(result_file("maqc_roc_curves.png"), "png")

print(data_summary)
print(pattern_counts)
print(roc_summary)
