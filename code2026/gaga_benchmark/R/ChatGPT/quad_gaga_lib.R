script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- args[grepl("^--file=", args)]
  if (length(file_arg) > 0) {
    return(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE))
  }
  normalizePath(getwd(), mustWork = TRUE)
}

chatgpt_dir <- function() {
  normalizePath(dirname(script_path()), mustWork = TRUE)
}

benchmark_dir <- function() {
  normalizePath(file.path(chatgpt_dir(), "..", ".."), mustWork = TRUE)
}

chatgpt_result_file <- function(name) {
  out_dir <- file.path(chatgpt_dir(), "results")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  file.path(out_dir, name)
}

setup_paths <- function() {
  bench <- benchmark_dir()
  paths <- .libPaths()
  system_paths <- paths[grepl("^/opt/homebrew", paths)]
  local_paths <- paths[!grepl("^/opt/homebrew", paths)]
  .libPaths(unique(c(
    system_paths,
    file.path(bench, "r-lib-gcrma"),
    file.path(bench, "r-lib"),
    file.path(bench, "..", "r-lib"),
    local_paths
  )))
  invisible(bench)
}

require_pkgs <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop("Missing packages: ", paste(missing, collapse = ", "), call. = FALSE)
  }
}

env_int <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) {
    return(default)
  }
  as.integer(value)
}

env_num <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) {
    return(default)
  }
  as.numeric(value)
}

safe_div <- function(num, den) {
  ifelse(den == 0, 0, num / den)
}

softmax_log <- function(logv) {
  m <- max(logv, na.rm = TRUE)
  z <- exp(logv - m)
  z / sum(z)
}

gamma_const <- -digamma(1)

rPIG_ERGamma <- function(n, cc, trunc) {
  require_pkgs("GIGrvg")
  output <- rep(0, n)
  if (trunc > 0) {
    for (i in seq_len(trunc)) {
      if (cc == 0) {
        output <- output + 1 / stats::rgamma(n, shape = 3 / 2) / (4 * i^2)
      } else {
        output <- output + GIGrvg::rgig(
          n = n,
          lambda = -3 / 2,
          chi = 1 / (2 * i^2),
          psi = 2 * cc^2
        )
      }
    }
  }
  if (cc == 0) {
    temp1 <- trigamma(1 + trunc)
    temp2 <- psigamma(1 + trunc, deriv = 2)
    rate <- -4 * temp1 / temp2
    shape <- 0.5 * temp1 * rate
  } else {
    temp1 <- digamma(1 + trunc + cc) - digamma(1 + trunc)
    temp2 <- temp1 - cc * trigamma(1 + trunc + cc)
    rate <- 2 * cc^2 * temp1 / temp2
    shape <- rate / (2 * cc) * temp1
  }
  output + stats::rgamma(n, shape = shape, rate = rate)
}

rPTN <- function(power, quad, lin, max_try = 5000) {
  if (!is.finite(power) || power <= 0 || !is.finite(quad) || quad <= 0) {
    stop("Invalid PTN parameters.")
  }
  if (!is.finite(lin) || abs(lin) < 1e-12) {
    lin <- 1e-12
  }
  tau <- if (lin > 0) {
    0.5 + sqrt(0.25 + 2 * quad * power / lin^2)
  } else {
    -0.5 + sqrt(0.25 + 2 * quad * power / lin^2)
  }
  rate <- tau * abs(lin) - lin
  center <- tau * abs(lin) / (2 * quad)
  for (i in seq_len(max_try)) {
    x <- stats::rgamma(1, shape = power, rate = rate)
    if (stats::runif(1) <= exp(-quad * (x - center)^2)) {
      return(x)
    }
  }
  NA_real_
}

pig_gas_sampler <- function(beta, c_lin, n_draws = 1000, burnin = 500,
                            trunc = 200, init = 5) {
  beta <- as.numeric(beta)
  if (any(beta <= 0)) {
    stop("beta entries must be positive.")
  }
  s <- length(beta)
  total <- n_draws + burnin
  alpha <- init
  out <- numeric(n_draws)
  for (iter in seq_len(total)) {
    tau <- stats::rgamma(s, shape = alpha + beta, rate = 1)
    omega <- rPIG_ERGamma(s, alpha, trunc)
    quad <- sum(omega)
    lin <- s * gamma_const + c_lin + sum(log(tau))
    alpha <- rPTN(s + 1, quad, lin)
    if (!is.finite(alpha)) {
      stop("PTN rejection sampler failed.")
    }
    if (iter > burnin) {
      out[iter - burnin] <- alpha
    }
  }
  out
}

load_armstrong_schliep <- function(data_file, exclude_last_mll = TRUE) {
  raw <- read.delim(data_file, check.names = FALSE, stringsAsFactors = FALSE)
  groups <- as.character(unlist(raw[1, -1]))
  x <- as.matrix(data.frame(lapply(raw[-1, -1, drop = FALSE], as.numeric), check.names = FALSE))
  rownames(x) <- raw[-1, 1]
  sample_ids <- paste0(groups, "_", ave(seq_along(groups), groups, FUN = seq_along))
  colnames(x) <- sample_ids
  keep <- groups %in% c("ALL", "MLL")
  if (exclude_last_mll) {
    mll_idx <- which(groups == "MLL")
    keep[tail(mll_idx, 2)] <- FALSE
  }
  list(x = x[, keep, drop = FALSE], groups = factor(groups[keep], levels = c("ALL", "MLL")))
}

load_armstrong_orange <- function(data_file, exclude_last_mll = TRUE) {
  raw <- read.delim(data_file, check.names = FALSE, stringsAsFactors = FALSE)
  groups <- raw$class[-c(1, 2)]
  xdf <- raw[-c(1, 2), -1, drop = FALSE]
  x <- t(apply(xdf, 2, as.numeric))
  rownames(x) <- colnames(raw)[-1]
  sample_ids <- paste0(groups, "_", ave(seq_along(groups), groups, FUN = seq_along))
  colnames(x) <- sample_ids
  keep <- groups %in% c("ALL", "MLL")
  if (exclude_last_mll) {
    mll_idx <- which(groups == "MLL")
    keep[tail(mll_idx, 2)] <- FALSE
  }
  list(x = x[, keep, drop = FALSE], groups = factor(groups[keep], levels = c("ALL", "MLL")))
}

transform_expr <- function(x, transform) {
  if (transform == "log") {
    return(log(pmax(x, 1)))
  }
  if (transform == "log1p") {
    return(log(pmax(x, 0) + 1))
  }
  if (transform == "log2_floor20") {
    return(log2(pmax(x, 20)))
  }
  if (transform == "log2p1") {
    return(log2(pmax(x, 0) + 1))
  }
  if (transform == "none") {
    if (min(x, na.rm = TRUE) <= 0) {
      stop("Expression matrix has non-positive values.")
    }
    return(x)
  }
  stop("Unknown transform: ", transform)
}

gaga_patterns_two_group <- function(groups) {
  patterns <- matrix(c(0, 0, 0, 1), 2, 2)
  colnames(patterns) <- levels(factor(groups))
  patterns
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

fit_gaga_original <- function(x, groups, patterns, nclust = 1,
                              method = "quickEM", fdr = 0.05) {
  require_pkgs("gaga")
  groups_fit <- as.character(groups)
  fit <- gaga::fitGG(
    x,
    groups_fit,
    patterns = patterns,
    equalcv = TRUE,
    nclust = nclust,
    method = method,
    trace = FALSE
  )
  gaga::parest(fit, x = x, groups = groups_fit, alpha = fdr)
}

block_stats <- function(xrow, groups, pattern_row) {
  labels <- pattern_row[as.character(groups)]
  blocks <- split(seq_along(groups), labels)
  lapply(blocks, function(idx) {
    vals <- pmax(xrow[idx], .Machine$double.eps)
    list(n = length(vals), sum = sum(vals), log_sum = sum(log(vals)))
  })
}

collapsed_alpha_log_kernel <- function(alpha, stats, par) {
  if (!is.finite(alpha) || alpha <= 0) {
    return(-Inf)
  }
  a0 <- as.numeric(par$a0[1])
  nu <- as.numeric(par$nu[1])
  beta <- as.numeric(par$balpha[1])
  mu <- as.numeric(par$nualpha[1])
  rate_alpha <- beta / mu
  r <- a0 / nu
  out <- beta * log(rate_alpha) - lgamma(beta) +
    (beta - 1) * log(alpha) - rate_alpha * alpha
  for (st in stats) {
    nb <- st$n
    sx <- st$sum
    lx <- st$log_sum
    out <- out +
      a0 * log(r) - lgamma(a0) +
      nb * alpha * log(alpha) -
      nb * lgamma(alpha) +
      (alpha - 1) * lx +
      lgamma(a0 + nb * alpha) -
      (a0 + nb * alpha) * log(r + alpha * sx)
  }
  out
}

log_integral_alpha <- function(stats, par, lower = -8, upper = log(5000)) {
  logf <- function(t) {
    vapply(t, function(tt) {
      alpha <- exp(tt)
      collapsed_alpha_log_kernel(alpha, stats, par) + tt
    }, numeric(1))
  }
  opt <- stats::optimize(function(t) -logf(t), interval = c(lower, upper))
  mode_t <- opt$minimum
  log_scale <- logf(mode_t)
  value <- tryCatch({
    stats::integrate(
      function(t) exp(logf(t) - log_scale),
      lower,
      upper,
      subdivisions = 200,
      rel.tol = 1e-5
    )$value
  }, error = function(e) NA_real_)
  if (!is.finite(value) || value <= 0) {
    return(NA_real_)
  }
  log(value) + log_scale
}

quad_gaga_pp <- function(x, groups, patterns, fit, max_genes = 0,
                         progress_every = 500) {
  par <- gaga::getpar(fit)
  probpat <- as.numeric(par$probpat)
  if (length(probpat) != nrow(patterns)) {
    probpat <- rep(1 / nrow(patterns), nrow(patterns))
  }
  if (max_genes > 0 && max_genes < nrow(x)) {
    x <- x[seq_len(max_genes), , drop = FALSE]
  }
  out <- matrix(NA_real_, nrow = nrow(x), ncol = nrow(patterns))
  rownames(out) <- rownames(x)
  colnames(out) <- paste0("pattern", seq_len(nrow(patterns)) - 1)
  log_marg <- out
  for (i in seq_len(nrow(x))) {
    if (progress_every > 0 && i %% progress_every == 0) {
      cat("  Quad-GaGa scoring gene ", i, "/", nrow(x), "\n", sep = "")
    }
    logp <- numeric(nrow(patterns))
    for (h in seq_len(nrow(patterns))) {
      stats <- block_stats(x[i, ], groups, patterns[h, ])
      li <- log_integral_alpha(stats, par)
      logp[h] <- log(probpat[h]) + li
      log_marg[i, h] <- li
    }
    out[i, ] <- softmax_log(logp)
  }
  list(pp = out, log_marginal = log_marg, par = par)
}

pig_gaga_pp <- function(...) {
  stop(
    "Deprecated: this empirical path uses adaptive quadrature, not PIG augmentation. ",
    "Use quad_gaga_pp() and label the method Quad_GaGa.",
    call. = FALSE
  )
}

find_quad_calls <- function(pp, fdr = 0.05) {
  null <- pp[, 1]
  ord <- order(null)
  fdr_path <- cumsum(null[ord]) / seq_along(ord)
  keep_ord <- which(fdr_path <= fdr)
  calls <- rep(FALSE, nrow(pp))
  if (length(keep_ord) > 0) {
    calls[ord[seq_len(max(keep_ord))]] <- TRUE
  }
  calls
}

find_pig_calls <- function(...) {
  stop(
    "Deprecated: use find_quad_calls(); the old name mislabeled quadrature results as PIG.",
    call. = FALSE
  )
}

score_by_symbol <- function(scores, affy_map) {
  idx <- match(affy_map$probe_id, names(scores))
  map <- affy_map[!is.na(idx), , drop = FALSE]
  map$score <- scores[idx[!is.na(idx)]]
  stats::aggregate(score ~ symbol, data = map, FUN = max, na.rm = TRUE)
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

read_text_lines <- function(file) {
  con <- if (grepl("\\.gz$", file)) gzfile(file, "rt") else file(file, "rt")
  on.exit(close(con), add = TRUE)
  readLines(con, warn = FALSE)
}

read_series_metadata <- function(file) {
  lines <- read_text_lines(file)
  geo_accession <- sub('^!Sample_geo_accession\\t"', "", grep("^!Sample_geo_accession", lines, value = TRUE))
  title <- sub('^!Sample_title\\t"', "", grep("^!Sample_title", lines, value = TRUE))
  gsm <- strsplit(gsub('"', "", geo_accession), "\t", fixed = TRUE)[[1]]
  ttl <- strsplit(gsub('"', "", title), "\t", fixed = TRUE)[[1]]
  data.frame(gsm = gsm, title = ttl, stringsAsFactors = FALSE)
}

read_series_matrix <- function(file) {
  lines <- read_text_lines(file)
  start <- grep("!series_matrix_table_begin", lines)
  end <- grep("!series_matrix_table_end", lines)
  txt <- lines[(start + 1):(end - 1)]
  tab <- read.delim(text = paste(txt, collapse = "\n"), check.names = FALSE,
                    stringsAsFactors = FALSE)
  ids <- tab[[1]]
  mat <- as.matrix(data.frame(lapply(tab[-1], as.numeric), check.names = FALSE))
  rownames(mat) <- ids
  mat
}

read_platform_table <- function(file, gz = TRUE) {
  con <- if (gz) gzfile(file, "rt") else file(file, "rt")
  on.exit(close(con))
  lines <- readLines(con, warn = FALSE)
  start <- grep("!platform_table_begin", lines)
  end <- grep("!platform_table_end", lines)
  read.delim(text = paste(lines[(start + 1):(end - 1)], collapse = "\n"),
             check.names = FALSE, stringsAsFactors = FALSE, quote = "")
}

split_symbols <- function(x) {
  lapply(strsplit(as.character(x), "///|//|;|,|\\|"), function(z) {
    z <- toupper(trimws(z))
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

base_f_pvalue <- function(x, groups) {
  groups <- factor(groups)
  design <- stats::model.matrix(~ groups)
  p <- apply(log2(pmax(x, .Machine$double.eps)), 1, function(y) {
    fit <- stats::lm.fit(design, y)
    rss1 <- sum(fit$residuals^2)
    rss0 <- sum((y - mean(y))^2)
    df1 <- ncol(design) - 1
    df2 <- length(y) - ncol(design)
    if (df2 <= 0 || rss1 <= 0) {
      return(NA_real_)
    }
    f <- ((rss0 - rss1) / df1) / (rss1 / df2)
    stats::pf(f, df1 = df1, df2 = df2, lower.tail = FALSE)
  })
  names(p) <- rownames(x)
  p
}

limma_f_pvalue <- function(x, groups) {
  if (!requireNamespace("limma", quietly = TRUE)) {
    return(base_f_pvalue(x, groups))
  }
  groups <- factor(groups)
  design <- stats::model.matrix(~ groups)
  fit <- limma::lmFit(log2(pmax(x, .Machine$double.eps)), design)
  fit <- limma::eBayes(fit)
  tt <- limma::topTable(fit, coef = seq_len(ncol(design))[-1], number = Inf,
                        sort.by = "none")
  p <- tt$P.Value
  names(p) <- rownames(tt)
  p
}

limma_qpcr_validation <- function(qpcr, groups) {
  limma_f_pvalue(qpcr, groups)
}

group_from_title <- function(title) {
  sub("^MAQC_[A-Z]+_[0-9]+_([ABCD]).*$", "\\1", title)
}

first_site_filter <- function(meta) {
  grepl("^MAQC_AFX_1_[ABCD][1-5]$", meta$title)
}
