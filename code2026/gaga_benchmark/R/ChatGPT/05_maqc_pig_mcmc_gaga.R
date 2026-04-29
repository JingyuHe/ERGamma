args <- commandArgs(FALSE)
this_file <- sub("^--file=", "", args[grepl("^--file=", args)][1])
src_dir <- dirname(normalizePath(this_file, mustWork = TRUE))
source(file.path(src_dir, "quad_gaga_lib.R"))
source(file.path(src_dir, "pig_mcmc_gaga_lib.R"))

bench <- setup_paths()
require_pkgs(c("gaga", "GIGrvg", "Biobase"))

fit_gaga_fallback <- function(x, groups, patterns, method) {
  methods <- unique(c(method, "quickEM", "EM"))
  errors <- character()
  for (m in methods) {
    fit <- tryCatch({
      fit_gaga_original(x, groups, patterns, nclust = 1, method = m, fdr = 0.05)
    }, error = function(e) {
      errors <<- c(errors, paste0(m, ": ", conditionMessage(e)))
      NULL
    })
    if (!is.null(fit)) {
      fit$method_used <- m
      return(fit)
    }
  }
  stop("GaGa fit failed: ", paste(errors, collapse = " | "))
}

plot_validation_curves <- function(curves, out_pdf, out_png) {
  methods <- unique(curves$method)
  cols <- c(GaGa = "#1B9E77", PIG_MCMC_GaGa = "#D95F02", limma_F = "#7570B3")
  draw <- function() {
    plot(NA, NA, xlim = c(0, max(curves$prop_nonvalidated, na.rm = TRUE)),
         ylim = c(0, max(curves$prop_validated, na.rm = TRUE)),
         xlab = "Proportion non-validated qPCR assays",
         ylab = "Proportion validated qPCR assays")
    abline(0, 1, col = "grey80", lty = 2)
    for (m in methods) {
      cur <- curves[curves$method == m, ]
      lines(cur$prop_nonvalidated, cur$prop_validated,
            col = cols[[m]], lwd = 2)
    }
    legend("bottomright", legend = methods, col = cols[methods], lwd = 2, bty = "n")
  }
  grDevices::pdf(out_pdf, width = 6, height = 5)
  draw()
  grDevices::dev.off()
  grDevices::png(out_png, width = 1500, height = 1200, res = 220)
  draw()
  grDevices::dev.off()
}

method <- Sys.getenv("PIG_MCMC_MAQC_FIT_METHOD", unset = "EM")
max_probes <- env_int("PIG_MCMC_MAQC_MAX_PROBES", 1000)
fdr <- env_num("PIG_MCMC_MAQC_FDR", 0.05)
progress_every <- env_int("PIG_MCMC_MAQC_PROGRESS_EVERY", 250)
seed <- env_int("PIG_MCMC_MAQC_SEED", 2026)
n_iter <- env_int("PIG_MCMC_ITER", env_int("PIG_MCMC_MAQC_ITER", 500))
burnin <- env_int("PIG_MCMC_BURNIN", env_int("PIG_MCMC_MAQC_BURNIN", 200))
thin <- env_int("PIG_MCMC_THIN", env_int("PIG_MCMC_MAQC_THIN", 1))
trunc <- env_int("PIG_MCMC_TRUNC", env_int("PIG_MCMC_MAQC_TRUNC", 60))
alpha_steps <- env_int("PIG_MCMC_ALPHA_STEPS", env_int("PIG_MCMC_MAQC_ALPHA_STEPS", 1))
alpha_kernel <- Sys.getenv("PIG_MCMC_ALPHA_KERNEL", unset = "slice")
slice_width <- env_num("PIG_MCMC_SLICE_WIDTH", env_num("PIG_MCMC_MAQC_SLICE_WIDTH", 0.5))
slice_max_steps <- env_int("PIG_MCMC_SLICE_MAX_STEPS", env_int("PIG_MCMC_MAQC_SLICE_MAX_STEPS", 50))
set.seed(seed)

data_dir <- file.path(bench, "data")
affy_file <- file.path(data_dir, "GSE5350-GPL570_series_matrix.txt.gz")
qpcr_file <- file.path(data_dir, "GSE5350-GPL4097_TaqMan_series_matrix.txt.gz")
gpl570_file <- file.path(data_dir, "GPL570.annot.gz")
gpl4097_file <- file.path(data_dir, "GPL4097_family.soft.gz")

cat("Reading MAQC matrices and annotations\n")
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
groups_affy <- factor(group_from_title(affy_meta_1$title),
                      levels = c("A", "C", "D", "B"))

qpcr_groups <- factor(group_from_title(qpcr_meta$title),
                      levels = c("A", "C", "D", "B"))
colnames(qpcr) <- qpcr_meta$title
qpcr_p <- limma_qpcr_validation(qpcr, qpcr_groups)
qpcr_annot <- gpl4097[, c("ID", "ORF", "AssayID"), drop = FALSE]
names(qpcr_annot) <- c("qpcr_id", "symbol", "assay_id")
qpcr_annot$symbol <- toupper(trimws(qpcr_annot$symbol))
qpcr_info <- merge(
  data.frame(qpcr_id = rownames(qpcr), qpcr_p = qpcr_p,
             stringsAsFactors = FALSE),
  qpcr_annot,
  by = "qpcr_id",
  all.x = TRUE,
  sort = FALSE
)
qpcr_info$validated <- qpcr_info$qpcr_p < 0.05
qpcr_symbols <- unique(qpcr_info$symbol[nzchar(qpcr_info$symbol)])

affy_map <- make_affy_symbol_map(gpl570)
mapped_probes <- unique(affy_map$probe_id[affy_map$symbol %in% qpcr_symbols])
mapped_probes <- intersect(mapped_probes, rownames(affy_1))
affy_fit <- affy_1[mapped_probes, , drop = FALSE]

if (max_probes > 0 && max_probes < nrow(affy_fit)) {
  vars <- apply(log2(pmax(affy_fit, .Machine$double.eps)), 1, stats::var)
  keep_probes <- order(vars, decreasing = TRUE)[seq_len(max_probes)]
  affy_fit <- affy_fit[keep_probes, , drop = FALSE]
}

affy_map_fit <- affy_map[affy_map$probe_id %in% rownames(affy_fit), , drop = FALSE]
patterns <- maqc_patterns()

cat("MAQC PIG-MCMC-GaGa benchmark\n")
cat("  probes=", nrow(affy_fit),
    " arrays=", ncol(affy_fit),
    " mapped_qpcr_symbols=", length(qpcr_symbols),
    " iter=", n_iter,
    " burnin=", burnin,
    " trunc=", trunc, "\n", sep = "")
cat("  alpha_kernel=", alpha_kernel,
    " slice_width=", slice_width,
    " slice_max_steps=", slice_max_steps, "\n", sep = "")
cat("  qPCR validation p-values from ",
    if (requireNamespace("limma", quietly = TRUE)) "limma" else "base F test",
    "\n", sep = "")

elapsed_gaga <- system.time({
  gg_fit <- fit_gaga_fallback(affy_fit, groups_affy, patterns, method = method)
})[["elapsed"]]

elapsed_pig <- system.time({
  pig <- pig_mcmc_gaga_pp(
    affy_fit, groups_affy, patterns, gg_fit,
    n_iter = n_iter,
    burnin = burnin,
    thin = thin,
    trunc = trunc,
    alpha_steps = alpha_steps,
    alpha_kernel = alpha_kernel,
    slice_width = slice_width,
    slice_max_steps = slice_max_steps,
    progress_every = progress_every,
    seed = seed
  )
})[["elapsed"]]

gaga_score <- 1 - gg_fit$pp[, 1]
pig_score <- 1 - pig$pp[, 1]
limma_p <- limma_f_pvalue(affy_fit, groups_affy)
limma_score <- -log10(pmax(limma_p, 1e-300))
names(gaga_score) <- rownames(affy_fit)
names(pig_score) <- rownames(affy_fit)

pig_calls <- find_pig_mcmc_calls(pig$pp, fdr = fdr)
gaga_calls <- as.logical(gaga::findgenes(gg_fit, affy_fit, as.character(groups_affy),
                                         fdrmax = fdr, parametric = TRUE)$d)

probe_scores <- data.frame(
  probe_id = rownames(affy_fit),
  gaga_score = gaga_score,
  pig_mcmc_score = pig_score,
  limma_score = limma_score,
  gaga_pp_null = gg_fit$pp[, 1],
  pig_mcmc_pp_null = pig$pp[, 1],
  gaga_call = gaga_calls,
  pig_mcmc_call = pig_calls,
  pig_alpha_mean = pig$alpha_mean,
  pig_accept_rate = pig$accept_rate,
  stringsAsFactors = FALSE
)

pattern_counts <- rbind(
  data.frame(method = "GaGa", pattern = 0:4,
             count = as.integer(tabulate(max.col(gg_fit$pp, ties.method = "first"), nbins = 5))),
  data.frame(method = "PIG_MCMC_GaGa", pattern = 0:4,
             count = as.integer(tabulate(max.col(pig$pp, ties.method = "first"), nbins = 5)))
)

method_scores <- list(
  GaGa = score_by_symbol(gaga_score, affy_map_fit),
  PIG_MCMC_GaGa = score_by_symbol(pig_score, affy_map_fit),
  limma_F = score_by_symbol(limma_score, affy_map_fit)
)

curve_rows <- list()
score_rows <- list()
summary_rows <- list()
for (m in names(method_scores)) {
  dat <- merge(qpcr_info, method_scores[[m]], by = "symbol",
               all.x = FALSE, sort = FALSE)
  curve <- validation_curve(dat$score, dat$validated)
  curve$method <- m
  curve_rows[[m]] <- curve
  dat$method <- m
  score_rows[[m]] <- dat
  final <- curve[nrow(curve), ]
  summary_rows[[m]] <- data.frame(
    method = m,
    n_fit_probes = nrow(affy_fit),
    n_qpcr_assays = nrow(qpcr_info),
    n_mapped_assays = nrow(dat),
    n_validated_mapped = sum(dat$validated),
    auc_scaled = final$auc_scaled,
    auc_standard = final$auc_standard,
    spearman_probe_score_vs_gaga = if (m == "GaGa") 1 else
      stats::cor(probe_scores$gaga_score,
                 if (m == "PIG_MCMC_GaGa") probe_scores$pig_mcmc_score else probe_scores$limma_score,
                 method = "spearman", use = "pairwise.complete.obs"),
    n_de_at_fdr = if (m == "GaGa") sum(gaga_calls) else if (m == "PIG_MCMC_GaGa") sum(pig_calls) else NA_integer_,
    stringsAsFactors = FALSE
  )
}

curves <- do.call(rbind, curve_rows)
validation_scores <- do.call(rbind, score_rows)
summary <- do.call(rbind, summary_rows)
summary$gaga_fit_method <- gg_fit$method_used
summary$elapsed_gaga_seconds <- elapsed_gaga
summary$elapsed_pig_seconds <- elapsed_pig
summary$pig_accept_rate_mean <- mean(pig$accept_rate, na.rm = TRUE)
summary$pig_accept_rate_median <- stats::median(pig$accept_rate, na.rm = TRUE)
summary$n_iter <- n_iter
summary$burnin <- burnin
summary$thin <- thin
summary$trunc <- trunc
summary$alpha_kernel <- alpha_kernel
summary$slice_width <- slice_width
summary$slice_max_steps <- slice_max_steps
summary$qPCR_validation_engine <- if (requireNamespace("limma", quietly = TRUE)) "limma" else "base_F"

write.csv(summary, chatgpt_result_file("maqc_pig_mcmc_gaga_validation_summary.csv"),
          row.names = FALSE)
write.csv(curves, chatgpt_result_file("maqc_pig_mcmc_gaga_validation_curve.csv"),
          row.names = FALSE)
write.csv(validation_scores, chatgpt_result_file("maqc_pig_mcmc_gaga_validation_scores.csv"),
          row.names = FALSE)
write.csv(probe_scores, chatgpt_result_file("maqc_pig_mcmc_gaga_probe_scores.csv"),
          row.names = FALSE)
write.csv(pattern_counts, chatgpt_result_file("maqc_pig_mcmc_gaga_pattern_counts.csv"),
          row.names = FALSE)
saveRDS(
  list(config = list(max_probes = max_probes, fdr = fdr,
                     fit_method = method, seed = seed,
                     n_iter = n_iter, burnin = burnin, thin = thin,
                     trunc = trunc, alpha_steps = alpha_steps,
                     alpha_kernel = alpha_kernel,
                     slice_width = slice_width,
                     slice_max_steps = slice_max_steps),
       summary = summary, curves = curves,
       validation_scores = validation_scores,
       probe_scores = probe_scores,
       pattern_counts = pattern_counts,
       fit = gg_fit, pig = pig),
  chatgpt_result_file("maqc_pig_mcmc_gaga.rds")
)

plot_validation_curves(
  curves,
  out_pdf = chatgpt_result_file("maqc_pig_mcmc_gaga_validation_curve.pdf"),
  out_png = chatgpt_result_file("maqc_pig_mcmc_gaga_validation_curve.png")
)

print(summary)
print(pattern_counts)
