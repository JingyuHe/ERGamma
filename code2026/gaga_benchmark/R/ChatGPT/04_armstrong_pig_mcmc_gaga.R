args <- commandArgs(FALSE)
this_file <- sub("^--file=", "", args[grepl("^--file=", args)][1])
src_dir <- dirname(normalizePath(this_file, mustWork = TRUE))
source(file.path(src_dir, "quad_gaga_lib.R"))
source(file.path(src_dir, "pig_mcmc_gaga_lib.R"))

bench <- setup_paths()
require_pkgs(c("gaga", "GIGrvg", "Biobase"))

fit_gaga_calls <- function(x, groups, patterns, fdr, method) {
  fit <- fit_gaga_original(x, groups, patterns, nclust = 1,
                           method = method, fdr = fdr)
  fg <- gaga::findgenes(fit, x, as.character(groups), fdrmax = fdr,
                        parametric = TRUE)
  list(
    fit = fit,
    calls = as.logical(fg$d),
    score = 1 - fit$pp[, 1],
    threshold = fg$threshold,
    fdr = fg$fdr
  )
}

fit_pig_mcmc_calls <- function(x, groups, patterns, gaga_fit, fdr,
                               n_iter, burnin, thin, trunc,
                               alpha_steps, alpha_kernel, slice_width,
                               slice_max_steps, progress_every, seed) {
  pig <- pig_mcmc_gaga_pp(
    x, groups, patterns, gaga_fit,
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
  list(
    fit = pig,
    calls = find_pig_mcmc_calls(pig$pp, fdr = fdr),
    score = 1 - pig$pp[, 1]
  )
}

summarize_pair <- function(gaga, pig, data_choice, transform_choice,
                           elapsed_gaga, elapsed_pig) {
  gaga_calls <- gaga$calls
  pig_calls <- pig$calls
  gaga_score <- gaga$score
  pig_score <- pig$score
  overlap <- sum(gaga_calls & pig_calls)
  union <- sum(gaga_calls | pig_calls)
  data.frame(
    method = c("GaGa", "PIG_MCMC_GaGa"),
    n_de = c(sum(gaga_calls), sum(pig_calls)),
    overlap_with_other = c(overlap, overlap),
    jaccard_with_other = c(safe_div(overlap, union), safe_div(overlap, union)),
    spearman_score_vs_other = c(
      stats::cor(gaga_score, pig_score, method = "spearman"),
      stats::cor(pig_score, gaga_score, method = "spearman")
    ),
    mean_abs_score_diff = c(
      mean(abs(gaga_score - pig_score)),
      mean(abs(pig_score - gaga_score))
    ),
    elapsed_seconds = c(elapsed_gaga, elapsed_pig),
    pig_accept_rate_mean = c(NA_real_, mean(pig$fit$accept_rate, na.rm = TRUE)),
    pig_accept_rate_median = c(NA_real_, stats::median(pig$fit$accept_rate, na.rm = TRUE)),
    data = data_choice,
    transform = transform_choice,
    stringsAsFactors = FALSE
  )
}

data_choice <- Sys.getenv("PIG_MCMC_ARM_DATA", unset = "schliep_filtered")
transform_choice <- Sys.getenv("PIG_MCMC_ARM_TRANSFORM", unset = "log")
exclude_last_mll <- Sys.getenv("PIG_MCMC_ARM_EXCLUDE_LAST_MLL", unset = "1") != "0"
fdr <- env_num("PIG_MCMC_ARM_FDR", 0.05)
method <- Sys.getenv("PIG_MCMC_ARM_FIT_METHOD", unset = "quickEM")
seed <- env_int("PIG_MCMC_ARM_SEED", 2026)
max_genes <- env_int("PIG_MCMC_ARM_MAX_GENES", 0)
n_iter <- env_int("PIG_MCMC_ITER", env_int("PIG_MCMC_ARM_ITER", 600))
burnin <- env_int("PIG_MCMC_BURNIN", env_int("PIG_MCMC_ARM_BURNIN", 200))
thin <- env_int("PIG_MCMC_THIN", env_int("PIG_MCMC_ARM_THIN", 1))
trunc <- env_int("PIG_MCMC_TRUNC", env_int("PIG_MCMC_ARM_TRUNC", 80))
alpha_steps <- env_int("PIG_MCMC_ALPHA_STEPS", env_int("PIG_MCMC_ARM_ALPHA_STEPS", 1))
alpha_kernel <- Sys.getenv("PIG_MCMC_ALPHA_KERNEL", unset = "slice")
slice_width <- env_num("PIG_MCMC_SLICE_WIDTH", env_num("PIG_MCMC_ARM_SLICE_WIDTH", 0.5))
slice_max_steps <- env_int("PIG_MCMC_SLICE_MAX_STEPS", env_int("PIG_MCMC_ARM_SLICE_MAX_STEPS", 50))
progress_every <- env_int("PIG_MCMC_ARM_PROGRESS_EVERY", 250)

data_dir <- file.path(bench, "data")
if (data_choice == "schliep_filtered") {
  data_file <- file.path(data_dir, "armstrong-2002-v2_database.txt")
  dat <- load_armstrong_schliep(data_file, exclude_last_mll = exclude_last_mll)
} else if (data_choice == "orange_full") {
  data_file <- file.path(data_dir, "MLL_orange_12533.tab")
  dat <- load_armstrong_orange(data_file, exclude_last_mll = exclude_last_mll)
} else {
  stop("Unknown PIG_MCMC_ARM_DATA: ", data_choice)
}

x <- transform_expr(dat$x, transform_choice)
if (max_genes > 0 && max_genes < nrow(x)) {
  vars <- apply(x, 1, stats::var)
  keep <- order(vars, decreasing = TRUE)[seq_len(max_genes)]
  x <- x[keep, , drop = FALSE]
}
groups <- droplevels(dat$groups)
patterns <- gaga_patterns_two_group(groups)

cat("Armstrong PIG-MCMC-GaGa benchmark\n")
cat("  data=", data_choice,
    " genes=", nrow(x),
    " samples=", ncol(x),
    " iter=", n_iter,
    " burnin=", burnin,
    " trunc=", trunc,
    " alpha_kernel=", alpha_kernel,
    " fdr=", fdr, "\n", sep = "")

elapsed_gaga <- system.time({
  gaga_res <- fit_gaga_calls(x, groups, patterns, fdr = fdr, method = method)
})[["elapsed"]]

elapsed_pig <- system.time({
  pig_res <- fit_pig_mcmc_calls(
    x, groups, patterns, gaga_res$fit, fdr = fdr,
    n_iter = n_iter, burnin = burnin, thin = thin, trunc = trunc,
    alpha_steps = alpha_steps, alpha_kernel = alpha_kernel,
    slice_width = slice_width, slice_max_steps = slice_max_steps,
    progress_every = progress_every, seed = seed
  )
})[["elapsed"]]

names(gaga_res$calls) <- rownames(x)
names(pig_res$calls) <- rownames(x)
names(gaga_res$score) <- rownames(x)
names(pig_res$score) <- rownames(x)

summary <- summarize_pair(
  gaga_res, pig_res, data_choice, transform_choice,
  elapsed_gaga = elapsed_gaga, elapsed_pig = elapsed_pig
)
summary$n_genes <- nrow(x)
summary$n_samples <- ncol(x)
summary$n_iter <- n_iter
summary$burnin <- burnin
summary$thin <- thin
summary$trunc <- trunc
summary$alpha_kernel <- alpha_kernel
summary$slice_width <- slice_width
summary$slice_max_steps <- slice_max_steps

scores <- data.frame(
  gene = rownames(x),
  gaga_score = gaga_res$score,
  pig_mcmc_score = pig_res$score,
  gaga_pp_null = gaga_res$fit$pp[, 1],
  pig_mcmc_pp_null = pig_res$fit$pp[, 1],
  gaga_call = gaga_res$calls,
  pig_mcmc_call = pig_res$calls,
  pig_alpha_mean = pig_res$fit$alpha_mean,
  pig_accept_rate = pig_res$fit$accept_rate,
  stringsAsFactors = FALSE
)
scores$abs_score_diff <- abs(scores$gaga_score - scores$pig_mcmc_score)

write.csv(summary, chatgpt_result_file("armstrong_pig_mcmc_gaga_summary.csv"),
          row.names = FALSE)
write.csv(scores, chatgpt_result_file("armstrong_pig_mcmc_gaga_scores.csv"),
          row.names = FALSE)
write.csv(scores[order(scores$abs_score_diff, decreasing = TRUE), ][seq_len(min(50, nrow(scores))), ],
          chatgpt_result_file("armstrong_pig_mcmc_gaga_top_disagreements.csv"),
          row.names = FALSE)
saveRDS(
  list(config = list(data = data_choice, transform = transform_choice,
                     fdr = fdr, method = method, seed = seed,
                     n_iter = n_iter, burnin = burnin, thin = thin,
                     trunc = trunc, alpha_steps = alpha_steps,
                     alpha_kernel = alpha_kernel,
                     slice_width = slice_width,
                     slice_max_steps = slice_max_steps),
       gaga = gaga_res, pig = pig_res, summary = summary, scores = scores),
  chatgpt_result_file("armstrong_pig_mcmc_gaga.rds")
)

print(summary)
