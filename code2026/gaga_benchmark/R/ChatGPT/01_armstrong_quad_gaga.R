args <- commandArgs(FALSE)
this_file <- sub("^--file=", "", args[grepl("^--file=", args)][1])
source(file.path(dirname(normalizePath(this_file, mustWork = TRUE)), "quad_gaga_lib.R"))

bench <- setup_paths()
require_pkgs(c("gaga", "Biobase"))

fit_gaga_calls <- function(x, groups, patterns, fdr, method) {
  fit <- fit_gaga_original(x, groups, patterns, nclust = 1, method = method, fdr = fdr)
  fg <- gaga::findgenes(fit, x, as.character(groups), fdrmax = fdr, parametric = TRUE)
  list(
    fit = fit,
    calls = as.logical(fg$d),
    score = 1 - fit$pp[, 1],
    threshold = fg$threshold,
    fdr = fg$fdr
  )
}

fit_quad_calls <- function(x, groups, patterns, gaga_fit, fdr, progress_every) {
  quad <- quad_gaga_pp(x, groups, patterns, gaga_fit, progress_every = progress_every)
  list(
    fit = quad,
    calls = find_quad_calls(quad$pp, fdr = fdr),
    score = 1 - quad$pp[, 1]
  )
}

run_one_dataset <- function(x, groups, fdr, method, progress_every) {
  patterns <- gaga_patterns_two_group(groups)
  gg <- fit_gaga_calls(x, groups, patterns, fdr = fdr, method = method)
  quad <- fit_quad_calls(x, groups, patterns, gg$fit, fdr = fdr,
                         progress_every = progress_every)
  names(gg$calls) <- rownames(x)
  names(quad$calls) <- rownames(x)
  names(gg$score) <- rownames(x)
  names(quad$score) <- rownames(x)
  list(gaga = gg, quad = quad, patterns = patterns)
}

summarize_full <- function(res) {
  gaga_calls <- res$gaga$calls
  quad_calls <- res$quad$calls
  gaga_score <- res$gaga$score
  quad_score <- res$quad$score
  overlap <- sum(gaga_calls & quad_calls)
  union <- sum(gaga_calls | quad_calls)
  data.frame(
    method = c("GaGa", "Quad_GaGa"),
    n_de = c(sum(gaga_calls), sum(quad_calls)),
    overlap_with_other = c(overlap, overlap),
    jaccard_with_other = c(safe_div(overlap, union), safe_div(overlap, union)),
    spearman_score_vs_other = c(
      stats::cor(gaga_score, quad_score, method = "spearman"),
      stats::cor(quad_score, gaga_score, method = "spearman")
    ),
    mean_abs_score_diff = c(
      mean(abs(gaga_score - quad_score)),
      mean(abs(quad_score - gaga_score))
    ),
    stringsAsFactors = FALSE
  )
}

run_reproducibility <- function(x, groups, full, reps, seed, fdr, method,
                                sample_sizes, progress_every) {
  if (reps <= 0) {
    return(data.frame())
  }
  set.seed(seed)
  all_idx <- which(groups == "ALL")
  mll_idx <- which(groups == "MLL")
  rows <- list()
  row_id <- 1
  for (rep_id in seq_len(reps)) {
    all_perm <- sample(all_idx)
    mll_perm <- sample(mll_idx)
    for (ss in sample_sizes) {
      if (length(all_perm) < ss || length(mll_perm) < ss) {
        next
      }
      cols <- c(all_perm[seq_len(ss)], mll_perm[seq_len(ss)])
      sub_x <- x[, cols, drop = FALSE]
      sub_groups <- droplevels(groups[cols])
      elapsed <- system.time({
        sub <- run_one_dataset(sub_x, sub_groups, fdr, method, progress_every = 0)
      })[["elapsed"]]
      for (method_name in c("GaGa", "Quad_GaGa")) {
        sub_calls <- if (method_name == "GaGa") sub$gaga$calls else sub$quad$calls
        full_calls <- if (method_name == "GaGa") full$gaga$calls else full$quad$calls
        rows[[row_id]] <- data.frame(
          method = method_name,
          rep = rep_id,
          n_per_group = ss,
          n_de = sum(sub_calls),
          reproducibility = if (sum(sub_calls) == 0) NA_real_ else
            mean(names(sub_calls)[sub_calls] %in% names(full_calls)[full_calls]),
          elapsed_pair_seconds = elapsed,
          stringsAsFactors = FALSE
        )
        row_id <- row_id + 1
      }
    }
  }
  do.call(rbind, rows)
}

data_choice <- Sys.getenv("QUAD_ARM_DATA", unset = "schliep_filtered")
transform_choice <- Sys.getenv("QUAD_ARM_TRANSFORM", unset = "log")
exclude_last_mll <- Sys.getenv("QUAD_ARM_EXCLUDE_LAST_MLL", unset = "1") != "0"
fdr <- env_num("QUAD_ARM_FDR", 0.05)
method <- Sys.getenv("QUAD_ARM_FIT_METHOD", unset = "quickEM")
seed <- env_int("QUAD_ARM_SEED", 2026)
max_genes <- env_int("QUAD_ARM_MAX_GENES", 0)
reps <- env_int("QUAD_ARM_REPS", 2)
sample_sizes <- as.integer(strsplit(Sys.getenv("QUAD_ARM_SAMPLE_SIZES", unset = "5,10,15"), ",")[[1]])
progress_every <- env_int("QUAD_ARM_PROGRESS_EVERY", 500)

data_dir <- file.path(bench, "data")
if (data_choice == "schliep_filtered") {
  data_file <- file.path(data_dir, "armstrong-2002-v2_database.txt")
  dat <- load_armstrong_schliep(data_file, exclude_last_mll = exclude_last_mll)
} else if (data_choice == "orange_full") {
  data_file <- file.path(data_dir, "MLL_orange_12533.tab")
  dat <- load_armstrong_orange(data_file, exclude_last_mll = exclude_last_mll)
} else {
  stop("Unknown QUAD_ARM_DATA: ", data_choice)
}

x_raw <- dat$x
x <- transform_expr(x_raw, transform_choice)
if (max_genes > 0 && max_genes < nrow(x)) {
  vars <- apply(x, 1, stats::var)
  keep <- order(vars, decreasing = TRUE)[seq_len(max_genes)]
  x <- x[keep, , drop = FALSE]
}
groups <- droplevels(dat$groups)

cat("Armstrong Quad-GaGa benchmark\n")
cat("  data=", data_choice, " genes=", nrow(x), " samples=", ncol(x),
    " fdr=", fdr, "\n", sep = "")

elapsed_full <- system.time({
  full <- run_one_dataset(x, groups, fdr, method, progress_every = progress_every)
})[["elapsed"]]

full_summary <- summarize_full(full)
full_summary$data <- data_choice
full_summary$transform <- transform_choice
full_summary$n_genes <- nrow(x)
full_summary$n_samples <- ncol(x)
full_summary$elapsed_full_seconds <- elapsed_full

scores <- data.frame(
  gene = rownames(x),
  gaga_score = full$gaga$score,
  quad_score = full$quad$score,
  gaga_pp_null = full$gaga$fit$pp[, 1],
  quad_pp_null = full$quad$fit$pp[, 1],
  gaga_call = full$gaga$calls,
  quad_call = full$quad$calls,
  stringsAsFactors = FALSE
)
scores$abs_score_diff <- abs(scores$gaga_score - scores$quad_score)

repro <- run_reproducibility(
  x, groups, full, reps = reps, seed = seed, fdr = fdr, method = method,
  sample_sizes = sample_sizes, progress_every = progress_every
)
if (nrow(repro) > 0) {
  repro_summary <- aggregate(
    cbind(n_de, reproducibility, elapsed_pair_seconds) ~ method + n_per_group,
    data = repro,
    FUN = function(v) mean(v, na.rm = TRUE)
  )
} else {
  repro_summary <- data.frame()
}

write.csv(full_summary, chatgpt_result_file("armstrong_quad_gaga_full_summary.csv"),
          row.names = FALSE)
write.csv(scores, chatgpt_result_file("armstrong_quad_gaga_scores.csv"),
          row.names = FALSE)
write.csv(scores[order(scores$abs_score_diff, decreasing = TRUE), ][seq_len(min(50, nrow(scores))), ],
          chatgpt_result_file("armstrong_quad_gaga_top_disagreements.csv"),
          row.names = FALSE)
write.csv(repro, chatgpt_result_file("armstrong_quad_gaga_reproducibility_raw.csv"),
          row.names = FALSE)
write.csv(repro_summary, chatgpt_result_file("armstrong_quad_gaga_reproducibility_summary.csv"),
          row.names = FALSE)
saveRDS(list(config = list(data = data_choice, transform = transform_choice,
                           fdr = fdr, method = method, reps = reps),
             full = full, full_summary = full_summary,
             scores = scores, reproducibility = repro,
             reproducibility_summary = repro_summary),
        chatgpt_result_file("armstrong_quad_gaga.rds"))

print(full_summary)
if (nrow(repro_summary) > 0) {
  print(repro_summary)
}
