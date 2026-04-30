args <- commandArgs(FALSE)
this_file <- sub("^--file=", "", args[grepl("^--file=", args)][1])
src_dir <- dirname(normalizePath(this_file, mustWork = TRUE))
source(file.path(src_dir, "quad_gaga_lib.R"))
source(file.path(src_dir, "pig_mcmc_gaga_lib.R"))

bench <- setup_paths()
require_pkgs(c("gaga", "GIGrvg", "Biobase"))

auc_binary <- function(score, truth) {
  truth <- as.logical(truth)
  if (sum(truth) == 0 || sum(!truth) == 0) {
    return(NA_real_)
  }
  ranks <- rank(score, ties.method = "average")
  n_pos <- sum(truth)
  n_neg <- sum(!truth)
  (sum(ranks[truth]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

eval_calls <- function(calls, truth, score) {
  truth <- as.logical(truth)
  calls <- as.logical(calls)
  data.frame(
    n_de = sum(calls),
    empirical_fdr = if (sum(calls) > 0) mean(!truth[calls]) else 0,
    power = if (sum(truth) > 0) sum(calls & truth) / sum(truth) else NA_real_,
    auc = auc_binary(score, truth),
    stringsAsFactors = FALSE
  )
}

fit_gaga_em <- function(x, groups, patterns, method, fdr) {
  fit <- gaga::fitGG(
    x,
    groups,
    patterns = patterns,
    equalcv = TRUE,
    nclust = 1,
    method = method,
    trace = FALSE
  )
  gaga::parest(fit, x = x, groups = groups, alpha = fdr)
}

find_gaga_calls <- function(fit, x, groups, fdr) {
  as.logical(gaga::findgenes(
    fit,
    x,
    groups,
    fdrmax = fdr,
    parametric = TRUE
  )$d)
}

summarize_pair <- function(gaga_calls, pig_calls, gaga_score, pig_score) {
  overlap <- sum(gaga_calls & pig_calls)
  union <- sum(gaga_calls | pig_calls)
  data.frame(
    overlap = overlap,
    union = union,
    jaccard = safe_div(overlap, union),
    n_gaga_only = sum(gaga_calls & !pig_calls),
    n_pig_only = sum(pig_calls & !gaga_calls),
    spearman_score = stats::cor(
      gaga_score,
      pig_score,
      method = "spearman",
      use = "pairwise.complete.obs"
    ),
    mean_abs_score_diff = mean(abs(gaga_score - pig_score)),
    stringsAsFactors = FALSE
  )
}

simulate_truth <- function(xsim) {
  fd <- Biobase::fData(xsim)
  mean_cols <- grep("^mean\\.", colnames(fd), value = TRUE)
  if (length(mean_cols) >= 2) {
    return(abs(fd[[mean_cols[1]]] - fd[[mean_cols[2]]]) >
             sqrt(.Machine$double.eps))
  }
  stop("Could not infer DE truth from simGG featureData.")
}

n_reps <- env_int("SIM_PIG_REPS", 3)
n_genes <- env_int("SIM_PIG_N_GENES", 200)
n_per_group <- env_int("SIM_PIG_N_PER_GROUP", 10)
p_de <- env_num("SIM_PIG_P_DE", 0.1)
fdr <- env_num("SIM_PIG_FDR", 0.05)
seed <- env_int("SIM_PIG_SEED", 20260430)
method <- Sys.getenv("SIM_PIG_GAGA_METHOD", unset = "EM")

a0 <- env_num("SIM_PIG_A0", 25.5)
nu <- env_num("SIM_PIG_NU", 0.109)
balpha <- env_num("SIM_PIG_BALPHA", 1.183)
nualpha <- env_num("SIM_PIG_NUALPHA", 1683)

n_iter <- env_int("PIG_MCMC_ITER", env_int("SIM_PIG_MCMC_ITER", 250))
burnin <- env_int("PIG_MCMC_BURNIN", env_int("SIM_PIG_MCMC_BURNIN", 100))
thin <- env_int("PIG_MCMC_THIN", env_int("SIM_PIG_MCMC_THIN", 1))
trunc <- env_int("PIG_MCMC_TRUNC", env_int("SIM_PIG_MCMC_TRUNC", 60))
alpha_steps <- env_int("PIG_MCMC_ALPHA_STEPS", env_int("SIM_PIG_ALPHA_STEPS", 1))
alpha_kernel <- Sys.getenv("PIG_MCMC_ALPHA_KERNEL", unset = "slice")
slice_width <- env_num("PIG_MCMC_SLICE_WIDTH", env_num("SIM_PIG_SLICE_WIDTH", 0.5))
slice_max_steps <- env_int("PIG_MCMC_SLICE_MAX_STEPS",
                           env_int("SIM_PIG_SLICE_MAX_STEPS", 50))
progress_every <- env_int("SIM_PIG_PROGRESS_EVERY", 0)

patterns <- matrix(c(0, 0, 0, 1), 2, 2)
colnames(patterns) <- c("group 1", "group 2")

cat("Stirling-vs-PIG simulation\n")
cat("  reps=", n_reps,
    " genes=", n_genes,
    " n_per_group=", n_per_group,
    " p_de=", p_de,
    " gaga_method=", method,
    " pig_iter=", n_iter,
    " pig_burnin=", burnin,
    " pig_trunc=", trunc,
    "\n", sep = "")

method_rows <- list()
pair_rows <- list()
config_rows <- list()

for (rep_id in seq_len(n_reps)) {
  cat("Simulation replicate ", rep_id, "/", n_reps, "\n", sep = "")
  sim_seed <- seed + rep_id
  pig_seed <- seed + 100000L + rep_id
  set.seed(sim_seed)
  xsim <- gaga::simGG(
    n = n_genes,
    m = c(n_per_group, n_per_group),
    p.de = p_de,
    a0 = a0,
    nu = nu,
    balpha = balpha,
    nualpha = nualpha,
    equalcv = TRUE
  )
  x <- Biobase::exprs(xsim)
  groups <- as.character(Biobase::pData(xsim)$group)
  truth <- simulate_truth(xsim)

  elapsed_gaga <- system.time({
    fit <- fit_gaga_em(x, groups, patterns, method = method, fdr = fdr)
    gaga_calls <- find_gaga_calls(fit, x, groups, fdr = fdr)
    gaga_score <- 1 - fit$pp[, 1]
  })[["elapsed"]]

  elapsed_pig <- system.time({
    pig <- pig_mcmc_gaga_pp(
      x,
      groups,
      patterns,
      fit,
      n_iter = n_iter,
      burnin = burnin,
      thin = thin,
      trunc = trunc,
      alpha_steps = alpha_steps,
      alpha_kernel = alpha_kernel,
      slice_width = slice_width,
      slice_max_steps = slice_max_steps,
      progress_every = progress_every,
      seed = pig_seed
    )
    pig_calls <- find_pig_mcmc_calls(pig$pp, fdr = fdr)
    pig_score <- 1 - pig$pp[, 1]
  })[["elapsed"]]

  gaga_eval <- eval_calls(gaga_calls, truth, gaga_score)
  pig_eval <- eval_calls(pig_calls, truth, pig_score)
  method_rows[[length(method_rows) + 1L]] <- data.frame(
    rep = rep_id,
    method = "GaGa_Stirling",
    gaga_eval,
    elapsed_seconds = elapsed_gaga,
    pig_movement_rate_mean = NA_real_,
    pig_movement_rate_median = NA_real_,
    stringsAsFactors = FALSE
  )
  method_rows[[length(method_rows) + 1L]] <- data.frame(
    rep = rep_id,
    method = "GaGa_PIG_MCMC",
    pig_eval,
    elapsed_seconds = elapsed_pig,
    pig_movement_rate_mean = mean(pig$accept_rate, na.rm = TRUE),
    pig_movement_rate_median = stats::median(pig$accept_rate, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  pair_rows[[length(pair_rows) + 1L]] <- data.frame(
    rep = rep_id,
    summarize_pair(gaga_calls, pig_calls, gaga_score, pig_score),
    stringsAsFactors = FALSE
  )
  par <- gaga::getpar(fit)
  config_rows[[length(config_rows) + 1L]] <- data.frame(
    rep = rep_id,
    sim_seed = sim_seed,
    pig_seed = pig_seed,
    truth_de = sum(truth),
    fit_a0 = as.numeric(par$a0[1]),
    fit_nu = as.numeric(par$nu[1]),
    fit_balpha = as.numeric(par$balpha[1]),
    fit_nualpha = as.numeric(par$nualpha[1]),
    fit_probpat0 = as.numeric(par$probpat[1]),
    fit_probpat1 = as.numeric(par$probpat[2]),
    stringsAsFactors = FALSE
  )
}

method_results <- do.call(rbind, method_rows)
pair_results <- do.call(rbind, pair_rows)
config_results <- do.call(rbind, config_rows)

method_summary <- aggregate(
  cbind(n_de, empirical_fdr, power, auc, elapsed_seconds) ~ method,
  method_results,
  function(z) mean(z, na.rm = TRUE)
)
names(method_summary)[-1] <- paste0(names(method_summary)[-1], "_mean")

pair_summary <- data.frame(
  n_reps = n_reps,
  jaccard_mean = mean(pair_results$jaccard, na.rm = TRUE),
  spearman_score_mean = mean(pair_results$spearman_score, na.rm = TRUE),
  mean_abs_score_diff_mean = mean(pair_results$mean_abs_score_diff, na.rm = TRUE),
  n_gaga_only_mean = mean(pair_results$n_gaga_only, na.rm = TRUE),
  n_pig_only_mean = mean(pair_results$n_pig_only, na.rm = TRUE)
)

write.csv(method_results,
          chatgpt_result_file("stirling_vs_pig_mcmc_sim_methods.csv"),
          row.names = FALSE)
write.csv(pair_results,
          chatgpt_result_file("stirling_vs_pig_mcmc_sim_pairwise.csv"),
          row.names = FALSE)
write.csv(config_results,
          chatgpt_result_file("stirling_vs_pig_mcmc_sim_config.csv"),
          row.names = FALSE)
write.csv(method_summary,
          chatgpt_result_file("stirling_vs_pig_mcmc_sim_method_summary.csv"),
          row.names = FALSE)
write.csv(pair_summary,
          chatgpt_result_file("stirling_vs_pig_mcmc_sim_pair_summary.csv"),
          row.names = FALSE)
saveRDS(
  list(
    config = list(
      n_reps = n_reps,
      n_genes = n_genes,
      n_per_group = n_per_group,
      p_de = p_de,
      fdr = fdr,
      seed = seed,
      method = method,
      sim_par = list(a0 = a0, nu = nu, balpha = balpha, nualpha = nualpha),
      pig = list(n_iter = n_iter, burnin = burnin, thin = thin, trunc = trunc,
                 alpha_steps = alpha_steps, alpha_kernel = alpha_kernel,
                 slice_width = slice_width,
                 slice_max_steps = slice_max_steps)
    ),
    methods = method_results,
    pairwise = pair_results,
    fit_config = config_results,
    method_summary = method_summary,
    pair_summary = pair_summary
  ),
  chatgpt_result_file("stirling_vs_pig_mcmc_simulation.rds")
)

print(method_summary)
print(pair_summary)
