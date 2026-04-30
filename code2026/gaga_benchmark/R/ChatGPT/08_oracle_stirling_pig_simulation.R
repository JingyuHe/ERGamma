args <- commandArgs(FALSE)
this_file <- sub("^--file=", "", args[grepl("^--file=", args)][1])
src_dir <- dirname(normalizePath(this_file, mustWork = TRUE))
source(file.path(src_dir, "quad_gaga_lib.R"))
source(file.path(src_dir, "pig_mcmc_gaga_lib.R"))

bench <- setup_paths()
require_pkgs(c("gaga", "GIGrvg"))
if (env_int("ORACLE_SIM_WARN", 0) > 0) {
  options(warn = 1)
}

posterior_fdr_calls <- function(pp, fdr = 0.05) {
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

average_precision <- function(score, truth) {
  truth <- as.logical(truth)
  if (sum(truth) == 0) {
    return(NA_real_)
  }
  ord <- order(score, decreasing = TRUE)
  y <- truth[ord]
  precision <- cumsum(y) / seq_along(y)
  sum(precision[y]) / sum(y)
}

safe_cor <- function(x, y, method) {
  out <- suppressWarnings(stats::cor(x, y, method = method,
                                     use = "pairwise.complete.obs"))
  ifelse(is.finite(out), out, NA_real_)
}

eval_calls <- function(pp, truth, fdr) {
  score <- 1 - pp[, 1]
  calls <- posterior_fdr_calls(pp, fdr)
  data.frame(
    fdr_target = fdr,
    n_selected = sum(calls),
    empirical_fdr = if (sum(calls) > 0) mean(!truth[calls]) else 0,
    power = if (sum(truth) > 0) mean(calls[truth]) else NA_real_,
    auc = auc_binary(score, truth),
    avg_precision = average_precision(score, truth),
    stringsAsFactors = FALSE
  )
}

score_error <- function(pp, pp_exact) {
  d0 <- pp[, 1] - pp_exact[, 1]
  s <- (1 - pp[, 1]) - (1 - pp_exact[, 1])
  data.frame(
    mae_pp0 = mean(abs(d0), na.rm = TRUE),
    rmse_pp0 = sqrt(mean(d0^2, na.rm = TRUE)),
    q95_abs_pp0 = as.numeric(stats::quantile(abs(d0), 0.95, na.rm = TRUE)),
    max_abs_pp0 = max(abs(d0), na.rm = TRUE),
    mae_score = mean(abs(s), na.rm = TRUE),
    spearman_score = safe_cor(1 - pp[, 1], 1 - pp_exact[, 1], "spearman"),
    pearson_score = safe_cor(1 - pp[, 1], 1 - pp_exact[, 1], "pearson"),
    stringsAsFactors = FALSE
  )
}

jaccard_calls <- function(pp, pp_exact, fdr) {
  calls <- posterior_fdr_calls(pp, fdr)
  calls_exact <- posterior_fdr_calls(pp_exact, fdr)
  overlap <- sum(calls & calls_exact)
  union <- sum(calls | calls_exact)
  data.frame(
    fdr_target = fdr,
    overlap = overlap,
    union = union,
    jaccard_vs_exact = safe_div(overlap, union),
    exact_only = sum(calls_exact & !calls),
    method_only = sum(calls & !calls_exact),
    stringsAsFactors = FALSE
  )
}

make_oracle_fit <- function(patterns, par) {
  fit <- list(
    parest = c(
      a0 = par$a0,
      nu = par$nu,
      balpha = par$balpha,
      nualpha = par$nualpha,
      probclus = 1,
      probpat = par$probpat
    ),
    mcmc = coda::as.mcmc(NA),
    lhood = NA_real_,
    equalcv = TRUE,
    nclust = 1L,
    patterns = patterns,
    method = "EM"
  )
  class(fit) <- "gagafit"
  fit
}

simulate_gaga_two_group <- function(n_genes, n_per_group, par, seed) {
  set.seed(seed)
  groups <- rep(c("G1", "G2"), each = n_per_group)
  x <- matrix(NA_real_, n_genes, length(groups))
  rownames(x) <- paste0("gene", seq_len(n_genes))
  colnames(x) <- paste0(groups, "_", ave(seq_along(groups), groups,
                                          FUN = seq_along))

  z <- rbinom(n_genes, size = 1, prob = par$probpat[2])
  alpha <- stats::rgamma(n_genes, shape = par$balpha,
                         rate = par$balpha / par$nualpha)
  alpha <- pmax(alpha, .Machine$double.eps)
  lambda <- matrix(NA_real_, n_genes, 2)

  for (i in seq_len(n_genes)) {
    if (z[i] == 0L) {
      q <- stats::rgamma(1, shape = par$a0, rate = par$a0 / par$nu)
      lambda[i, ] <- 1 / q
    } else {
      q <- stats::rgamma(2, shape = par$a0, rate = par$a0 / par$nu)
      lambda[i, ] <- 1 / q
    }
    x[i, groups == "G1"] <- stats::rgamma(
      n_per_group,
      shape = alpha[i],
      rate = alpha[i] / lambda[i, 1]
    )
    x[i, groups == "G2"] <- stats::rgamma(
      n_per_group,
      shape = alpha[i],
      rate = alpha[i] / lambda[i, 2]
    )
  }

  list(
    x = pmax(x, .Machine$double.xmin),
    groups = groups,
    truth_de = z == 1L,
    z = z,
    alpha = alpha,
    lambda = lambda
  )
}

log_integral_alpha_oracle <- function(n_c, S_c, L_c, hyper,
                                      lower, upper, rel_tol) {
  logf <- function(t) {
    vapply(t, function(tt) {
      alpha <- exp(tt)
      pig_log_alpha_target(alpha, n_c, S_c, L_c, hyper) + tt
    }, numeric(1))
  }
  opt <- tryCatch(
    stats::optimize(function(t) -logf(t), interval = c(lower, upper)),
    error = function(e) NULL
  )
  if (is.null(opt)) {
    return(NA_real_)
  }
  mode_t <- opt$minimum
  log_scale <- logf(mode_t)
  if (!is.finite(log_scale)) {
    return(NA_real_)
  }
  value <- tryCatch({
    stats::integrate(
      function(t) exp(logf(t) - log_scale),
      lower,
      upper,
      subdivisions = 300,
      rel.tol = rel_tol,
      stop.on.error = FALSE
    )$value
  }, error = function(e) NA_real_)
  if (!is.finite(value) || value <= 0) {
    return(NA_real_)
  }
  log(value) + log_scale
}

exact_oracle_pp <- function(x, groups, patterns, par,
                            lower = -12, upper = NULL,
                            rel_tol = 1e-5, progress_every = 0) {
  hyper <- list(
    alpha0 = par$a0,
    nu = par$nu,
    b_alpha = par$balpha,
    nu_alpha = par$nualpha
  )
  if (is.null(upper)) {
    alpha_upper <- stats::qgamma(
      1 - 1e-10,
      shape = par$balpha,
      rate = par$balpha / par$nualpha
    )
    upper <- log(max(5000, alpha_upper, 10 * par$nualpha))
  }
  probpat <- pmax(par$probpat, .Machine$double.eps)
  probpat <- probpat / sum(probpat)
  st_all <- pig_pattern_stats(x, groups, patterns)
  pp <- matrix(NA_real_, nrow = nrow(x), ncol = nrow(patterns))
  log_m <- pp
  rownames(pp) <- rownames(x)
  colnames(pp) <- paste0("pattern", seq_len(nrow(patterns)) - 1)
  rownames(log_m) <- rownames(x)
  colnames(log_m) <- colnames(pp)

  for (i in seq_len(nrow(x))) {
    if (progress_every > 0 && (i == 1 || i %% progress_every == 0)) {
      cat("  Exact quadrature gene ", i, "/", nrow(x), "\n", sep = "")
    }
    logp <- numeric(nrow(patterns))
    for (h in seq_len(nrow(patterns))) {
      st <- st_all[[h]]
      li <- log_integral_alpha_oracle(
        st$n_c,
        st$sums[i, ],
        st$log_sums[i, ],
        hyper,
        lower = lower,
        upper = upper,
        rel_tol = rel_tol
      )
      log_m[i, h] <- li
      logp[h] <- log(probpat[h]) + li
    }
    pp[i, ] <- softmax_log(logp)
  }
  list(pp = pp, log_marginal = log_m, lower = lower, upper = upper)
}

stirling_oracle_pp <- function(x, groups, patterns, par) {
  out <- getFromNamespace("ppGG", "gaga")(
    x,
    groups,
    a0 = par$a0,
    nu = par$nu,
    balpha = par$balpha,
    nualpha = par$nualpha,
    equalcv = TRUE,
    probclus = 1,
    probpat = par$probpat,
    patterns = patterns
  )$pp
  rownames(out) <- rownames(x)
  colnames(out) <- paste0("pattern", seq_len(ncol(out)) - 1)
  out
}

alpha_bin_errors <- function(pp, pp_exact, alpha, scenario, rep_id, method) {
  bins <- cut(
    alpha,
    breaks = c(-Inf, 1, 3, 10, Inf),
    labels = c("<1", "1-3", "3-10", ">=10"),
    right = FALSE
  )
  d <- abs(pp[, 1] - pp_exact[, 1])
  rows <- lapply(levels(bins), function(bin) {
    idx <- bins == bin
    data.frame(
      scenario = scenario,
      rep = rep_id,
      method = method,
      alpha_bin = bin,
      n = sum(idx, na.rm = TRUE),
      mae_pp0 = if (any(idx, na.rm = TRUE)) mean(d[idx], na.rm = TRUE) else NA_real_,
      q95_abs_pp0 = if (any(idx, na.rm = TRUE)) {
        as.numeric(stats::quantile(d[idx], 0.95, na.rm = TRUE))
      } else {
        NA_real_
      },
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

scenario_table <- function() {
  data.frame(
    scenario = c(
      "S0_paper_like",
      "S1_benign",
      "S2_moderate_smalln",
      "S3_low_alpha",
      "S4_severe_stress",
      "S5_weak_signal"
    ),
    a0 = c(25.5, 15, 8, 5, 1.5, 50),
    nu = c(0.109, 1, 1, 1, 1, 1),
    balpha = c(1.183, 8, 3, 1.5, 0.8, 3),
    nualpha = c(1683, 50, 8, 1.5, 0.7, 5),
    n_per_group = c(10L, 10L, 5L, 5L, 3L, 5L),
    p_de = c(0.10, 0.10, 0.10, 0.10, 0.10, 0.05),
    stringsAsFactors = FALSE
  )
}

parse_scenarios <- function(scenarios, selected) {
  if (!nzchar(selected) || selected == "all") {
    return(scenarios)
  }
  keep <- trimws(strsplit(selected, ",", fixed = TRUE)[[1]])
  out <- scenarios[scenarios$scenario %in% keep, , drop = FALSE]
  if (nrow(out) == 0) {
    stop("No scenarios matched SIM_SCENARIOS=", selected)
  }
  out
}

bind_rows <- function(rows) {
  if (!length(rows)) {
    return(data.frame())
  }
  do.call(rbind, rows)
}

paper_summary <- function(data, by, vars) {
  key <- interaction(data[by], drop = TRUE, lex.order = TRUE, sep = "\r")
  pieces <- lapply(vars, function(v) {
    rows <- lapply(split(seq_len(nrow(data)), key), function(idx) {
      z <- data[[v]][idx]
      z <- z[is.finite(z)]
      data.frame(
        data[idx[1], by, drop = FALSE],
        n = length(z),
        mean = if (length(z)) mean(z) else NA_real_,
        sd = if (length(z) > 1) stats::sd(z) else NA_real_,
        se = if (length(z) > 1) stats::sd(z) / sqrt(length(z)) else NA_real_,
        row.names = NULL
      )
    })
    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    data.frame(out[, by, drop = FALSE], metric = v,
               out[, c("n", "mean", "sd", "se"), drop = FALSE],
               stringsAsFactors = FALSE)
  })
  do.call(rbind, pieces)
}

n_reps <- env_int("ORACLE_SIM_REPS", 20)
n_genes <- env_int("ORACLE_SIM_GENES", 500)
seed <- env_int("ORACLE_SIM_SEED", 20260430)
fdr_grid <- as.numeric(strsplit(Sys.getenv("ORACLE_SIM_FDR", "0.01,0.05,0.10"),
                                ",", fixed = TRUE)[[1]])
exact_rel_tol <- env_num("ORACLE_SIM_EXACT_RELTOL", 1e-5)
exact_lower <- env_num("ORACLE_SIM_EXACT_LOWER", -12)
exact_progress <- env_int("ORACLE_SIM_EXACT_PROGRESS", 0)
selected_scenarios <- Sys.getenv("ORACLE_SIM_SCENARIOS", unset = "all")

n_iter <- env_int("PIG_MCMC_ITER", env_int("ORACLE_SIM_PIG_ITER", 800))
burnin <- env_int("PIG_MCMC_BURNIN", env_int("ORACLE_SIM_PIG_BURNIN", 300))
thin <- env_int("PIG_MCMC_THIN", env_int("ORACLE_SIM_PIG_THIN", 2))
trunc <- env_int("PIG_MCMC_TRUNC", env_int("ORACLE_SIM_PIG_TRUNC", 60))
alpha_steps <- env_int("PIG_MCMC_ALPHA_STEPS",
                       env_int("ORACLE_SIM_PIG_ALPHA_STEPS", 1))
alpha_kernel <- Sys.getenv("PIG_MCMC_ALPHA_KERNEL", unset = "slice")
slice_width <- env_num("PIG_MCMC_SLICE_WIDTH",
                       env_num("ORACLE_SIM_PIG_SLICE_WIDTH", 0.5))
slice_max_steps <- env_int("PIG_MCMC_SLICE_MAX_STEPS",
                           env_int("ORACLE_SIM_PIG_SLICE_MAX_STEPS", 50))
pig_progress <- env_int("ORACLE_SIM_PIG_PROGRESS", 0)

scenarios <- parse_scenarios(scenario_table(), selected_scenarios)
patterns <- matrix(c(0, 0, 0, 1), 2, 2, byrow = TRUE)
colnames(patterns) <- c("G1", "G2")

cat("Oracle GaGa Stirling vs PIG-MCMC simulation\n")
cat("  scenarios=", paste(scenarios$scenario, collapse = ","),
    " reps=", n_reps,
    " genes=", n_genes,
    " pig_iter=", n_iter,
    " pig_burnin=", burnin,
    " pig_thin=", thin,
    " pig_trunc=", trunc,
    "\n", sep = "")

approx_rows <- list()
decision_rows <- list()
jaccard_rows <- list()
alpha_rows <- list()
config_rows <- list()
timing_rows <- list()

for (sidx in seq_len(nrow(scenarios))) {
  sc <- scenarios[sidx, ]
  par <- list(
    a0 = sc$a0,
    nu = sc$nu,
    balpha = sc$balpha,
    nualpha = sc$nualpha,
    probpat = c(1 - sc$p_de, sc$p_de)
  )
  fit <- make_oracle_fit(patterns, par)
  for (rep_id in seq_len(n_reps)) {
    rep_seed <- seed + 10000L * sidx + rep_id
    pig_seed <- seed + 900000L + 10000L * sidx + rep_id
    cat("Scenario ", sc$scenario, " replicate ", rep_id, "/", n_reps, "\n",
        sep = "")

    sim <- simulate_gaga_two_group(
      n_genes = n_genes,
      n_per_group = sc$n_per_group,
      par = par,
      seed = rep_seed
    )

    t_exact <- system.time({
      exact <- exact_oracle_pp(
        sim$x,
        sim$groups,
        patterns,
        par,
        lower = exact_lower,
        rel_tol = exact_rel_tol,
        progress_every = exact_progress
      )
    })[["elapsed"]]

    t_stirling <- system.time({
      pp_stirling <- stirling_oracle_pp(sim$x, sim$groups, patterns, par)
    })[["elapsed"]]

    t_pig <- system.time({
      pig <- pig_mcmc_gaga_pp(
        sim$x,
        sim$groups,
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
        progress_every = pig_progress,
        seed = pig_seed
      )
    })[["elapsed"]]

    methods <- list(
      Exact = exact$pp,
      Stirling = pp_stirling,
      PIG_MCMC = pig$pp
    )
    for (method_name in c("Stirling", "PIG_MCMC")) {
      approx_rows[[length(approx_rows) + 1L]] <- data.frame(
        scenario = sc$scenario,
        rep = rep_id,
        method = method_name,
        score_error(methods[[method_name]], exact$pp),
        stringsAsFactors = FALSE
      )
      alpha_rows[[length(alpha_rows) + 1L]] <- alpha_bin_errors(
        methods[[method_name]],
        exact$pp,
        sim$alpha,
        scenario = sc$scenario,
        rep_id = rep_id,
        method = method_name
      )
    }

    for (method_name in names(methods)) {
      for (fdr in fdr_grid) {
        decision_rows[[length(decision_rows) + 1L]] <- data.frame(
          scenario = sc$scenario,
          rep = rep_id,
          method = method_name,
          eval_calls(methods[[method_name]], sim$truth_de, fdr),
          stringsAsFactors = FALSE
        )
        if (method_name != "Exact") {
          jaccard_rows[[length(jaccard_rows) + 1L]] <- data.frame(
            scenario = sc$scenario,
            rep = rep_id,
            method = method_name,
            jaccard_calls(methods[[method_name]], exact$pp, fdr),
            stringsAsFactors = FALSE
          )
        }
      }
    }

    config_rows[[length(config_rows) + 1L]] <- data.frame(
      scenario = sc$scenario,
      rep = rep_id,
      seed = rep_seed,
      pig_seed = pig_seed,
      n_genes = n_genes,
      n_per_group = sc$n_per_group,
      truth_de = sum(sim$truth_de),
      a0 = sc$a0,
      nu = sc$nu,
      balpha = sc$balpha,
      nualpha = sc$nualpha,
      p_de = sc$p_de,
      exact_lower = exact$lower,
      exact_upper = exact$upper,
      pig_movement_rate_mean = mean(pig$accept_rate, na.rm = TRUE),
      pig_movement_rate_median = stats::median(pig$accept_rate, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    timing_rows[[length(timing_rows) + 1L]] <- data.frame(
      scenario = sc$scenario,
      rep = rep_id,
      exact_seconds = t_exact,
      stirling_seconds = t_stirling,
      pig_seconds = t_pig,
      stringsAsFactors = FALSE
    )
  }
}

approx_results <- bind_rows(approx_rows)
decision_results <- bind_rows(decision_rows)
jaccard_results <- bind_rows(jaccard_rows)
alpha_results <- bind_rows(alpha_rows)
config_results <- bind_rows(config_rows)
timing_results <- bind_rows(timing_rows)

approx_summary <- stats::aggregate(
  cbind(mae_pp0, rmse_pp0, q95_abs_pp0, max_abs_pp0, mae_score,
        spearman_score, pearson_score) ~ scenario + method,
  approx_results,
  function(z) mean(z, na.rm = TRUE)
)
names(approx_summary)[-(1:2)] <- paste0(names(approx_summary)[-(1:2)], "_mean")

decision_summary <- stats::aggregate(
  cbind(n_selected, empirical_fdr, power, auc, avg_precision) ~
    scenario + method + fdr_target,
  decision_results,
  function(z) mean(z, na.rm = TRUE)
)
names(decision_summary)[-(1:3)] <- paste0(names(decision_summary)[-(1:3)],
                                          "_mean")

jaccard_summary <- stats::aggregate(
  cbind(jaccard_vs_exact, exact_only, method_only) ~
    scenario + method + fdr_target,
  jaccard_results,
  function(z) mean(z, na.rm = TRUE)
)
names(jaccard_summary)[-(1:3)] <- paste0(names(jaccard_summary)[-(1:3)],
                                         "_mean")

alpha_summary <- stats::aggregate(
  cbind(n, mae_pp0, q95_abs_pp0) ~ scenario + method + alpha_bin,
  alpha_results,
  function(z) mean(z, na.rm = TRUE)
)
names(alpha_summary)[-(1:3)] <- paste0(names(alpha_summary)[-(1:3)], "_mean")

timing_summary <- stats::aggregate(
  cbind(exact_seconds, stirling_seconds, pig_seconds) ~ scenario,
  timing_results,
  function(z) mean(z, na.rm = TRUE)
)
names(timing_summary)[-1] <- paste0(names(timing_summary)[-1], "_mean")

prefix <- Sys.getenv("ORACLE_SIM_PREFIX", unset = "oracle_stirling_pig")
out <- function(suffix) chatgpt_result_file(paste0(prefix, "_", suffix))

write.csv(approx_results, out("approx_by_rep.csv"), row.names = FALSE)
write.csv(decision_results, out("decision_by_rep.csv"), row.names = FALSE)
write.csv(jaccard_results, out("jaccard_by_rep.csv"), row.names = FALSE)
write.csv(alpha_results, out("alpha_bin_by_rep.csv"), row.names = FALSE)
write.csv(config_results, out("config_by_rep.csv"), row.names = FALSE)
write.csv(timing_results, out("timing_by_rep.csv"), row.names = FALSE)
write.csv(approx_summary, out("approx_summary.csv"), row.names = FALSE)
write.csv(decision_summary, out("decision_summary.csv"), row.names = FALSE)
write.csv(jaccard_summary, out("jaccard_summary.csv"), row.names = FALSE)
write.csv(alpha_summary, out("alpha_bin_summary.csv"), row.names = FALSE)
write.csv(timing_summary, out("timing_summary.csv"), row.names = FALSE)
write.csv(
  paper_summary(
    approx_results,
    by = c("scenario", "method"),
    vars = c("mae_pp0", "rmse_pp0", "q95_abs_pp0", "max_abs_pp0",
             "spearman_score", "pearson_score")
  ),
  out("approx_paper_summary.csv"),
  row.names = FALSE
)
write.csv(
  paper_summary(
    decision_results,
    by = c("scenario", "method", "fdr_target"),
    vars = c("n_selected", "empirical_fdr", "power", "auc",
             "avg_precision")
  ),
  out("decision_paper_summary.csv"),
  row.names = FALSE
)
write.csv(
  paper_summary(
    jaccard_results,
    by = c("scenario", "method", "fdr_target"),
    vars = c("jaccard_vs_exact", "exact_only", "method_only")
  ),
  out("jaccard_paper_summary.csv"),
  row.names = FALSE
)
write.csv(
  paper_summary(
    alpha_results,
    by = c("scenario", "method", "alpha_bin"),
    vars = c("mae_pp0", "q95_abs_pp0")
  ),
  out("alpha_bin_paper_summary.csv"),
  row.names = FALSE
)
saveRDS(
  list(
    scenarios = scenarios,
    config = list(
      n_reps = n_reps,
      n_genes = n_genes,
      seed = seed,
      fdr_grid = fdr_grid,
      exact_rel_tol = exact_rel_tol,
      exact_lower = exact_lower,
      pig = list(n_iter = n_iter, burnin = burnin, thin = thin,
                 trunc = trunc, alpha_steps = alpha_steps,
                 alpha_kernel = alpha_kernel,
                 slice_width = slice_width,
                 slice_max_steps = slice_max_steps)
    ),
    approx = approx_results,
    decision = decision_results,
    jaccard = jaccard_results,
    alpha_bin = alpha_results,
    timing = timing_results,
    summaries = list(
      approx = approx_summary,
      decision = decision_summary,
      jaccard = jaccard_summary,
      alpha_bin = alpha_summary,
      timing = timing_summary
    )
  ),
  out("full_results.rds")
)

pdf(out("summary_plots.pdf"), width = 10, height = 7)
op <- par(no.readonly = TRUE)
on.exit(par(op), add = TRUE)
par(mfrow = c(2, 2), mar = c(8, 4, 3, 1))

for (metric in c("mae_pp0_mean", "q95_abs_pp0_mean")) {
  mat <- xtabs(as.formula(paste(metric, "~ method + scenario")),
               data = approx_summary)
  barplot(mat, beside = TRUE, las = 2, ylab = metric,
          main = paste("Approximation error:", metric),
          col = c("gray55", "steelblue"))
  legend("topright", legend = rownames(mat), fill = c("gray55", "steelblue"),
         bty = "n")
}

dec05 <- decision_summary[abs(decision_summary$fdr_target - 0.05) < 1e-12, ]
mat_power <- xtabs(power_mean ~ method + scenario, data = dec05)
barplot(mat_power, beside = TRUE, las = 2, ylab = "mean power",
        main = "Power at posterior FDR 0.05",
        col = c("black", "gray55", "steelblue"))
legend("topright", legend = rownames(mat_power),
       fill = c("black", "gray55", "steelblue"), bty = "n")

mat_fdr <- xtabs(empirical_fdr_mean ~ method + scenario, data = dec05)
barplot(mat_fdr, beside = TRUE, las = 2, ylab = "mean empirical FDR",
        main = "Empirical FDR at posterior FDR 0.05",
        col = c("black", "gray55", "steelblue"))
abline(h = 0.05, col = "red", lty = 2)
legend("topright", legend = rownames(mat_fdr),
       fill = c("black", "gray55", "steelblue"), bty = "n")
dev.off()

cat("\nApproximation summary:\n")
print(approx_summary)
cat("\nDecision summary:\n")
print(decision_summary)
cat("\nAlpha-bin summary:\n")
print(alpha_summary)
cat("\nTiming summary:\n")
print(timing_summary)
