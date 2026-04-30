script_path <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) == 0) {
    return(normalizePath("code2026/gamma_shape_inference/run_pig_exactness_chain_length.R"))
  }
  normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)
}

script_dir <- dirname(script_path())
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)
source(file.path(script_dir, "gamma_shape_lib.R"))
gamma_shape_init_libpaths(repo_root)

arg_value <- function(flag, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  hit <- which(args == flag)
  if (length(hit) > 0 && hit[1] < length(args)) {
    return(args[hit[1] + 1])
  }
  prefix <- paste0(flag, "=")
  hit2 <- grep(paste0("^", prefix), args, value = TRUE)
  if (length(hit2) > 0) {
    return(sub(paste0("^", prefix), "", hit2[[1]]))
  }
  default
}

mc_se <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) <= 1) return(NA_real_)
  stats::sd(x) / sqrt(length(x))
}

seed <- as.integer(arg_value("--seed", "20260430"))
reps <- as.integer(arg_value("--reps", "4"))
cores <- as.integer(arg_value("--cores", "4"))
chunk_size <- as.integer(arg_value("--chunk-size", as.character(max(cores, 1))))
grid_size <- as.integer(arg_value("--grid-size", "20001"))
tail_nats <- as.numeric(arg_value("--tail-nats", "50"))
pig_N <- as.integer(arg_value("--pig-N", "200"))
n_iter_values <- as.integer(strsplit(arg_value("--iter-values", "5000,20000,100000"), ",")[[1]])
burnin_frac <- as.numeric(arg_value("--burnin-frac", "0.2"))
out_dir <- arg_value(
  "--out-dir",
  file.path(repo_root, "code2026", "gamma_shape_inference", "results", "pig_exactness_chain_length")
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cells <- data.frame(
  cell_id = seq_len(3),
  cell = c("low-shape", "moderate-shape", "high-shape-large-n"),
  alpha_true = c(0.1, 5, 20),
  n = c(100, 50, 500),
  stringsAsFactors = FALSE
)

make_cell_target <- function(cell_row) {
  set.seed(seed + 10000L * cell_row$cell_id)
  stats <- simulate_gamma_stats(cell_row$alpha_true, cell_row$n)
  pp <- posterior_params_from_stats(
    n = cell_row$n,
    x_a = stats$x_a,
    x_g = stats$x_g,
    prior_delta = 0,
    prior_center = NA_real_
  )
  truth <- xi2_grid_truth(
    delta = pp$delta_post,
    mu_ratio = pp$mu_ratio,
    grid_size = grid_size,
    tail_nats = tail_nats
  )
  list(stats = stats, pp = pp, truth = truth)
}

targets <- lapply(split(cells, cells$cell_id), make_cell_target)

approx_rows <- list()
for (i in seq_len(nrow(cells))) {
  cell_row <- cells[i, ]
  target <- targets[[as.character(cell_row$cell_id)]]
  truth <- target$truth
  pp <- target$pp

  mp <- miller_gamma_params(pp$delta_post, pp$mu_ratio, init = truth$mode)
  approx_rows[[length(approx_rows) + 1]] <- cbind(
    cell_row,
    n_iter = NA_integer_,
    burnin = NA_integer_,
    rep = NA_integer_,
    x_a = target$stats$x_a,
    x_g = target$stats$x_g,
    delta_post = pp$delta_post,
    mu_ratio = pp$mu_ratio,
    gamma_approx_row(
      truth = truth,
      shape = mp$shape,
      rate = mp$rate,
      method = "Miller-Gamma",
      runtime_sec = NA_real_,
      alpha_true = cell_row$alpha_true,
      status = if (isTRUE(mp$ok)) "ok" else "failed",
      message = if (isTRUE(mp$ok)) "" else "Miller derivative matching failed"
    )
  )

  sp <- stirling_gamma_params(pp$delta_post, pp$mu_ratio)
  approx_rows[[length(approx_rows) + 1]] <- cbind(
    cell_row,
    n_iter = NA_integer_,
    burnin = NA_integer_,
    rep = NA_integer_,
    x_a = target$stats$x_a,
    x_g = target$stats$x_g,
    delta_post = pp$delta_post,
    mu_ratio = pp$mu_ratio,
    gamma_approx_row(
      truth = truth,
      shape = sp$shape,
      rate = sp$rate,
      method = "Stirling-Gamma",
      runtime_sec = NA_real_,
      alpha_true = cell_row$alpha_true,
      status = if (isTRUE(sp$ok)) "ok" else "failed",
      message = if (isTRUE(sp$ok)) "" else "Stirling approximation failed"
    )
  )
}
approx_raw <- do.call(rbind, approx_rows)

tasks <- expand.grid(
  cell_id = cells$cell_id,
  n_iter = n_iter_values,
  rep = seq_len(reps),
  KEEP.OUT.ATTRS = FALSE
)
tasks <- merge(tasks, cells, by = "cell_id", sort = FALSE)
tasks$task_id <- seq_len(nrow(tasks))
tasks <- tasks[order(tasks$task_id), ]

run_one <- function(task) {
  cell_row <- task[1, c("cell_id", "cell", "alpha_true", "n")]
  target <- targets[[as.character(task$cell_id)]]
  pp <- target$pp
  truth <- target$truth
  burnin <- as.integer(round(task$n_iter * burnin_frac))
  set.seed(seed + 100000L * task$cell_id + 1000L * match(task$n_iter, n_iter_values) + task$rep)

  t0 <- proc.time()[["elapsed"]]
  pig <- sample_pig_xi2(
    delta = pp$delta_post,
    mu_ratio = pp$mu_ratio,
    n_iter = task$n_iter,
    burnin = burnin,
    thin = 1,
    pig_N = pig_N,
    init = truth$mode,
    max_delta = Inf
  )
  runtime <- proc.time()[["elapsed"]] - t0
  cbind(
    cell_row,
    n_iter = task$n_iter,
    burnin = burnin,
    rep = task$rep,
    x_a = target$stats$x_a,
    x_g = target$stats$x_g,
    delta_post = pp$delta_post,
    mu_ratio = pp$mu_ratio,
    sample_summary_row(
      samples = pig$samples,
      truth = truth,
      method = "P-IG Gibbs",
      runtime_sec = runtime,
      n_iter = task$n_iter,
      alpha_true = task$alpha_true,
      pig_N = pig_N
    )
  )
}

message("P-IG exactness chain-length experiment")
message("  seed: ", seed)
message("  reps: ", reps)
message("  iter values: ", paste(n_iter_values, collapse = ", "))
message("  output: ", out_dir)
message("  cores/chunk_size: ", cores, "/", chunk_size)

start_time <- Sys.time()
raw_list <- list()
task_chunks <- split(seq_len(nrow(tasks)), ceiling(seq_len(nrow(tasks)) / chunk_size))
for (chunk_index in seq_along(task_chunks)) {
  idx <- task_chunks[[chunk_index]]
  chunk_tasks <- tasks[idx, , drop = FALSE]
  t0 <- proc.time()[["elapsed"]]
  chunk_list <- if (cores > 1 && .Platform$OS.type != "windows") {
    parallel::mclapply(
      split(chunk_tasks, seq_len(nrow(chunk_tasks))),
      run_one,
      mc.cores = cores,
      mc.set.seed = TRUE
    )
  } else {
    lapply(split(chunk_tasks, seq_len(nrow(chunk_tasks))), run_one)
  }
  raw_list <- c(raw_list, chunk_list)
  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  completed <- length(raw_list)
  rate <- completed / elapsed
  eta <- if (rate > 0) (nrow(tasks) - completed) / rate else NA_real_
  message(sprintf(
    "[%s] completed %d/%d (%.1f%%); chunk %.1fs; elapsed %.1f min; ETA %.1f min",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    completed,
    nrow(tasks),
    100 * completed / nrow(tasks),
    proc.time()[["elapsed"]] - t0,
    elapsed / 60,
    eta / 60
  ))
}
pig_raw <- do.call(rbind, raw_list)
rownames(pig_raw) <- NULL
raw <- rbind(approx_raw, pig_raw)

summary_groups <- split(
  raw[raw$status == "ok", , drop = FALSE],
  paste(raw$cell[raw$status == "ok"], raw$method[raw$status == "ok"], raw$n_iter[raw$status == "ok"], sep = "\r")
)
summary <- do.call(rbind, lapply(summary_groups, function(d) {
  data.frame(
    cell = d$cell[1],
    alpha_true = d$alpha_true[1],
    n = d$n[1],
    method = d$method[1],
    n_iter = d$n_iter[1],
    burnin = d$burnin[1],
    reps = nrow(d),
    ks_mean = mean(d$ks, na.rm = TRUE),
    ks_se = mc_se(d$ks),
    w1_mean = mean(d$w1, na.rm = TRUE),
    w1_se = mc_se(d$w1),
    coverage = mean(d$cover, na.rm = TRUE),
    ess_iter_median = stats::median(d$ess_per_iter, na.rm = TRUE),
    ess_sec_median = stats::median(d$ess_per_sec, na.rm = TRUE),
    runtime_median = stats::median(d$runtime_sec, na.rm = TRUE),
    runtime_p90 = as.numeric(stats::quantile(d$runtime_sec, 0.9, na.rm = TRUE, names = FALSE)),
    acf1_mean = mean(d$acf1, na.rm = TRUE),
    acf10_mean = mean(d$acf10, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))
rownames(summary) <- NULL
summary <- summary[order(summary$cell, summary$method, summary$n_iter), ]

raw_file <- file.path(out_dir, "pig_exactness_chain_length_raw.csv")
summary_file <- file.path(out_dir, "pig_exactness_chain_length_summary.csv")
rds_file <- file.path(out_dir, "pig_exactness_chain_length_results.rds")
write.csv(raw, raw_file, row.names = FALSE)
write.csv(summary, summary_file, row.names = FALSE)
saveRDS(
  list(
    raw = raw,
    summary = summary,
    tasks = tasks,
    cells = cells,
    settings = list(
      seed = seed,
      reps = reps,
      cores = cores,
      grid_size = grid_size,
      tail_nats = tail_nats,
      pig_N = pig_N,
      n_iter_values = n_iter_values,
      burnin_frac = burnin_frac
    ),
    elapsed_sec = as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  ),
  rds_file
)

if (requireNamespace("ggplot2", quietly = TRUE)) {
  fig_dir <- file.path(out_dir, "figures")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  pdat <- summary[summary$method == "P-IG Gibbs", , drop = FALSE]
  adat <- summary[summary$method %in% c("Miller-Gamma", "Stirling-Gamma"), , drop = FALSE]

  p_ks <- ggplot2::ggplot(pdat, ggplot2::aes(n_iter, ks_mean)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = pmax(0, ks_mean - ks_se), ymax = ks_mean + ks_se),
      fill = "#9ecae1",
      alpha = 0.35
    ) +
    ggplot2::geom_line(colour = "#08519c", linewidth = 0.7) +
    ggplot2::geom_point(colour = "#08519c", size = 1.8) +
    ggplot2::geom_hline(
      data = adat,
      ggplot2::aes(yintercept = ks_mean, colour = method, linetype = method),
      linewidth = 0.6
    ) +
    ggplot2::facet_wrap(~ cell, scales = "free_y") +
    ggplot2::scale_x_log10(breaks = n_iter_values) +
    ggplot2::labs(x = "P-IG iterations", y = "KS distance to grid", colour = NULL, linetype = NULL) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
  ggplot2::ggsave(file.path(fig_dir, "pig_exactness_ks_vs_iterations.pdf"), p_ks, width = 7.5, height = 4)

  p_w1 <- ggplot2::ggplot(pdat, ggplot2::aes(n_iter, w1_mean)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = pmax(0, w1_mean - w1_se), ymax = w1_mean + w1_se),
      fill = "#a1d99b",
      alpha = 0.35
    ) +
    ggplot2::geom_line(colour = "#006d2c", linewidth = 0.7) +
    ggplot2::geom_point(colour = "#006d2c", size = 1.8) +
    ggplot2::geom_hline(
      data = adat,
      ggplot2::aes(yintercept = w1_mean, colour = method, linetype = method),
      linewidth = 0.6
    ) +
    ggplot2::facet_wrap(~ cell, scales = "free_y") +
    ggplot2::scale_x_log10(breaks = n_iter_values) +
    ggplot2::labs(x = "P-IG iterations", y = "Wasserstein-1 distance to grid", colour = NULL, linetype = NULL) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
  ggplot2::ggsave(file.path(fig_dir, "pig_exactness_w1_vs_iterations.pdf"), p_w1, width = 7.5, height = 4)
}

message("Wrote:")
message("  ", raw_file)
message("  ", summary_file)
message("  ", rds_file)
message("Elapsed minutes: ", round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 2))
