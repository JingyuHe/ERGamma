script_path <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) == 0) {
    return(normalizePath("code2026/gamma_shape_inference/run_gamma_shape_experiments_checkpoint.R"))
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

as_bool <- function(x) {
  if (is.logical(x)) return(x)
  tolower(as.character(x)) %in% c("1", "true", "t", "yes", "y")
}

as_num_or_inf <- function(x) {
  if (tolower(as.character(x)) %in% c("inf", "infinity")) return(Inf)
  as.numeric(x)
}

suite <- arg_value("--suite", "smoke")
profile <- arg_value("--profile", if (suite == "smoke") "smoke" else "smoke")
reps_arg <- arg_value("--reps", NULL)
reps <- if (is.null(reps_arg)) NULL else as.integer(reps_arg)

out_dir <- arg_value(
  "--out-dir",
  file.path(repo_root, "code2026", "gamma_shape_inference", "results", paste0(suite, "_", profile, "_checkpoint"))
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

settings <- list(
  seed = as.integer(arg_value("--seed", "20260430")),
  grid_size = as.integer(arg_value("--grid-size", if (profile == "full") "20001" else "4001")),
  tail_nats = as.numeric(arg_value("--tail-nats", "50")),
  n_iter = as.integer(arg_value("--iter", if (profile == "full") "5000" else "250")),
  burnin = as.integer(arg_value("--burnin", if (profile == "full") "1000" else "50")),
  thin = as.integer(arg_value("--thin", "1")),
  pig_N = as.integer(arg_value("--pig-N", if (profile == "full") "200" else "20")),
  pig_max_delta = as_num_or_inf(arg_value("--pig-max-delta", if (profile == "full") "50" else "30")),
  truncation_N = as.integer(strsplit(arg_value("--truncation-N", if (profile == "full") "50,100,200,500,1000" else "10,20"), ",")[[1]]),
  run_pig = as_bool(arg_value("--run-pig", "true")),
  run_slice = as_bool(arg_value("--run-slice", "true")),
  run_approximations = as_bool(arg_value("--run-approximations", "true")),
  run_miller_imh = as_bool(arg_value("--run-miller-imh", "true")),
  slice_w = as.numeric(arg_value("--slice-w", "1")),
  slice_m = as.integer(arg_value("--slice-m", "100"))
)

if (settings$n_iter <= settings$burnin) {
  stop("--iter must exceed --burnin", call. = FALSE)
}

cores <- as.integer(arg_value("--cores", "1"))
chunk_size <- as.integer(arg_value("--chunk-size", as.character(max(cores, 1) * 4)))
resume <- as_bool(arg_value("--resume", "true"))

tasks <- make_tasks(suite = suite, profile = profile, reps = reps)
raw_file <- file.path(out_dir, "gamma_shape_raw.csv")
summary_file <- file.path(out_dir, "gamma_shape_summary.csv")
rds_file <- file.path(out_dir, "gamma_shape_results.rds")
progress_file <- file.path(out_dir, "progress.csv")

done_ids <- integer()
raw_existing <- NULL
if (resume && file.exists(raw_file)) {
  raw_existing <- read.csv(raw_file, stringsAsFactors = FALSE)
  done_ids <- sort(unique(raw_existing$task_id))
}

todo <- tasks[!(tasks$task_id %in% done_ids), , drop = FALSE]
start_time <- Sys.time()
initial_done <- length(done_ids)
total <- nrow(tasks)

message("Gamma shape checkpoint experiment")
message("  repo: ", repo_root)
message("  suite/profile: ", suite, "/", profile)
message("  tasks: ", total, " total, ", initial_done, " already done, ", nrow(todo), " remaining")
message("  output: ", out_dir)
message("  cores/chunk_size: ", cores, "/", chunk_size)
message("  P-IG: ", settings$run_pig, ", pig_N=", settings$pig_N,
        ", pig_max_delta=", settings$pig_max_delta)
message("  MCMC iter/burnin/thin: ", settings$n_iter, "/", settings$burnin, "/", settings$thin)

append_raw <- function(chunk_raw) {
  write.table(
    chunk_raw,
    file = raw_file,
    sep = ",",
    row.names = FALSE,
    col.names = !file.exists(raw_file),
    append = file.exists(raw_file)
  )
}

write_progress <- function(completed, chunk_elapsed, last_task_id) {
  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  newly_done <- completed - initial_done
  rate <- if (newly_done > 0) newly_done / elapsed else NA_real_
  remaining <- total - completed
  eta_sec <- if (is.finite(rate) && rate > 0) remaining / rate else NA_real_
  row <- data.frame(
    time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    completed = completed,
    total = total,
    remaining = remaining,
    percent = 100 * completed / total,
    elapsed_sec = elapsed,
    eta_sec = eta_sec,
    chunk_elapsed_sec = chunk_elapsed,
    last_task_id = last_task_id
  )
  write.table(
    row,
    file = progress_file,
    sep = ",",
    row.names = FALSE,
    col.names = !file.exists(progress_file),
    append = file.exists(progress_file)
  )
  message(sprintf(
    "[%s] completed %d/%d (%.1f%%); chunk %.1fs; elapsed %.1f min; ETA %.1f min",
    row$time, completed, total, row$percent, chunk_elapsed, elapsed / 60, eta_sec / 60
  ))
}

run_one <- function(task_row) {
  tryCatch(
    run_gamma_shape_task(task_row, settings),
    error = function(e) {
      cbind(
        task_row[rep(1, 1), c(
          "task_id", "suite", "scenario", "rep", "alpha_true", "n",
          "prior_delta", "prior_center", "fixed_stats"
        )],
        x_a = NA_real_, x_g = NA_real_,
        delta_post = NA_real_, mu_ratio = NA_real_, truth_mode = NA_real_,
        skip_row("TASK-FAILED", e$message),
        stringsAsFactors = FALSE
      )
    }
  )
}

if (nrow(todo) > 0) {
  chunks <- split(seq_len(nrow(todo)), ceiling(seq_len(nrow(todo)) / chunk_size))
  for (idx in chunks) {
    chunk_tasks <- todo[idx, , drop = FALSE]
    t0 <- proc.time()[["elapsed"]]
    if (cores > 1 && .Platform$OS.type != "windows") {
      chunk_results <- parallel::mclapply(
        split(chunk_tasks, seq_len(nrow(chunk_tasks))),
        run_one,
        mc.cores = cores,
        mc.set.seed = TRUE
      )
    } else {
      chunk_results <- lapply(split(chunk_tasks, seq_len(nrow(chunk_tasks))), run_one)
    }
    chunk_raw <- do.call(rbind, chunk_results)
    rownames(chunk_raw) <- NULL
    append_raw(chunk_raw)
    done_ids <- c(done_ids, chunk_tasks$task_id)
    write_progress(
      completed = length(unique(done_ids)),
      chunk_elapsed = proc.time()[["elapsed"]] - t0,
      last_task_id = max(chunk_tasks$task_id)
    )
  }
}

raw <- read.csv(raw_file, stringsAsFactors = FALSE)
summary <- summarise_gamma_shape_results(raw)
write.csv(summary, summary_file, row.names = FALSE)
saveRDS(
  list(
    raw = raw,
    summary = summary,
    tasks = tasks,
    settings = settings,
    suite = suite,
    profile = profile,
    repo_root = repo_root
  ),
  rds_file
)

message("Wrote:")
message("  ", raw_file)
message("  ", summary_file)
message("  ", rds_file)
