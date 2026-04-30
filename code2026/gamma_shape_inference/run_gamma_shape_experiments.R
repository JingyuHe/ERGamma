script_path <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) == 0) {
    return(normalizePath("code2026/gamma_shape_inference/run_gamma_shape_experiments.R"))
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
  if (tolower(as.character(x)) %in% c("inf", "infinity")) {
    return(Inf)
  }
  as.numeric(x)
}

suite <- arg_value("--suite", "smoke")
profile <- arg_value("--profile", if (suite == "smoke") "smoke" else "smoke")
reps_arg <- arg_value("--reps", NULL)
reps <- if (is.null(reps_arg)) NULL else as.integer(reps_arg)

out_dir <- arg_value(
  "--out-dir",
  file.path(repo_root, "code2026", "gamma_shape_inference", "results", paste0(suite, "_", profile))
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
tasks <- make_tasks(suite = suite, profile = profile, reps = reps)

message("Gamma shape experiment")
message("  repo: ", repo_root)
message("  suite/profile: ", suite, "/", profile)
message("  tasks: ", nrow(tasks))
message("  output: ", out_dir)
message("  P-IG: ", settings$run_pig, ", pig_N=", settings$pig_N,
        ", pig_max_delta=", settings$pig_max_delta)
message("  MCMC iter/burnin/thin: ", settings$n_iter, "/", settings$burnin, "/", settings$thin)

run_indices <- seq_len(nrow(tasks))
run_one <- function(i) {
  run_gamma_shape_task(tasks[i, , drop = FALSE], settings)
}

if (cores > 1 && .Platform$OS.type != "windows") {
  results <- parallel::mclapply(run_indices, run_one, mc.cores = cores, mc.set.seed = TRUE)
} else {
  results <- lapply(run_indices, run_one)
}

raw <- do.call(rbind, results)
rownames(raw) <- NULL
summary <- summarise_gamma_shape_results(raw)

raw_file <- file.path(out_dir, "gamma_shape_raw.csv")
summary_file <- file.path(out_dir, "gamma_shape_summary.csv")
rds_file <- file.path(out_dir, "gamma_shape_results.rds")

write.csv(raw, raw_file, row.names = FALSE)
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
