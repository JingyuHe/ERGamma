script_path <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) == 0) {
    return(normalizePath("code2026/gamma_shape_inference/make_gamma_shape_manuscript_tables.R"))
  }
  normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)
}

script_dir <- dirname(script_path())
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)
source(file.path(script_dir, "gamma_shape_lib.R"))

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

split_csv_arg <- function(x) {
  x <- trimws(unlist(strsplit(x, ",")))
  x[nzchar(x)]
}

fmt_num <- function(x, digits = 3) {
  ifelse(
    is.na(x),
    "",
    ifelse(abs(x) > 0 & abs(x) < 10^(-digits),
           formatC(x, format = "e", digits = 2),
           formatC(x, format = "f", digits = digits))
  )
}

fmt_int <- function(x) {
  ifelse(is.na(x), "", formatC(round(x), format = "d", big.mark = ","))
}

read_raw_dir <- function(path) {
  raw_file <- file.path(path, "gamma_shape_raw.csv")
  if (!file.exists(raw_file)) {
    stop("Missing raw result file: ", raw_file, call. = FALSE)
  }
  out <- read.csv(raw_file, stringsAsFactors = FALSE)
  out$source_dir <- basename(normalizePath(path, mustWork = TRUE))
  out
}

write_table <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(x, path, row.names = FALSE)
  message("  ", path)
  invisible(x)
}

write_markdown_table <- function(x, path, digits = 3, max_rows = Inf) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (nrow(x) > max_rows) {
    x <- head(x, max_rows)
  }
  y <- x
  for (nm in names(y)) {
    if (is.numeric(y[[nm]])) {
      y[[nm]] <- fmt_num(y[[nm]], digits = digits)
    }
    y[[nm]][is.na(y[[nm]])] <- ""
  }
  lines <- c(
    paste0("| ", paste(names(y), collapse = " | "), " |"),
    paste0("| ", paste(rep("---", ncol(y)), collapse = " | "), " |")
  )
  if (nrow(y) > 0) {
    body <- apply(y, 1, function(row) paste0("| ", paste(row, collapse = " | "), " |"))
    lines <- c(lines, body)
  }
  writeLines(lines, path)
  message("  ", path)
  invisible(path)
}

mc_se <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) <= 1) {
    return(NA_real_)
  }
  stats::sd(x) / sqrt(length(x))
}

summarise_by <- function(d, group_cols) {
  if (nrow(d) == 0) {
    empty <- as.data.frame(setNames(replicate(length(group_cols), logical(0), simplify = FALSE), group_cols))
    empty$reps <- integer(0)
    empty$rows <- integer(0)
    empty$ks_mean <- numeric(0)
    empty$ks_se <- numeric(0)
    empty$w1_mean <- numeric(0)
    empty$w1_se <- numeric(0)
    empty$coverage <- numeric(0)
    empty$mean_error <- numeric(0)
    empty$ess_iter_median <- numeric(0)
    empty$ess_sec_median <- numeric(0)
    empty$runtime_median <- numeric(0)
    empty$runtime_p90 <- numeric(0)
    empty$accept_median <- numeric(0)
    empty$acf1_mean <- numeric(0)
    empty$acf10_mean <- numeric(0)
    return(empty)
  }
  key <- d[group_cols]
  for (nm in names(key)) {
    key[[nm]] <- ifelse(is.na(key[[nm]]), "<NA>", as.character(key[[nm]]))
  }
  groups <- split(d, do.call(paste, c(key, sep = "\r")))
  rows <- lapply(groups, function(g) {
    num <- function(nm) suppressWarnings(as.numeric(g[[nm]]))
    cover_num <- suppressWarnings(as.numeric(g$cover))
    data.frame(
      g[1, group_cols, drop = FALSE],
      reps = length(unique(g$task_id)),
      rows = nrow(g),
      ks_mean = mean(num("ks"), na.rm = TRUE),
      ks_se = mc_se(num("ks")),
      w1_mean = mean(num("w1"), na.rm = TRUE),
      w1_se = mc_se(num("w1")),
      coverage = mean(cover_num, na.rm = TRUE),
      mean_error = mean(num("mean") - num("alpha_true"), na.rm = TRUE),
      ess_iter_median = stats::median(num("ess_per_iter"), na.rm = TRUE),
      ess_sec_median = stats::median(num("ess_per_sec"), na.rm = TRUE),
      runtime_median = stats::median(num("runtime_sec"), na.rm = TRUE),
      runtime_p90 = as.numeric(stats::quantile(num("runtime_sec"), 0.9, na.rm = TRUE, names = FALSE)),
      accept_median = stats::median(num("accept_rate"), na.rm = TRUE),
      acf1_mean = mean(num("acf1"), na.rm = TRUE),
      acf10_mean = mean(num("acf10"), na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

method_rank <- function(x) {
  match(x, c(
    "P-IG Gibbs", "Slice-log-alpha", "Miller-IMH",
    "Miller-Gamma", "Stirling-Gamma", "Numerical grid"
  ))
}

round_report <- function(x, digits = 3) {
  as.numeric(formatC(x, format = "fg", digits = digits))
}

make_pairwise_pig_table <- function(main_summary) {
  exact <- main_summary[
    main_summary$method %in% c("P-IG Gibbs", "Slice-log-alpha", "Miller-IMH"),
    ,
    drop = FALSE
  ]
  if (nrow(exact) == 0 || !("P-IG Gibbs" %in% exact$method)) {
    return(data.frame())
  }
  split_key <- paste(exact$alpha_true, exact$n, sep = "\r")
  by_cell <- split(exact, split_key)
  rows <- lapply(by_cell, function(g) {
    pig <- g[g$method == "P-IG Gibbs", , drop = FALSE]
    comp <- g[g$method != "P-IG Gibbs", , drop = FALSE]
    if (nrow(pig) == 0 || nrow(comp) == 0) {
      return(NULL)
    }
    best_ks <- comp[which.min(comp$ks_mean), , drop = FALSE]
    best_ess <- comp[which.max(comp$ess_sec_median), , drop = FALSE]
    data.frame(
      alpha_true = pig$alpha_true[1],
      n = pig$n[1],
      pig_ks = pig$ks_mean[1],
      best_competitor_ks = best_ks$ks_mean[1],
      best_ks_method = best_ks$method[1],
      pig_ess_sec = pig$ess_sec_median[1],
      best_competitor_ess_sec = best_ess$ess_sec_median[1],
      best_ess_method = best_ess$method[1],
      ess_sec_ratio_pig_to_best = pig$ess_sec_median[1] / best_ess$ess_sec_median[1],
      pig_runtime_median = pig$runtime_median[1],
      best_competitor_runtime_median = best_ess$runtime_median[1],
      stringsAsFactors = FALSE
    )
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) {
    return(data.frame(
      alpha_true = numeric(0),
      n = numeric(0),
      pig_ks = numeric(0),
      best_competitor_ks = numeric(0),
      best_ks_method = character(0),
      pig_ess_sec = numeric(0),
      best_competitor_ess_sec = numeric(0),
      best_ess_method = character(0),
      ess_sec_ratio_pig_to_best = numeric(0),
      pig_runtime_median = numeric(0),
      best_competitor_runtime_median = numeric(0)
    ))
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out[order(out$alpha_true, out$n), , drop = FALSE]
}

make_report <- function(tables, raw, out_file) {
  ok <- raw[raw$status == "ok", , drop = FALSE]
  skipped <- raw[raw$status == "skipped", , drop = FALSE]
  main <- tables$main_accuracy
  main_exact <- main[main$method %in% c("P-IG Gibbs", "Slice-log-alpha", "Miller-IMH"), , drop = FALSE]
  main_approx <- main[main$method %in% c("Miller-Gamma", "Stirling-Gamma"), , drop = FALSE]
  main_pig <- main[main$method == "P-IG Gibbs", , drop = FALSE]
  main_miller_imh <- main[main$method == "Miller-IMH", , drop = FALSE]
  main_stirling <- main[main$method == "Stirling-Gamma", , drop = FALSE]
  main_miller_gamma <- main[main$method == "Miller-Gamma", , drop = FALSE]
  extreme <- tables$extreme_accuracy
  prior_status <- tables$noninteger_prior_status
  trunc <- tables$truncation_pig
  pair <- tables$pig_vs_best_exact_main

  lines <- c(
    "# Gamma Shape Inference Simulation Results",
    "",
    paste0("Generated at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    "",
    "## Scope",
    "",
    paste0("- Raw result rows: ", fmt_int(nrow(raw)), "."),
    paste0("- Successful result rows: ", fmt_int(nrow(ok)), "."),
    paste0("- Skipped result rows: ", fmt_int(nrow(skipped)), "."),
    "- The numerical grid is used as the one-dimensional ground truth, not as a sampling competitor.",
    "- P-IG Gibbs is run only where the posterior delta is integer under the current construction.",
    "",
    "## Main Cross-Regime Simulation",
    ""
  )

  if (nrow(main_pig) > 0) {
    lines <- c(lines, paste0(
      "- P-IG Gibbs mean KS over the 18 main cells is ",
      fmt_num(mean(main_pig$ks_mean, na.rm = TRUE), 4),
      "; median ESS/sec is ",
      fmt_num(stats::median(main_pig$ess_sec_median, na.rm = TRUE), 1),
      "."
    ))
  }
  if (nrow(main_miller_imh) > 0) {
    lines <- c(lines, paste0(
      "- Miller-IMH mean KS over the main cells is ",
      fmt_num(mean(main_miller_imh$ks_mean, na.rm = TRUE), 4),
      "; median ESS/sec is ",
      fmt_num(stats::median(main_miller_imh$ess_sec_median, na.rm = TRUE), 1),
      "."
    ))
  }
  if (nrow(main_miller_gamma) > 0 && nrow(main_stirling) > 0) {
    lines <- c(lines, paste0(
      "- Among deterministic approximations, Miller-Gamma has mean KS ",
      fmt_num(mean(main_miller_gamma$ks_mean, na.rm = TRUE), 4),
      " versus ",
      fmt_num(mean(main_stirling$ks_mean, na.rm = TRUE), 4),
      " for Stirling-Gamma."
    ))
  }
  if (nrow(pair) > 0) {
    lines <- c(lines, paste0(
      "- P-IG Gibbs is the lowest-KS exact sampler in ",
      sum(pair$pig_ks <= pair$best_competitor_ks, na.rm = TRUE),
      " of ",
      nrow(pair),
      " main cells, and its median ESS/sec ratio to the fastest exact competitor is ",
      fmt_num(stats::median(pair$ess_sec_ratio_pig_to_best, na.rm = TRUE), 4),
      "."
    ))
  }

  lines <- c(lines, "", "## Extreme Low-Shape Stress Test", "")
  if (nrow(extreme) > 0) {
    extreme_wide <- do.call(rbind, lapply(split(extreme, extreme$method), function(e) {
      data.frame(
        method = e$method[1],
        ks_mean = mean(e$ks_mean, na.rm = TRUE),
        w1_mean = mean(e$w1_mean, na.rm = TRUE),
        coverage = mean(e$coverage, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }))
    for (m in c("P-IG Gibbs", "Miller-IMH", "Miller-Gamma", "Stirling-Gamma", "Slice-log-alpha")) {
      e <- extreme_wide[extreme_wide$method == m, , drop = FALSE]
      if (nrow(e) > 0) {
        lines <- c(lines, paste0(
          "- ", m, ": mean KS ",
          fmt_num(e$ks_mean[1], 4),
          ", mean W1 ",
          fmt_num(e$w1_mean[1], 4),
          ", coverage ",
          fmt_num(e$coverage[1], 3),
          "."
        ))
      }
    }
  }

  lines <- c(lines, "", "## Prior Robustness", "")
  if (nrow(prior_status) > 0) {
    pig_skip <- prior_status[
      grepl("^P-IG", prior_status$method) & prior_status$status == "skipped",
      ,
      drop = FALSE
    ]
    if (nrow(pig_skip) > 0) {
      lines <- c(lines, paste0(
        "- P-IG Gibbs skipped ",
        sum(pig_skip$rows),
        " prior rows because the current integer-delta augmentation does not cover noninteger posterior delta."
      ))
    }
  }

  lines <- c(lines, "", "## Truncation Sensitivity", "")
  if (nrow(trunc) > 0) {
    n200 <- trunc[trunc$pig_N == 200, , drop = FALSE]
    n1000 <- trunc[trunc$pig_N == 1000, , drop = FALSE]
    if (nrow(n200) > 0 && nrow(n1000) > 0) {
      lines <- c(lines, paste0(
        "- With N=200, mean truncation-suite KS is ",
        fmt_num(mean(n200$ks_mean, na.rm = TRUE), 4),
        "; with N=1000 it is ",
        fmt_num(mean(n1000$ks_mean, na.rm = TRUE), 4),
        ". Median runtime increases from ",
        fmt_num(stats::median(n200$runtime_median, na.rm = TRUE), 2),
        "s to ",
        fmt_num(stats::median(n1000$runtime_median, na.rm = TRUE), 2),
        "s."
      ))
    }
  }

  lines <- c(lines, "", "## Manuscript-Ready Interpretation", "")
  lines <- c(lines,
    "The simulation now separates two questions that were conflated in the earlier draft: numerical correctness of the target distribution and computational efficiency of competing samplers. The grid calculation supplies a stable one-dimensional reference for each generated data set, so KS and Wasserstein errors are interpretable as errors relative to the actual posterior rather than discrepancies among Monte Carlo methods.",
    "",
    "The exact P-IG Gibbs sampler is robust in the sense that it targets the grid posterior across small, moderate, and large shapes when the integer posterior-delta condition is satisfied. However, the cost grows directly with posterior delta because each iteration samples a sum of delta P-IG variables. The comparison with slice sampling and Miller-based independence MH should therefore be reported in ESS/sec, not only ESS/iteration.",
    "",
    "Miller-Gamma is the strongest deterministic approximation in these experiments and Miller-IMH is the most important exact competitor because it uses that approximation only as a proposal while preserving the target. Stirling-Gamma is useful as a historical baseline but should not be presented as competitive in the low-shape regimes if its KS/Wasserstein errors are visibly larger.",
    "",
    "The prior robustness experiment exposes the current limitation of the augmentation: noninteger posterior delta is not covered by the present Gibbs construction. The paper should either state the integer-delta condition explicitly or add a separate noninteger construction before claiming full robustness to arbitrary informative priors."
  )

  writeLines(lines, out_file)
  message("  ", out_file)
  invisible(lines)
}

main <- function() {
  default_bench <- file.path(repo_root, "code2026", "gamma_shape_inference", "results", "all_full_benchmarks_no_pig")
  default_pig <- file.path(repo_root, "code2026", "gamma_shape_inference", "results", "all_full_pig_only")
  result_dirs <- split_csv_arg(arg_value("--result-dirs", paste(default_bench, default_pig, sep = ",")))
  out_dir <- arg_value(
    "--out-dir",
    file.path(repo_root, "code2026", "gamma_shape_inference", "results", "manuscript_tables")
  )
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  raw <- do.call(rbind, lapply(result_dirs, read_raw_dir))
  raw$row_id <- seq_len(nrow(raw))
  raw_ok <- raw[raw$status == "ok" & raw$method != "Numerical grid", , drop = FALSE]
  combined_summary <- summarise_gamma_shape_results(raw)

  main_accuracy <- summarise_by(
    raw_ok[raw_ok$suite == "main", , drop = FALSE],
    c("alpha_true", "n", "method")
  )
  if (nrow(main_accuracy) > 0) {
    main_accuracy <- main_accuracy[order(main_accuracy$alpha_true, main_accuracy$n, method_rank(main_accuracy$method)), ]
  }

  main_efficiency <- main_accuracy[
    main_accuracy$method %in% c("P-IG Gibbs", "Slice-log-alpha", "Miller-IMH"),
    c("alpha_true", "n", "method", "reps", "ess_iter_median", "ess_sec_median",
      "runtime_median", "runtime_p90", "accept_median", "acf1_mean", "acf10_mean"),
    drop = FALSE
  ]

  damsleth_accuracy <- summarise_by(
    raw_ok[raw_ok$suite == "damsleth", , drop = FALSE],
    c("n", "method")
  )
  if (nrow(damsleth_accuracy) > 0) {
    damsleth_accuracy <- damsleth_accuracy[order(damsleth_accuracy$n, method_rank(damsleth_accuracy$method)), ]
  }

  extreme_accuracy <- summarise_by(
    raw_ok[raw_ok$suite == "extreme", , drop = FALSE],
    c("alpha_true", "n", "method")
  )
  if (nrow(extreme_accuracy) > 0) {
    extreme_accuracy <- extreme_accuracy[order(extreme_accuracy$alpha_true, method_rank(extreme_accuracy$method)), ]
  }

  prior_robustness <- summarise_by(
    raw_ok[raw_ok$suite == "prior", , drop = FALSE],
    c("prior_delta", "prior_center", "method")
  )
  if (nrow(prior_robustness) > 0) {
    prior_robustness <- prior_robustness[
      order(prior_robustness$prior_delta, prior_robustness$prior_center, method_rank(prior_robustness$method)),
      ,
      drop = FALSE
    ]
  }

  truncation_pig <- summarise_by(
    raw_ok[raw_ok$suite == "truncation" & grepl("^P-IG Gibbs", raw_ok$method), , drop = FALSE],
    c("alpha_true", "n", "pig_N", "method")
  )
  if (nrow(truncation_pig) > 0) {
    truncation_pig <- truncation_pig[order(truncation_pig$alpha_true, truncation_pig$n, truncation_pig$pig_N), ]
  }

  method_overall <- summarise_by(
    raw_ok[raw_ok$suite %in% c("main", "extreme"), , drop = FALSE],
    c("suite", "method")
  )
  if (nrow(method_overall) > 0) {
    method_overall <- method_overall[order(method_overall$suite, method_rank(method_overall$method)), ]
  }

  noninteger_prior_status <- summarise_by(
    raw[raw$suite == "prior", , drop = FALSE],
    c("prior_delta", "prior_center", "method", "status")
  )
  if (nrow(noninteger_prior_status) > 0) {
    noninteger_prior_status <- noninteger_prior_status[
      order(noninteger_prior_status$prior_delta, noninteger_prior_status$prior_center,
            method_rank(noninteger_prior_status$method), noninteger_prior_status$status),
      ,
      drop = FALSE
    ]
  }

  pig_vs_best_exact_main <- make_pairwise_pig_table(main_accuracy)

  tables <- list(
    damsleth_accuracy = damsleth_accuracy,
    main_accuracy = main_accuracy,
    main_efficiency = main_efficiency,
    pig_vs_best_exact_main = pig_vs_best_exact_main,
    extreme_accuracy = extreme_accuracy,
    prior_robustness = prior_robustness,
    truncation_pig = truncation_pig,
    method_overall = method_overall,
    noninteger_prior_status = noninteger_prior_status
  )

  message("Writing manuscript tables:")
  for (nm in names(tables)) {
    write_table(tables[[nm]], file.path(out_dir, paste0(nm, ".csv")))
  }
  write_markdown_table(
    main_accuracy,
    file.path(out_dir, "main_accuracy.md"),
    digits = 4
  )
  write_markdown_table(
    main_efficiency,
    file.path(out_dir, "main_efficiency.md"),
    digits = 3
  )
  write_markdown_table(
    pig_vs_best_exact_main,
    file.path(out_dir, "pig_vs_best_exact_main.md"),
    digits = 4
  )
  write_markdown_table(
    extreme_accuracy,
    file.path(out_dir, "extreme_accuracy.md"),
    digits = 4
  )
  make_report(tables, raw, file.path(out_dir, "simulation_results_summary.md"))

  saveRDS(
    list(raw = raw, tables = tables, result_dirs = result_dirs),
    file.path(out_dir, "gamma_shape_manuscript_tables.rds")
  )
  message("  ", file.path(out_dir, "gamma_shape_manuscript_tables.rds"))
  saveRDS(
    list(raw = raw, summary = combined_summary, result_dirs = result_dirs),
    file.path(out_dir, "gamma_shape_combined_results.rds")
  )
  message("  ", file.path(out_dir, "gamma_shape_combined_results.rds"))
}

main()
