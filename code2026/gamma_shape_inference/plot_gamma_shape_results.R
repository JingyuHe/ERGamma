script_path <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) == 0) {
    return(normalizePath("code2026/gamma_shape_inference/plot_gamma_shape_results.R"))
  }
  normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)
}

script_dir <- dirname(script_path())
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)
source(file.path(script_dir, "gamma_shape_lib.R"))
gamma_shape_init_libpaths(repo_root)
require_pkg("ggplot2")

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

result_rds <- arg_value(
  "--result-rds",
  file.path(repo_root, "code2026", "gamma_shape_inference", "results", "smoke_smoke", "gamma_shape_results.rds")
)
res <- readRDS(result_rds)
raw <- res$raw
summary <- res$summary

fig_dir <- arg_value("--fig-dir", file.path(dirname(result_rds), "figures"))
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

theme_gamma_shape <- function() {
  ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold"),
      strip.background = ggplot2::element_rect(fill = "grey95", colour = "grey70")
    )
}

plot_accuracy_heatmap <- function(summary, fig_dir) {
  dat <- summary[summary$suite == "main" & summary$method != "Numerical grid", , drop = FALSE]
  if (!nrow(dat) || !any(dat$method == "P-IG Gibbs")) {
    return(NULL)
  }
  pig <- dat[dat$method == "P-IG Gibbs", c("alpha_true", "n", "ks_mean")]
  names(pig)[3] <- "pig_ks"
  comp <- dat[dat$method %in% c("Slice-log-alpha", "Miller-IMH", "Miller-Gamma", "Stirling-Gamma"),
              c("alpha_true", "n", "ks_mean"), drop = FALSE]
  if (!nrow(comp)) return(NULL)
  best <- aggregate(ks_mean ~ alpha_true + n, comp, min, na.rm = TRUE)
  names(best)[3] <- "best_comp_ks"
  rel <- merge(pig, best, by = c("alpha_true", "n"), all = TRUE)
  rel$ratio <- rel$best_comp_ks / rel$pig_ks
  rel$label <- ifelse(is.finite(rel$ratio), sprintf("%.2f", rel$ratio), "")

  p <- ggplot2::ggplot(rel, ggplot2::aes(factor(n), factor(alpha_true), fill = ratio)) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::geom_text(ggplot2::aes(label = label), size = 3) +
    ggplot2::scale_fill_gradient2(
      low = "#b2182b", mid = "white", high = "#2166ac",
      midpoint = 1, na.value = "grey90",
      name = "Best competitor KS / P-IG KS"
    ) +
    ggplot2::labs(
      x = "n",
      y = expression(alpha[true]),
      title = "P-IG Accuracy Relative to Best Competitor"
    ) +
    theme_gamma_shape()
  file <- file.path(fig_dir, "gamma_shape_main_ks_heatmap.pdf")
  ggplot2::ggsave(file, p, width = 6.5, height = 4.5)
  file
}

plot_ess_heatmap <- function(summary, fig_dir) {
  dat <- summary[summary$suite == "main" &
                   summary$method %in% c("P-IG Gibbs", "Slice-log-alpha", "Miller-IMH"), , drop = FALSE]
  if (!nrow(dat) || !all(c("P-IG Gibbs", "Slice-log-alpha") %in% dat$method)) {
    return(NULL)
  }
  pig <- dat[dat$method == "P-IG Gibbs", c("alpha_true", "n", "ess_sec_median")]
  names(pig)[3] <- "pig_ess_sec"
  sl <- dat[dat$method == "Slice-log-alpha", c("alpha_true", "n", "ess_sec_median")]
  names(sl)[3] <- "slice_ess_sec"
  rel <- merge(pig, sl, by = c("alpha_true", "n"), all = TRUE)
  rel$ratio <- rel$pig_ess_sec / rel$slice_ess_sec
  rel$label <- ifelse(is.finite(rel$ratio), sprintf("%.2f", rel$ratio), "")

  p <- ggplot2::ggplot(rel, ggplot2::aes(factor(n), factor(alpha_true), fill = ratio)) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::geom_text(ggplot2::aes(label = label), size = 3) +
    ggplot2::scale_fill_gradient2(
      low = "#b2182b", mid = "white", high = "#2166ac",
      midpoint = 1, na.value = "grey90",
      name = "P-IG ESS/sec / slice ESS/sec"
    ) +
    ggplot2::labs(
      x = "n",
      y = expression(alpha[true]),
      title = "P-IG Efficiency Relative to Slice Sampling"
    ) +
    theme_gamma_shape()
  file <- file.path(fig_dir, "gamma_shape_main_ess_sec_heatmap.pdf")
  ggplot2::ggsave(file, p, width = 6.5, height = 4.5)
  file
}

plot_coverage <- function(summary, fig_dir) {
  dat <- summary[summary$suite == "main" &
                   summary$method %in% c("Numerical grid", "P-IG Gibbs", "Slice-log-alpha",
                                         "Miller-IMH", "Miller-Gamma", "Stirling-Gamma"), , drop = FALSE]
  if (!nrow(dat)) return(NULL)
  dat$cell <- paste0("a=", dat$alpha_true, ", n=", dat$n)
  p <- ggplot2::ggplot(dat, ggplot2::aes(cell, coverage, colour = method, group = method)) +
    ggplot2::geom_hline(yintercept = 0.95, linetype = 2, colour = "grey45") +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(x = "simulation cell", y = "empirical 95% coverage", title = "Coverage by Method") +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  file <- file.path(fig_dir, "gamma_shape_main_coverage.pdf")
  ggplot2::ggsave(file, p, width = 9, height = 4.5)
  file
}

plot_extreme_density <- function(raw, fig_dir) {
  extreme_reps <- raw$rep[raw$suite == "extreme"]
  if (!length(extreme_reps)) return(NULL)
  dat <- raw[raw$suite == "extreme" & raw$rep == min(extreme_reps), , drop = FALSE]
  if (!nrow(dat)) return(NULL)
  first_rows <- dat[!duplicated(dat$alpha_true), , drop = FALSE]
  density_rows <- list()
  for (i in seq_len(nrow(first_rows))) {
    row <- first_rows[i, ]
    truth <- xi2_grid_truth(row$delta_post, row$mu_ratio, grid_size = 4001)
    alpha <- truth$alpha
    keep <- alpha <= stats::quantile(alpha, 0.995)
    alpha <- alpha[keep]
    density_rows[[length(density_rows) + 1]] <- data.frame(
      alpha_true = row$alpha_true,
      alpha = alpha,
      density = truth$density[keep],
      method = "Numerical grid"
    )
    mp <- miller_gamma_params(row$delta_post, row$mu_ratio, init = truth$mode)
    if (isTRUE(mp$ok)) {
      density_rows[[length(density_rows) + 1]] <- data.frame(
        alpha_true = row$alpha_true,
        alpha = alpha,
        density = stats::dgamma(alpha, shape = mp$shape, rate = mp$rate),
        method = "Miller-Gamma"
      )
    }
    sp <- stirling_gamma_params(row$delta_post, row$mu_ratio)
    if (isTRUE(sp$ok)) {
      density_rows[[length(density_rows) + 1]] <- data.frame(
        alpha_true = row$alpha_true,
        alpha = alpha,
        density = stats::dgamma(alpha, shape = sp$shape, rate = sp$rate),
        method = "Stirling-Gamma"
      )
    }
  }
  dd <- do.call(rbind, density_rows)
  p <- ggplot2::ggplot(dd, ggplot2::aes(alpha, density, colour = method, linetype = method)) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::facet_wrap(~ alpha_true, scales = "free", labeller = ggplot2::label_bquote(alpha[true] == .(alpha_true))) +
    ggplot2::labs(x = expression(alpha), y = "density", title = "Extreme Low-Shape Approximation Stress Test") +
    theme_gamma_shape()
  file <- file.path(fig_dir, "gamma_shape_extreme_density.pdf")
  ggplot2::ggsave(file, p, width = 7, height = 3.8)
  file
}

files <- c(
  plot_accuracy_heatmap(summary, fig_dir),
  plot_ess_heatmap(summary, fig_dir),
  plot_coverage(summary, fig_dir),
  plot_extreme_density(raw, fig_dir)
)
files <- files[!vapply(files, is.null, logical(1))]

message("Wrote figures:")
for (f in files) message("  ", f)
