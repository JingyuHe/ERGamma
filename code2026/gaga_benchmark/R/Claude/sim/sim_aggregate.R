#!/usr/bin/env Rscript
# Aggregate sim results into the paper-ready summary tables and figures.

src_dir <- normalizePath(dirname(sub("^--file=", "",
  commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][1])))
results_dir <- file.path(src_dir, "results")
agg_dir <- file.path(src_dir, "summary")
dir.create(agg_dir, recursive = TRUE, showWarnings = FALSE)
af <- function(f) file.path(agg_dir, f)

scenario_files <- list.files(results_dir, pattern = "^scenario_[A-D]_rep_",
                             full.names = TRUE)
sweep_files    <- list.files(results_dir, pattern = "^sweep_nua_",
                             full.names = TRUE)

# ---------------- main scenario aggregation -----------------------------
rows <- list(); cross <- list(); reps <- list()
for (f in scenario_files) {
  obj <- readRDS(f)
  pm <- obj$metrics$per_method
  pm$scenario <- obj$scenario; pm$rep <- obj$rep
  pm$elapsed_sec <- obj$elapsed_sec
  rows[[length(rows) + 1]] <- pm
  cr <- obj$metrics$cross
  cr$scenario <- obj$scenario; cr$rep <- obj$rep
  cross[[length(cross) + 1]] <- cr
  reps[[length(reps) + 1]] <- data.frame(
    scenario = obj$scenario, rep = obj$rep,
    spec_n = obj$spec$n, spec_G = obj$spec$G,
    spec_b_alpha = obj$spec$b_alpha, spec_nu_alpha = obj$spec$nu_alpha,
    spec_a0 = obj$spec$a0, spec_nu = obj$spec$nu,
    spec_p_de = obj$spec$p_de,
    est_a0 = obj$hyper_estimated$alpha0,
    est_nu = obj$hyper_estimated$nu,
    est_b_alpha = obj$hyper_estimated$b_alpha,
    est_nu_alpha = obj$hyper_estimated$nu_alpha,
    n_de_truth = sum(obj$truth),
    elapsed_sec = obj$elapsed_sec, stringsAsFactors = FALSE)
}
per_method <- do.call(rbind, rows)
cross_all <- do.call(rbind, cross)
spec_all <- do.call(rbind, reps)
write.csv(per_method, af("scenarios_per_method_raw.csv"), row.names = FALSE)
write.csv(cross_all,  af("scenarios_cross_raw.csv"),     row.names = FALSE)
write.csv(spec_all,   af("scenarios_specs_raw.csv"),    row.names = FALSE)

mean_se <- function(x) c(mean = mean(x, na.rm = TRUE),
                         se = sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))

# Per-method agg: power, FDR, AUC by scenario × method
agg_per_method <- aggregate(
  cbind(n_de, tp, fp, empirical_fdr, power, auc) ~ scenario + method,
  data = per_method, FUN = mean)
write.csv(agg_per_method, af("scenarios_per_method_mean.csv"), row.names = FALSE)

# Per-method standard errors
agg_per_method_se <- aggregate(
  cbind(empirical_fdr, power, auc) ~ scenario + method,
  data = per_method,
  FUN = function(v) sd(v) / sqrt(length(v)))
write.csv(agg_per_method_se, af("scenarios_per_method_se.csv"), row.names = FALSE)

# Cross-method agg
agg_cross <- aggregate(
  cbind(mean_abs_pp_diff, max_abs_pp_diff, spearman_score, jaccard_calls) ~
    scenario + pair, data = cross_all, FUN = mean)
write.csv(agg_cross, af("scenarios_cross_mean.csv"), row.names = FALSE)
agg_cross_se <- aggregate(
  cbind(mean_abs_pp_diff, spearman_score) ~ scenario + pair,
  data = cross_all, FUN = function(v) sd(v) / sqrt(length(v)))
write.csv(agg_cross_se, af("scenarios_cross_se.csv"), row.names = FALSE)

# Hyperparameter recovery: oracle vs estimated
hyper_recovery <- aggregate(
  cbind(est_a0, est_nu, est_b_alpha, est_nu_alpha) ~ scenario,
  data = spec_all, FUN = function(v) round(mean(v), 3))
hyper_recovery <- merge(hyper_recovery,
  unique(spec_all[, c("scenario", "spec_a0", "spec_nu",
                      "spec_b_alpha", "spec_nu_alpha")]),
  by = "scenario")
write.csv(hyper_recovery, af("hyperparam_recovery.csv"), row.names = FALSE)

# Headline summary table
headline <- data.frame()
for (s in c("A", "B", "C", "D")) {
  pm <- per_method[per_method$scenario == s, ]
  cr <- cross_all[cross_all$scenario == s, ]
  spec_s <- unique(spec_all[spec_all$scenario == s,
                            c("spec_n", "spec_nu_alpha", "spec_a0")])[1, ]
  rb <- function(method) {
    sub <- pm[pm$method == method, ]
    paste0(sprintf("%.3f", mean(sub$power)), " (",
           sprintf("%.3f", sd(sub$power) / sqrt(nrow(sub))), ")")
  }
  headline <- rbind(headline, data.frame(
    scenario = s,
    n_per_grp = spec_s$spec_n,
    nu_alpha = spec_s$spec_nu_alpha,
    a0 = spec_s$spec_a0,
    gold_pwr = rb("gold"),
    gaga_pwr = rb("gaga"),
    pig_pwr = rb("pig"),
    gaga_vs_gold_meanAbsDpp = sprintf("%.4g",
      mean(cr$mean_abs_pp_diff[cr$pair == "gaga_vs_gold"])),
    pig_vs_gold_meanAbsDpp = sprintf("%.4g",
      mean(cr$mean_abs_pp_diff[cr$pair == "pig_vs_gold"])),
    gaga_spearman_gold = sprintf("%.4f",
      mean(cr$spearman_score[cr$pair == "gaga_vs_gold"])),
    pig_spearman_gold = sprintf("%.4f",
      mean(cr$spearman_score[cr$pair == "pig_vs_gold"])),
    stringsAsFactors = FALSE
  ))
}
write.csv(headline, af("headline_table.csv"), row.names = FALSE)

# Format an FDR/power/auc table (mean ± se) per scenario × method
fmt_msd <- function(m, se) sprintf("%.3f ± %.3f", m, se)
fdr_pwr_auc <- data.frame()
for (s in c("A", "B", "C", "D")) {
  for (m in c("gold", "gaga", "pig")) {
    sub <- per_method[per_method$scenario == s & per_method$method == m, ]
    fdr_pwr_auc <- rbind(fdr_pwr_auc, data.frame(
      scenario = s, method = m,
      empirical_fdr = fmt_msd(mean(sub$empirical_fdr),
                              sd(sub$empirical_fdr) / sqrt(nrow(sub))),
      power = fmt_msd(mean(sub$power), sd(sub$power) / sqrt(nrow(sub))),
      auc = fmt_msd(mean(sub$auc, na.rm = TRUE),
                    sd(sub$auc, na.rm = TRUE) /
                      sqrt(sum(!is.na(sub$auc)))),
      stringsAsFactors = FALSE
    ))
  }
}
write.csv(fdr_pwr_auc, af("fdr_power_auc_table.csv"), row.names = FALSE)

# ---------------- alpha sweep aggregation -------------------------------
sweep_rows <- list(); sweep_cross <- list()
for (f in sweep_files) {
  obj <- readRDS(f)
  cr <- obj$metrics$cross
  cr$nu_alpha <- obj$nu_alpha; cr$rep <- obj$rep
  sweep_cross[[length(sweep_cross) + 1]] <- cr
  pm <- obj$metrics$per_method
  pm$nu_alpha <- obj$nu_alpha; pm$rep <- obj$rep
  sweep_rows[[length(sweep_rows) + 1]] <- pm
}
sweep_pm <- do.call(rbind, sweep_rows)
sweep_cr <- do.call(rbind, sweep_cross)
write.csv(sweep_pm, af("sweep_per_method_raw.csv"), row.names = FALSE)
write.csv(sweep_cr, af("sweep_cross_raw.csv"),     row.names = FALSE)

agg_sweep_cross <- aggregate(
  cbind(mean_abs_pp_diff, spearman_score, jaccard_calls) ~ nu_alpha + pair,
  data = sweep_cr, FUN = mean)
write.csv(agg_sweep_cross, af("sweep_cross_mean.csv"), row.names = FALSE)

agg_sweep_pm <- aggregate(
  cbind(empirical_fdr, power, auc) ~ nu_alpha + method,
  data = sweep_pm, FUN = mean)
write.csv(agg_sweep_pm, af("sweep_per_method_mean.csv"), row.names = FALSE)

# ---------------- plots --------------------------------------------------
plot_scatter_pp <- function() {
  pdf(af("scatter_pp_vs_gold.pdf"), width = 9, height = 7)
  par(mfrow = c(2, 4), mar = c(4.2, 4.2, 2.5, 1))
  for (s in c("A", "B", "C", "D")) {
    files <- list.files(results_dir, pattern = paste0("^scenario_", s, "_rep_"),
                        full.names = TRUE)
    pp_g <- pp_gg <- pp_p <- numeric(0)
    for (f in files) {
      obj <- readRDS(f)
      pp_g <- c(pp_g, obj$pp_gold[, 1])
      pp_gg <- c(pp_gg, obj$pp_gaga[, 1])
      pp_p <- c(pp_p, obj$pp_pig[, 1])
    }
    plot(pp_g, pp_gg, pch = 16, cex = 0.4,
         col = adjustcolor("#D95F02", 0.4),
         main = paste0("Scenario ", s, ": gaga vs gold"),
         xlab = "pp_gold(EE)", ylab = "pp_gaga(EE)",
         xlim = c(0, 1), ylim = c(0, 1))
    abline(0, 1, lwd = 1.5, lty = 2, col = "grey40")
    plot(pp_g, pp_p, pch = 16, cex = 0.4,
         col = adjustcolor("#1B9E77", 0.4),
         main = paste0("Scenario ", s, ": PIG-MCMC vs gold"),
         xlab = "pp_gold(EE)", ylab = "pp_pig(EE)",
         xlim = c(0, 1), ylim = c(0, 1))
    abline(0, 1, lwd = 1.5, lty = 2, col = "grey40")
  }
  dev.off()
}
plot_scatter_pp()

plot_alpha_sweep <- function() {
  ag <- agg_sweep_cross
  pdf(af("alpha_sweep_pp_diff.pdf"), width = 7, height = 5)
  par(mar = c(4.5, 4.5, 2, 1))
  ag_g <- ag[ag$pair == "gaga_vs_gold", ]
  ag_p <- ag[ag$pair == "pig_vs_gold", ]
  ymax <- max(c(ag_g$mean_abs_pp_diff, ag_p$mean_abs_pp_diff))
  plot(ag_g$nu_alpha, ag_g$mean_abs_pp_diff, log = "x", type = "b",
       col = "#D95F02", lwd = 2, pch = 16, ylim = c(0, ymax * 1.1),
       xlab = expression(nu[alpha]),
       ylab = "mean | pp(EE) - gold |",
       main = "Stirling vs PIG-MCMC bias against gold standard")
  lines(ag_p$nu_alpha, ag_p$mean_abs_pp_diff, type = "b", lwd = 2,
        col = "#1B9E77", pch = 17)
  legend("topright",
         legend = c("gaga (Stirling)", "PIG-MCMC (slice)"),
         col = c("#D95F02", "#1B9E77"), pch = c(16, 17), lwd = 2, bty = "n")
  dev.off()
  png(af("alpha_sweep_pp_diff.png"), width = 1400, height = 1000, res = 200)
  par(mar = c(4.5, 4.5, 2, 1))
  plot(ag_g$nu_alpha, ag_g$mean_abs_pp_diff, log = "x", type = "b",
       col = "#D95F02", lwd = 2, pch = 16, ylim = c(0, ymax * 1.1),
       xlab = expression(nu[alpha]),
       ylab = "mean | pp(EE) - gold |",
       main = "Stirling vs PIG-MCMC bias against gold standard")
  lines(ag_p$nu_alpha, ag_p$mean_abs_pp_diff, type = "b", lwd = 2,
        col = "#1B9E77", pch = 17)
  legend("topright",
         legend = c("gaga (Stirling)", "PIG-MCMC (slice)"),
         col = c("#D95F02", "#1B9E77"), pch = c(16, 17), lwd = 2, bty = "n")
  dev.off()
}
plot_alpha_sweep()

# ---------------- print headline -----------------------------------------
cat("\n=== HEADLINE TABLE ===\n")
print(headline, row.names = FALSE)
cat("\n=== HYPERPARAMETER RECOVERY (estimated vs oracle) ===\n")
print(hyper_recovery, row.names = FALSE)
cat("\n=== FDR / POWER / AUC (mean ± SE across reps) ===\n")
print(fdr_pwr_auc, row.names = FALSE)
cat("\n=== ALPHA SWEEP: gaga and PIG vs gold (mean over reps) ===\n")
sweep_pretty <- agg_sweep_cross[agg_sweep_cross$pair %in%
                                c("gaga_vs_gold", "pig_vs_gold"), ]
print(sweep_pretty[order(sweep_pretty$nu_alpha, sweep_pretty$pair), ],
      row.names = FALSE, digits = 4)
cat("\nFigures saved:\n")
cat("  ", af("scatter_pp_vs_gold.pdf"), "\n")
cat("  ", af("alpha_sweep_pp_diff.pdf"), "\n")
cat("  ", af("alpha_sweep_pp_diff.png"), "\n")
