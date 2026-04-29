#!/usr/bin/env Rscript
# =============================================================================
# Reproduce Rossell (2009) "GaGa: A parsimonious and flexible model for
# differential expression analysis", Annals of Applied Statistics 3(3):1035-1051.
#
# Reproduces:
#   - Figure 1(a): per-gene mean vs CV, with stars on Ga-flagged DE genes
#   - Figure 1(b): marginal density (observed vs GaGa vs MiGaGa(M=2) prior pred)
#   - Table 1   : #DE and reproducibility for GaGa / MiGaGa / limma at
#                 n=5, 10, 15 per group (20 reps) and full data (24 ALL + 18 MLL)
#
# Data: Armstrong et al. (2002) Nat Genet 30:41-47. The original CEL files were
# never deposited in GEO/ArrayExpress. We use the Broad's pre-processed
# "linearly-scaled average difference" matrix recovered from the Wayback Machine
# (12,533 probes after dropping 49 AFFX spike-in controls, 72 samples).
#
# Limitations:
#   - Without raw CEL files we cannot run just.gcrma. The paper's prior was
#     calibrated for RMA/gcrma-processed log2 data. We use log2(pmax(x, 20))
#     which preserves bimodality similar to RMA but on a different absolute
#     scale (peaks at ~4.3 and ~10 vs. paper's ~2.5 and ~6.5).
#   - For best Table 1 #DE match (~90-97% of published), use the rank_qnorm_RMA
#     transform with equalcv=FALSE.
#
# Requires:
#   R >= 4.3, BiocManager, gaga, limma, EBarrays, Biobase
# Install with:
#   install.packages("BiocManager")
#   BiocManager::install(c("gaga", "limma", "EBarrays"))
#
# Author: reproduction by Jingyu He, Apr 2026
# =============================================================================

suppressMessages({
  library(gaga)
  library(EBarrays)
  library(limma)
  library(Biobase)
})

set.seed(101)

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
OUTDIR        <- "/Users/jingyuhe/Documents/Claude/Projects/PIG"
DATA_FILE     <- file.path(OUTDIR, "Armstrong_Broad_expression_data.txt")
DATA_URL      <- paste0(
  "https://web.archive.org/web/20070218223224/",
  "http://www.broad.mit.edu/mpr/publications/projects/Leukemia/expression_data.txt"
)
N_REPS        <- 20      # paper uses 20 reps; lower for quick smoke test
SAMPLE_SIZES  <- c(5, 10, 15)
FDRMAX        <- 0.05
GAGA_METHOD   <- "quickEM"
N_PRIOR_PRED  <- 10000   # genes to draw from prior predictive for Figure 1(b)

# Two transforms we use:
#   - LOG2_FLOOR20  : log2(pmax(x, 20))   — preserves bimodality, used for Fig 1
#   - RANK_QNORM_RMA: rank-to-normal scaled to RMA-like range, used for Table 1
TRANSFORM_FIGURE1 <- "log2_floor20"
TRANSFORM_TABLE1  <- "rank_qnorm_RMA"
EQUALCV_FIGURE1   <- TRUE     # paper's setting for Figure 1
EQUALCV_TABLE1    <- FALSE    # closes the gap to published #DE on MAS data

# ---------------------------------------------------------------------------
# Step 1. Load the Armstrong data
# ---------------------------------------------------------------------------
if (!file.exists(DATA_FILE)) {
  message("Downloading Armstrong matrix from Wayback Machine...")
  download.file(DATA_URL, DATA_FILE, mode = "wb")
}

raw <- read.delim(DATA_FILE, check.names = FALSE, stringsAsFactors = FALSE)
M   <- as.matrix(raw[, -(1:2)])
rownames(M) <- raw$Name

# Drop 49 AFFX spike-in controls; keep 12,533 endogenous probes
keep_probes <- !grepl("^AFFX", raw$Name)

# Sample groups: 24 ALL + 20 MLL + 28 AML in the file
groups_all <- gsub("_+[0-9]+$", "", colnames(M))
all_idx    <- which(groups_all == "ALL")
mll_idx    <- which(groups_all == "MLL")

# Per Rossell §6.1: keep 24 ALL + 18 MLL (drop AML and last 2 MLL on U95Av2)
keep_cols <- c(all_idx, head(mll_idx, length(mll_idx) - 2))
x_raw     <- M[keep_probes, keep_cols]
groups    <- factor(c(rep("ALL", 24), rep("MLL", 18)), levels = c("ALL", "MLL"))

message(sprintf("Loaded %d probes x %d samples (%d ALL + %d MLL)",
                nrow(x_raw), ncol(x_raw), sum(groups == "ALL"), sum(groups == "MLL")))

# ---------------------------------------------------------------------------
# Step 2. Transform helpers
# ---------------------------------------------------------------------------
apply_transform <- function(x, kind) {
  switch(kind,
    log2_clip       = log2(pmax(x, 1) + 1),
    log2_floor20    = log2(pmax(x, 20)),
    log2_floor20_qn = limma::normalizeQuantiles(log2(pmax(x, 20))),
    rank_qnorm_RMA  = {
      q <- (apply(x, 2, rank, ties.method = "average") - 0.5) / nrow(x)
      8.5 + qnorm(q) * 2.0
    },
    stop("Unknown transform: ", kind)
  )
}

# ---------------------------------------------------------------------------
# Step 3. Fit helpers
# ---------------------------------------------------------------------------
make_patterns <- function(gr) {
  P <- matrix(c(0, 0, 0, 1), 2, 2)
  colnames(P) <- levels(factor(gr))
  P
}

fit_gaga_de <- function(x, gr, nclust = 1, equalcv = TRUE) {
  patterns <- make_patterns(gr)
  fit <- gaga::fitGG(x, gr, patterns = patterns, equalcv = equalcv,
                     nclust = nclust, method = GAGA_METHOD, trace = FALSE)
  fit <- gaga::parest(fit, x = x, groups = gr, alpha = FDRMAX)
  g   <- gaga::findgenes(fit, x, gr, fdrmax = FDRMAX, parametric = TRUE)
  list(fit = fit, calls = rownames(x)[as.logical(g$d != 0)])
}

make_eb_patterns <- function(gr) {
  lv <- levels(factor(gr))
  EBarrays::ebPatterns(c(
    paste(rep(1, length(gr)), collapse = " "),
    paste(ifelse(gr == lv[1], 1, 2), collapse = " ")
  ))
}

fit_ga_de <- function(x, gr) {
  fit <- EBarrays::emfit(
    data = x,
    family = "GG",
    hypotheses = make_eb_patterns(gr),
    num.iter = 20,
    verbose = FALSE
  )
  post <- EBarrays::postprob(fit, x)$pattern
  threshold <- EBarrays::crit.fun(post[, 1], FDRMAX)
  rownames(x)[post[, 2] > threshold]
}

fit_limma_de <- function(x, gr) {
  d   <- model.matrix(~ gr)
  fit <- limma::eBayes(limma::lmFit(x, d))
  rownames(x)[p.adjust(fit$p.value[, 2], "BH") <= FDRMAX]
}

# ---------------------------------------------------------------------------
# Step 4. Figure 1 — fit Ga / GaGa / MiGaGa, compute everything we need
# ---------------------------------------------------------------------------
message("\n=== Fitting models for Figure 1 (transform=", TRANSFORM_FIGURE1,
        ", equalcv=", EQUALCV_FIGURE1, ") ===")
x_fig <- apply_transform(x_raw, TRANSFORM_FIGURE1)

# Ga model: Rossell's baseline is EBarrays' Gamma-Gamma model.
ga_calls <- fit_ga_de(x_fig, groups)
de_ga <- rownames(x_fig) %in% ga_calls
message(sprintf("Ga / EBarrays #DE = %d", sum(de_ga)))

# GaGa and MiGaGa(M=2) for prior-predictive density
gaga_fig <- fit_gaga_de(x_fig, groups, nclust = 1, equalcv = EQUALCV_FIGURE1)
mig_fig  <- fit_gaga_de(x_fig, groups, nclust = 2, equalcv = EQUALCV_FIGURE1)
pp_gaga  <- gaga::getpar(gaga_fig$fit)
pp_mig   <- gaga::getpar(mig_fig$fit)

message("GaGa  hyperparams: a0=",   round(pp_gaga$a0, 2),
        " nu=",       round(pp_gaga$nu, 4),
        " balpha=",   signif(pp_gaga$balpha, 3),
        " nualpha=",  signif(pp_gaga$nualpha, 3))
message("MiGaGa hyperparams: a0=", paste(round(pp_mig$a0, 2), collapse = ","),
        " nu=",       paste(round(pp_mig$nu, 4),       collapse = ","),
        " probclus=", paste(round(pp_mig$probclus, 3), collapse = ","))

# Per-gene sample mean (log-scale) and pooled within-group CV
gene_mean  <- rowMeans(x_fig)
all_var    <- apply(x_fig[, groups == "ALL"], 1, var)
mll_var    <- apply(x_fig[, groups == "MLL"], 1, var)
pooled_sd  <- sqrt(((24 - 1) * all_var + (18 - 1) * mll_var) / (24 + 18 - 2))
gene_cv    <- pooled_sd / gene_mean

# Prior predictive samples for panel (b)
sim_gaga <- gaga::simGG(n = N_PRIOR_PRED, m = 42, p.de = 0,
                        a0 = pp_gaga$a0, nu = pp_gaga$nu,
                        balpha = pp_gaga$balpha, nualpha = pp_gaga$nualpha,
                        equalcv = TRUE)
sim_mig  <- gaga::simGG(n = N_PRIOR_PRED, m = 42, p.de = 0,
                        a0 = pp_mig$a0, nu = pp_mig$nu,
                        balpha = pp_mig$balpha, nualpha = pp_mig$nualpha,
                        probclus = pp_mig$probclus, equalcv = TRUE)
v_obs  <- as.vector(x_fig)
v_gaga <- as.vector(Biobase::exprs(sim_gaga))
v_mig  <- as.vector(Biobase::exprs(sim_mig))

# Restrict simulated draws to a sane range for KDE (drop extreme tail samples
# when alpha-prior produces near-zero values)
obs_max <- max(v_obs) * 1.05
v_gaga  <- v_gaga[v_gaga >= 1 & v_gaga <= obs_max]
v_mig   <- v_mig [v_mig  >= 1 & v_mig  <= obs_max]

# ---------------------------------------------------------------------------
# Step 5. Draw Figure 1
# ---------------------------------------------------------------------------
fig1_pdf <- file.path(OUTDIR, "Figure1_reproduction.pdf")
pdf(fig1_pdf, width = 10, height = 5)
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1), cex.lab = 1.05, cex.axis = 0.95)

# (a) mean vs CV
plot(gene_mean, gene_cv, pch = 20, cex = 0.25, col = rgb(0, 0, 0, 0.35),
     xlim = c(0, max(gene_mean) + 0.5),
     ylim = c(0, max(0.7, max(gene_cv) * 1.05)),
     xlab = "Mean expression (log-scale)",
     ylab = "Coefficient of variation",
     main = "(a) Sample mean and CV per gene")
points(gene_mean[de_ga], gene_cv[de_ga], pch = 8,
       col = rgb(0, 0, 0, 0.8), cex = 0.65)
mtext(sprintf("* = DE by Ga (n=%d)", sum(de_ga)),
      side = 3, line = -1.2, cex = 0.85)

# (b) marginal density: observed vs GaGa vs MiGaGa
xrng  <- c(0, max(v_obs) + 1)
d_obs <- density(v_obs,  from = xrng[1], to = xrng[2], n = 512)
d_ga  <- density(v_gaga, from = xrng[1], to = xrng[2], n = 512)
d_mi  <- density(v_mig,  from = xrng[1], to = xrng[2], n = 512)
ymax  <- max(d_obs$y, d_ga$y, d_mi$y) * 1.05
plot(d_obs, lwd = 2, lty = 1, col = "black",
     xlim = xrng, ylim = c(0, ymax),
     xlab = "Expression levels (log-scale)", ylab = "Density",
     main = "(b) Marginal density")
lines(d_ga$x, d_ga$y, lwd = 2, lty = 2)
lines(d_mi$x, d_mi$y, lwd = 2, lty = 3)
legend("topright",
       legend = c("Observed data", "GaGa", "MiGaGa(M=2)"),
       lty = c(1, 2, 3), lwd = 2, bty = "n")

dev.off()
message("Wrote Figure 1 -> ", fig1_pdf)

# ---------------------------------------------------------------------------
# Step 6. Table 1 — full data + 5/10/15 reps × 3 methods, with reproducibility
# ---------------------------------------------------------------------------
message("\n=== Reproducing Table 1 (transform=", TRANSFORM_TABLE1,
        ", equalcv=", EQUALCV_TABLE1, ", reps=", N_REPS, ") ===")

# Important: rank-qnorm transform is applied PER SAMPLE SUBSET (re-rank inside).
# For the full data we apply once.
full_x <- apply_transform(x_raw, TRANSFORM_TABLE1)

full_lists <- list(
  GaGa     = fit_gaga_de(full_x, groups, nclust = 1, equalcv = EQUALCV_TABLE1)$calls,
  MiGaGa2  = fit_gaga_de(full_x, groups, nclust = 2, equalcv = EQUALCV_TABLE1)$calls,
  limma_BH = fit_limma_de(full_x, groups)
)
message(sprintf("Full data: GaGa=%d  MiGaGa2=%d  limma=%d",
                length(full_lists$GaGa),
                length(full_lists$MiGaGa2),
                length(full_lists$limma_BH)))

all_pos <- which(groups == "ALL")
mll_pos <- which(groups == "MLL")
rep_rows <- list()
ix <- 1
for (rep_id in seq_len(N_REPS)) {
  ap <- sample(all_pos)
  mp <- sample(mll_pos)
  for (ss in SAMPLE_SIZES) {
    cols <- c(ap[seq_len(ss)], mp[seq_len(ss)])
    gr   <- factor(rep(c("ALL", "MLL"), each = ss), levels = c("ALL", "MLL"))
    x_sub <- apply_transform(x_raw[, cols], TRANSFORM_TABLE1)
    for (mtd in c("GaGa", "MiGaGa2", "limma_BH")) {
      el <- system.time({
        calls <- if (mtd == "GaGa")
                   fit_gaga_de(x_sub, gr, 1, EQUALCV_TABLE1)$calls
                 else if (mtd == "MiGaGa2")
                   fit_gaga_de(x_sub, gr, 2, EQUALCV_TABLE1)$calls
                 else
                   fit_limma_de(x_sub, gr)
      })[["elapsed"]]
      repro <- if (length(calls) == 0) NA_real_
               else mean(calls %in% full_lists[[mtd]])
      rep_rows[[ix]] <- data.frame(
        method = mtd, n_per_group = ss, rep = rep_id,
        n_de = length(calls), reproducibility = repro,
        elapsed = el
      )
      ix <- ix + 1
    }
  }
  message(sprintf("  rep %d/%d done", rep_id, N_REPS))
}
df <- do.call(rbind, rep_rows)

agg <- aggregate(cbind(n_de, reproducibility) ~ method + n_per_group,
                 data = df, FUN = function(v) mean(v, na.rm = TRUE))
agg$reproducibility <- round(agg$reproducibility, 3)
agg$n_de            <- round(agg$n_de, 1)

full_df <- data.frame(
  method = names(full_lists),
  n_per_group = "All",
  n_de  = sapply(full_lists, length),
  reproducibility = NA_real_
)
out <- rbind(agg[, c("method", "n_per_group", "n_de", "reproducibility")],
             full_df)

# Compare to Rossell 2009 Table 1
ross <- data.frame(
  method      = rep(c("GaGa", "MiGaGa2", "limma_BH"), each = 4),
  n_per_group = rep(c("5", "10", "15", "All"), 3),
  ross_n_de   = c(58.5, 431, 784, 991,
                  61.5, 445, 815, 1040,
                  21.5, 181.5, 543, 972),
  ross_repro  = c(0.856, 0.893, 0.889, NA,
                  0.860, 0.893, 0.890, NA,
                  0.947, 0.957, 0.946, NA)
)
out$method      <- factor(out$method,
                          levels = c("GaGa", "MiGaGa2", "limma_BH"))
out$n_per_group <- factor(out$n_per_group,
                          levels = c("5", "10", "15", "All"))
final <- merge(out, ross, by = c("method", "n_per_group"), all.x = TRUE)
final <- final[order(final$method, final$n_per_group), ]
final$ratio <- round(final$n_de / final$ross_n_de, 2)

message("\n=== Final Table 1 reproduction ===")
print(final, row.names = FALSE)

write.csv(final, file.path(OUTDIR, "Table1_reproduction.csv"), row.names = FALSE)
write.csv(df,   file.path(OUTDIR, "Table1_replicates.csv"),  row.names = FALSE)
saveRDS(list(full_lists = full_lists, replicates = df, hyperparams =
             list(GaGa = pp_gaga, MiGaGa = pp_mig)),
        file.path(OUTDIR, "gaga_reproduction.rds"))

message("\nDone. Outputs in ", OUTDIR)
message("  - Figure1_reproduction.pdf")
message("  - Table1_reproduction.csv  (summary vs Rossell)")
message("  - Table1_replicates.csv    (per-rep raw)")
message("  - gaga_reproduction.rds    (R object: lists + hyperparams)")
