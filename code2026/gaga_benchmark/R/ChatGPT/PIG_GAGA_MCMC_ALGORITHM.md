# PIG-MCMC GaGa benchmark algorithm

This is the runnable algorithm implemented in `pig_mcmc_gaga_lib.R`.

## Inputs

```text
X         G x n positive expression matrix
group     sample group labels
patterns  H x K pattern matrix; row 1 is the all-equal/null pattern
fit       gaga::fitGG object after parest(), used only for hyperparameters
```

The benchmark fixes

```text
theta = (a0, nu, b_alpha, nu_alpha, pi)
```

at the `gaga::fitGG` estimates so the comparison isolates the posterior
shape/pattern computation.

## Per-gene sampler

For each gene `i`:

```text
Initialize alpha_i by method of moments.
Initialize z_i from Pr(z_i | alpha_i, X_i, theta).

For iter = 1, ..., n_iter:
  1. Compute collapsed pattern log-probabilities
       log pi_h + log m_ih(alpha_i),  h = 0, ..., H-1.

  2. Convert them to p_h = Pr(z_i = h | alpha_i, X_i, theta).

  3. Sample z_i ~ Categorical(p_0, ..., p_{H-1}).

  4. For the selected pattern, collect cluster stats
       n_c, S_c = sum x, L_c = sum log x.

  5. Draw PIG auxiliaries using the original `sampling.R` implementation
     with finite GIG terms plus the Gamma tail correction
       omega_s | alpha_i ~ PIG(alpha_i), s = 1, ..., sum_c n_c.

  6. Draw Damsleth auxiliaries
       tau_c | alpha_i ~ Gamma(a0 + n_c alpha_i, 1).

  7. Define y = log alpha and sample y by stepping-out slice sampling from
       ell(y) = (N + b_alpha)y - Omega exp(2y) + B exp(y) + R(exp(y)).

  8. Set alpha_i = exp(y).

  9. After burn-in, accumulate the Rao-Blackwellized pattern probability p_h.

Return
  pp_ih = retained average of p_h,
  z_counts for diagnostics,
  alpha_mean and alpha movement rate.
```

The default code path does not call numerical integration or Laplace
approximation.  A plain exact PTN independence kernel is kept as
`PIG_MCMC_ALPHA_KERNEL=ptn` for diagnostics, but the benchmark default is
`PIG_MCMC_ALPHA_KERNEL=slice` because the PTN independence kernel has near-zero
movement on the real data.

## Benchmark scripts

```text
04_armstrong_pig_mcmc_gaga.R
  Armstrong ALL vs MLL, reports DE counts, overlap/Jaccard, Spearman score
  correlation, and PIG acceptance diagnostics.

05_maqc_pig_mcmc_gaga.R
  MAQC first-site Affymetrix vs qPCR validation, reports ROC/AUC summaries,
  pattern counts, probe scores, and PIG acceptance diagnostics.

06_run_pig_mcmc_all.R
  Runs both entry scripts with the current environment settings.
```

## Commands used for the current fast benchmark

```sh
PIG_MCMC_ITER=150 PIG_MCMC_BURNIN=50 PIG_MCMC_TRUNC=40 \
  PIG_MCMC_ARM_PROGRESS_EVERY=500 \
  Rscript code2026/gaga_benchmark/R/ChatGPT/04_armstrong_pig_mcmc_gaga.R

PIG_MCMC_MAQC_MAX_PROBES=1000 PIG_MCMC_ITER=120 PIG_MCMC_BURNIN=40 \
  PIG_MCMC_TRUNC=30 PIG_MCMC_MAQC_PROGRESS_EVERY=250 \
  Rscript code2026/gaga_benchmark/R/ChatGPT/05_maqc_pig_mcmc_gaga.R
```

These are intentionally medium-length runs.  For final paper numbers, increase
`PIG_MCMC_ITER`, `PIG_MCMC_BURNIN`, and `PIG_MCMC_TRUNC`, and run MAQC with
`PIG_MCMC_MAQC_MAX_PROBES=0` to use all qPCR-mapped probes.
