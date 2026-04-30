# Simulation design: GaGa Stirling vs PIG-MCMC posterior scorer

Goal: isolate the effect of replacing GaGa's Stirling/gamma-shape posterior
scoring with PIG-MCMC, while keeping the empirical-Bayes GaGa fit unchanged.

For each replicate:

```text
1. Simulate two-group data from gaga::simGG under the original GaGa model.
2. Fit the original empirical-Bayes GaGa model:
     theta_hat = gaga::fitGG(..., method = "EM")
     fit       = gaga::parest(fit, ...)
   This is the same plug-in hyperparameter step as empirical GaGa.
3. Original scorer:
     pp_stirling = fit$pp
     calls_stirling = gaga::findgenes(fit, ...)
4. PIG scorer:
     keep the same theta_hat, same patterns, same data;
     run pig_mcmc_gaga_pp() to sample alpha_i and z_i with PIG augmentation;
     pp_pig = retained Rao-Blackwellized mean Pr(z_i | alpha_i, x_i, theta_hat).
5. Compare both scorers against simGG truth and against each other.
```

Primary outputs:

```text
stirling_vs_pig_mcmc_sim_methods.csv
  replicate-level FDR, power, AUC, #DE, elapsed time for each method

stirling_vs_pig_mcmc_sim_pairwise.csv
  replicate-level Jaccard, Spearman(score), mean absolute score difference,
  GaGa-only calls and PIG-only calls

stirling_vs_pig_mcmc_sim_method_summary.csv
stirling_vs_pig_mcmc_sim_pair_summary.csv
  averages across replicates
```

Default command for a small check:

```sh
Rscript code2026/gaga_benchmark/R/ChatGPT/07_simulate_stirling_vs_pig_mcmc.R
```

Paper-scale command:

```sh
SIM_PIG_REPS=100 SIM_PIG_N_GENES=2000 SIM_PIG_N_PER_GROUP=10 \
SIM_PIG_GAGA_METHOD=EM PIG_MCMC_ITER=1000 PIG_MCMC_BURNIN=400 \
PIG_MCMC_TRUNC=100 \
Rscript code2026/gaga_benchmark/R/ChatGPT/07_simulate_stirling_vs_pig_mcmc.R
```

Interpretation:

This simulation does not test a new hyperparameter estimator.  It tests the
posterior scoring layer:

```text
same simulated data
same fitGG EM theta_hat
same posterior FDR rule
Stirling scorer vs PIG-MCMC scorer
```

The stronger oracle simulation is implemented in
`08_oracle_stirling_pig_simulation.R`.  It generates data directly from the
GaGa model under fixed oracle parameters, scores each replicate with
high-accuracy log-alpha quadrature, original `gaga` Stirling scoring, and
PIG-MCMC scoring, and then reports approximation error against the exact
quadrature posterior.  `09_summarize_oracle_paper_simulation.R` combines the
per-scenario paper-scale runs into manuscript-ready tables.
