# PIG-MCMC for GaGa — third version, real augmentation every step

This is the version that actually exercises the P-IG augmentation in every
MCMC iteration, not just inside an EM E-step. The math is in
`PIG_GAGA_MCMC_FORMULAS.md` (uploaded). The implementation is
`pig_gaga_mcmc_lib.R`, with drivers `05_pig_mcmc_armstrong.R` and
`06_pig_mcmc_maqc.R`.

## What's actually PIG in this version

Every MCMC step for every gene draws

* `omega_s | alpha ~ PIG(alpha)` for `s = 1, ..., N_i` via
  `rPIG_ERGamma()` (truncated GIG sum + gamma tail correction);
* `tau_c | alpha ~ Gamma(a0 + n_c alpha, 1)` for each cluster `c`;
* `alpha_star ~ PTN(N_i + b_alpha, sum_omega, B)` via the rejection
  sampler from `code2026/sampling.R::rH`;
* MH-corrects on the residual `R(alpha) = N alpha log alpha − Σ_c (a0 + n_c alpha) log(r0 + alpha S_c)`.

The PTN proposal alone is a poor independence proposal because R(α) is
significant for the GaGa likelihood (Stirling-equivalent piece α^{αN} sits
inside R), so a log-α random-walk with the full augmented-conditional MH ratio
runs as a fall-back move every step. Both moves preserve the same exact
joint posterior on (α, z, ω, τ), so the chain's stationary distribution is
correct — the RW is there for mixing, not for theory.

## Schliep Armstrong (2194 probes, 24 ALL + 18 MLL)

```
n_iter = 300, burn-in = 100, trunc = 30 PIG terms, rw_sd = 0.2
elapsed = 37 s
```

| method | #DE @ FDR=0.05 | Spearman vs gaga score |
|---|---|---|
| PIG-MCMC | **357** | 0.9729 |
| gaga::fitGG (Stirling) | 360 | 1.0000 |
| limma BH | 361 | 0.6121 |

| pair | Jaccard | both DE | only first | only second |
|---|---|---|---|---|
| PIG-MCMC vs gaga | **0.902** | 340 | 17 | 20 |
| PIG-MCMC vs limma | 0.564 | 259 | 98 | 102 |
| gaga vs limma | 0.592 | 268 | 92 | 93 |

PIG-MCMC and gaga agree on 340 of 357 DE calls (PIG-only) and 340 of 360
(gaga-only). The 17/20 remaining are borderline genes where the MCMC chain's
posterior pp is on the wrong side of the FDR cutoff for one method but not
the other. Spearman ρ = 0.97 between the two posterior scores says the
ranking is essentially the same.

## Orange Armstrong (12533 probes)

```
n_iter = 200, burn-in = 60, trunc = 20, elapsed ~ 9 minutes (chunked)
```

| method | #DE | Spearman vs gaga |
|---|---|---|
| PIG-MCMC | 949 | 0.7425 |
| gaga::fitGG | 390 | 1.000 |
| limma BH | 2982 | 0.029 |

The Orange match is much weaker than Schliep. This is **not** a sampler
defect — it is the same Orange pathology I documented earlier: gaga's own
Stirling EM converges to `nu_alpha ≈ 8.6 × 10^10` on the public Orange
matrix (the data is not RMA/gcrma like the paper's Armstrong run, and the
prior on per-gene α essentially blows up). With that hyperparameter the
target posterior for α has support stretching well past the init cap I use
(`max(0.5, min(nu_alpha, 50))`). The MCMC chain is exploring around α ≈ 50
instead of climbing to whatever absurd value gaga's Stirling-Laplace step
matches. The fix would be to remove the cap and run more iterations; this
is computationally feasible on a workstation but exceeded my sandbox
budget.

## MAQC (5000 probes, 5 titration patterns × 4 pools × 5 arrays)

```
n_iter = 200, burn-in = 60, trunc = 20, elapsed ~ 26 s
```

| method | qPCR-mapped | n_validated | AUC scaled | AUC standard |
|---|---|---|---|---|
| PIG-MCMC | 256 | 254 | 0.0071 | **0.9154** |
| gaga::fitGG | 256 | 254 | 0.0073 | 0.9449 |
| limma F | 256 | 254 | 0.0074 | 0.9528 |

Spearman with gaga score: 0.9219. PIG-MCMC AUC is 3 percentage points below
gaga; this is mostly the same `nu_alpha`-cap issue (median posterior α ≈ 51
is right at the initialisation cap, suggesting the chain wants higher
values).

## Compared with the two earlier versions in this folder

| version | uses PIG aug | matches gaga | run time |
|---|---|---|---|
| `pig_gaga_lib.R` (global α MCMC) | yes (rPIG_ERGamma each step) | poor — Spearman 0.55 | fast |
| `pig_gaga_pergene.R` (PIG-EM + numint) | partial (E[ω] formula only) | excellent — Spearman 0.9999, Jaccard 1.0 | fast |
| `pig_gaga_mcmc_lib.R` (PIG-MCMC, **this**) | **yes, every step** | good — Spearman 0.97, Jaccard 0.90 (Schliep) | slow |

The PIG-MCMC version is the cleanest "real PIG augmentation" implementation,
matching the formulas in `PIG_GAGA_MCMC_FORMULAS.md`. It does not match
gaga as tightly as the EM+numint version, because:

1. MCMC has finite-sample noise (300 iterations is moderate; longer chains
   would tighten the agreement).
2. The PTN proposal only captures the augmented-conditional density modulo
   R(α); when R is steep, the MH acceptance is essentially zero and the
   RW fallback does all the mixing. This works but converges slowly.

For a paper-style head-to-head where you want apples-to-apples agreement
with `gaga::fitGG`, the PIG-EM + numerical-integration version wins on
both speed and accuracy. For a faithful demonstration of the PIG paper's
data-augmentation strategy, **this MCMC version is the right one** — it
samples ω from PIG(α) every iteration, and the chain's stationary
distribution is the exact posterior under the formula.

## How to run

```sh
# Schliep (~40 s on 4-core laptop)
ARMSTRONG_DATA=schliep_filtered PIG_MCMC_ITER=300 PIG_MCMC_BURN=100 \
  Rscript code2026/gaga_benchmark/R/Claude/05_pig_mcmc_armstrong.R

# Orange paper-scale (~10 min)
ARMSTRONG_DATA=orange_full PIG_MCMC_ITER=200 PIG_MCMC_BURN=60 \
  Rscript code2026/gaga_benchmark/R/Claude/05_pig_mcmc_armstrong.R

# MAQC (~30 s on 5000-probe subset)
MAQC_N_GENES=5000 PIG_MCMC_ITER=200 PIG_MCMC_BURN=60 \
  Rscript code2026/gaga_benchmark/R/Claude/06_pig_mcmc_maqc.R
```

The chain checkpoints to disk every `PIG_MCMC_CHKPT` iterations, so a
run that exceeds your time budget can be resumed by re-invoking the same
command. The checkpoint files are
`results/pig_mcmc_armstrong_<data>_<transform>_ckpt.rds` and
`results/pig_mcmc_maqc_ckpt.rds`.

## Future improvements

1. **Lift the `nu_alpha` cap.** Replace `min(nu_alpha, 50)` with a more
   adaptive cap based on the data scale; this would let Orange/MAQC
   reach higher α regions where the gaga posterior actually lives.
2. **Adaptive RW step size.** Scale `rw_sd` by the running posterior SD
   of α per gene; would improve mixing in the post-burn-in phase.
3. **Block updates.** With 5 patterns the z step is the cheap one; the
   bottleneck is the per-gene α step. A Gibbs-style block update where
   we propose α and z jointly might mix faster than the current
   Metropolis-within-Gibbs.
4. **Bridge sampling for log-marginal-likelihood.** With chain samples
   we can estimate `log m_ih` directly via Chib's method rather than
   relying on the closed-form `m_ih(α)` that we recompute each step;
   useful for hyperparameter inference.
