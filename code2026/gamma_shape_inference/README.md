# Gamma Shape Inference Experiments

This directory contains the replacement experiment framework for the Gamma
shape inference section of the P-IG manuscript.  It is intentionally separate
from the older `gamma_shape*.R` scripts because those scripts execute global
code on source and mix demonstration plots with method definitions.

## Methods

The scripts compare the following targets and methods:

- Numerical grid: one-dimensional ground truth for `xi_2`.
- P-IG Gibbs: exact sampler under the current integer posterior `delta`
  construction.
- Slice-log-alpha: exact generic baseline on `log(alpha)`.
- Miller-Gamma: derivative-matching gamma approximation from Miller (2019).
- Miller-IMH: exact independence-MH sampler using Miller-Gamma as the proposal.
- Stirling-Gamma: Rossell-style gamma approximation from Stirling's formula.

The P-IG implementation marks non-integer posterior `delta` as skipped.  This is
deliberate: the current manuscript derivation introduces `delta` independent
P-IG variables and is exact only when the posterior `delta` is an integer.

## Quick Check

Run a small smoke test:

```sh
Rscript code2026/gamma_shape_inference/run_gamma_shape_experiments.R --suite smoke --profile smoke
Rscript code2026/gamma_shape_inference/plot_gamma_shape_results.R \
  --result-rds code2026/gamma_shape_inference/results/smoke_smoke/gamma_shape_results.rds
```

The smoke run is only for code validation.  Do not use it as paper evidence.

## Full Experiments

For full runs, prefer the checkpoint runner.  It writes `gamma_shape_raw.csv`
after every chunk and appends progress/ETA rows to `progress.csv`.

Run all non-P-IG competitors first:

```sh
Rscript code2026/gamma_shape_inference/run_gamma_shape_experiments_checkpoint.R \
  --suite all --profile full --cores 8 --chunk-size 128 \
  --run-pig false \
  --out-dir code2026/gamma_shape_inference/results/all_full_benchmarks_no_pig
```

Then run a small P-IG timing probe:

```sh
Rscript code2026/gamma_shape_inference/run_gamma_shape_experiments_checkpoint.R \
  --suite all --profile full --reps 2 --cores 8 --chunk-size 16 \
  --run-slice false --run-approximations false --run-pig true \
  --pig-max-delta 500 \
  --out-dir code2026/gamma_shape_inference/results/pig_timing_probe_reps2_fixed
```

Finally run the P-IG-only full experiment:

```sh
Rscript code2026/gamma_shape_inference/run_gamma_shape_experiments_checkpoint.R \
  --suite all --profile full --cores 8 --chunk-size 16 \
  --run-slice false --run-approximations false --run-pig true \
  --pig-max-delta 500 \
  --out-dir code2026/gamma_shape_inference/results/all_full_pig_only_fixed
```

The manuscript tables combine the completed non-P-IG and P-IG outputs:

```sh
Rscript code2026/gamma_shape_inference/make_gamma_shape_manuscript_tables.R \
  --result-dirs code2026/gamma_shape_inference/results/all_full_benchmarks_no_pig,code2026/gamma_shape_inference/results/all_full_pig_only_fixed
```

The legacy runner remains useful for small focused checks:

Damsleth replication:

```sh
Rscript code2026/gamma_shape_inference/run_gamma_shape_experiments.R \
  --suite damsleth --profile full --reps 100 --cores 4
```

Main cross-regime simulation:

```sh
Rscript code2026/gamma_shape_inference/run_gamma_shape_experiments.R \
  --suite main --profile full --reps 200 --cores 4
```

Prior robustness:

```sh
Rscript code2026/gamma_shape_inference/run_gamma_shape_experiments.R \
  --suite prior --profile full --reps 100 --cores 4
```

P-IG truncation sensitivity:

```sh
Rscript code2026/gamma_shape_inference/run_gamma_shape_experiments.R \
  --suite truncation --profile full --reps 50 --cores 4
```

Extreme low-shape stress test:

```sh
Rscript code2026/gamma_shape_inference/run_gamma_shape_experiments.R \
  --suite extreme --profile full --reps 20 --cores 4
```

## Large `n` P-IG Runs

Exact P-IG Gibbs draws `delta` auxiliary P-IG variables at each iteration.  For
`n=500`, the default guard skips P-IG unless explicitly overridden:

```sh
Rscript code2026/gamma_shape_inference/run_gamma_shape_experiments.R \
  --suite main --profile full --reps 200 --cores 4 --pig-max-delta 500
```

This is computationally expensive and should be run separately from the generic
baselines if wall time matters.

## Outputs

Each run writes:

- `gamma_shape_raw.csv`: one row per task and method.
- `gamma_shape_summary.csv`: cell-level means, Monte Carlo SEs, coverage, and
  runtime summaries.
- `gamma_shape_results.rds`: raw results, summary, settings, and task metadata.

The plotting script writes figures under the run directory's `figures/`
subdirectory.
