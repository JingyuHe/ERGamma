# GaGa Reproduction and P-IG Benchmark

This directory is a clean benchmark area. It does not modify the original
project scripts. All generated outputs go to `results/`, and optional packages
are installed into `r-lib/`.

## What the Rossell GaGa paper reports

Paper: Rossell (2009), *GaGa: A parsimonious and flexible model for
differential expression analysis*, Annals of Applied Statistics 3(3), 1035-1051.

The paper has three reproducible empirical/simulation targets:

1. Armstrong leukemia data:
   - ALL vs MLL, excluding the two later U95Av2 MLL arrays.
   - 20 random subset repeats with 5, 10, and 15 arrays per group.
   - Table 1 reports average DE calls and reproducibility against the full data.

2. Armstrong-based simulations:
   - 200 parametric simulations from a GaGa posterior predictive fit.
   - 200 nonparametric bootstrap simulations from Armstrong arrays.
   - Table 2 reports FDR and power for 5, 10, 15, and 20 arrays per group.
   - Figure 2(a) reports ROC/AUC on the nonparametric simulations.

3. MAQC:
   - First-site Affymetrix hgu133plus2 arrays, four titration pools A/B/C/D.
   - Five patterns: all equal; A differs; C differs; D differs; all differ.
   - Figure 2(b) uses qPCR-confirmed genes as the validation set.

The most relevant benchmark for the P-IG paper is the gamma-shape step used by
GaGa. The `gaga` package approximates this distribution through `rcgamma()` and
`mcgamma()`. The scripts below compare that approximation with exact numerical
moments, and separately compare the existing P-IG sampler with Miller-style
gamma approximation on the special gamma-shape examples already used in this
project.

## Scripts

Run from the repository root:

```sh
Rscript code2026/gaga_benchmark/R/install_deps.R
Rscript code2026/gaga_benchmark/R/01_vignette_smoke.R
Rscript code2026/gaga_benchmark/R/02_gashape_benchmark.R
Rscript code2026/gaga_benchmark/R/03_armstrong_style_simulation.R
Rscript code2026/gaga_benchmark/R/04_special_pig_vs_miller.R
Rscript code2026/gaga_benchmark/R/05_gashape_stress_grid.R
Rscript code2026/gaga_benchmark/R/06_armstrong_table1.R
Rscript code2026/gaga_benchmark/R/07_maqc_reproduce.R
Rscript code2026/gaga_benchmark/R/08_armstrong_gcrma_preprocess.R
Rscript code2026/gaga_benchmark/R/09_armstrong_figure1.R
Rscript code2026/gaga_benchmark/R/10_armstrong_simulation62.R
```

For a quick check, the defaults are intentionally small. For a paper-scale run:

```sh
GAGA_REPS=200 GAGA_N_GENES=2000 Rscript code2026/gaga_benchmark/R/03_armstrong_style_simulation.R
GAGA_GAS_CASES=50 GAGA_GAS_DRAWS=20000 Rscript code2026/gaga_benchmark/R/02_gashape_benchmark.R
PIG_DRAWS=10000 PIG_BURNIN=5000 PIG_TRUNC=1000 Rscript code2026/gaga_benchmark/R/04_special_pig_vs_miller.R
GAGA_STRESS_CASES=500 Rscript code2026/gaga_benchmark/R/05_gashape_stress_grid.R
ARMSTRONG_REPS=20 ARMSTRONG_N_GENES=0 Rscript code2026/gaga_benchmark/R/06_armstrong_table1.R
ARM62_REPS=200 Rscript code2026/gaga_benchmark/R/10_armstrong_simulation62.R
```

## Outputs

- `results/vignette_smoke_summary.csv`
- `results/gas_approx_benchmark.csv`
- `results/armstrong_style_simulation.csv`
- `results/armstrong_style_summary.csv`
- `results/special_pig_vs_miller.csv`
- `results/gas_stress_grid.csv`
- `results/gas_stress_worst.csv`
- `results/armstrong_table1_summary.csv`
- `results/armstrong_table1_replicates.csv`
- `results/armstrong_table1_data_summary.csv`
- `results/maqc_data_summary.csv`
- `results/maqc_pattern_counts.csv`
- `results/maqc_roc_summary.csv`
- `results/maqc_roc_curves.csv`
- `results/maqc_qpcr_scores.csv`
- `results/maqc_top_pattern4_genes.csv`
- `results/maqc_roc_curves.pdf`
- `results/maqc_roc_curves.png`
- `results/armstrong_figure1_reproduction_<data>_<transform>.pdf`
- `results/armstrong_figure1_reproduction_<data>_<transform>.png`
- `results/armstrong_figure1_summary_<data>_<transform>.csv`
- `results/armstrong62_table2_summary_<data>_<transform>_reps<N>.csv`
- `results/armstrong62_roc_summary_<data>_<transform>_reps<N>.csv`
- `results/armstrong62_roc_curves_<data>_<transform>_reps<N>.csv`
- `results/armstrong62_roc_curves_<data>_<transform>_reps<N>.pdf`
- `results/armstrong62_roc_curves_<data>_<transform>_reps<N>.png`
- `results/armstrong62_simulation_<data>_<transform>_reps<N>.rds`

The Armstrong script supports two public mirrors:

- `orange_full`: 12533 probes and 72 samples from the Orange cancer
  projection mirror. This is the default because it keeps the full feature set,
  but it is not the exact `just.gcrma` matrix used by Rossell. Values include
  negatives, so the default transform is `log2(pmax(x, 1) + 1)`.
- `schliep_filtered`: 2194 positive probes from Schliep lab's filtered
  CompCancer mirror.

The Rossell paper uses background-corrected, normalized and summarized CEL files
via `just.gcrma`. The public mirror reproduction is therefore a procedural
replication of Table 1, not bitwise reproduction of the published counts.

`09_armstrong_figure1.R` reproduces Figure 1 procedurally. Its default uses the
Schliep filtered Armstrong matrix with a natural-log transform because that
public matrix is positive and has the same approximate 0-10 log-scale range as
the published figure. The Orange/Broad full matrix can be plotted with
`FIG1_DATA=orange_full FIG1_TRANSFORM=log2_floor20`, but its scale is shifted
because it is not the `just.gcrma` matrix used in the paper.

`10_armstrong_simulation62.R` reproduces the Section 6.2 simulation workflow:
parametric simulations are drawn with `gaga::simGG()` using full-data GaGa
posterior-mean hyperparameters, and nonparametric bootstraps use full-data GaGa
calls as truth. Both are evaluated by GaGa, MiGaGa2, EBarrays Ga, and limma BH.
The default is `ARM62_REPS=20` for a practical local check. Use `ARM62_REPS=200`
for the paper-scale Table 2 and Figure 2(a) run.

The Broad publication page archived by Wayback identifies the original raw CEL
archives as:

- `mll_scans_ALL1.tar.gz`
- `mll_scans_ALL2.tar.gz`
- `mll_scans_MLL1.tar.gz`
- `mll_scans_MLL2.tar.gz`
- `mll_scans_AML1.tar.gz`
- `mll_scans_AML2.tar.gz`

The archived page and metadata files are available, but the tar payloads were
not archived by Wayback and the current Broad URLs return 404. If the CEL
archives are obtained separately, place them under
`data/armstrong_broad_raw/` and run `08_armstrong_gcrma_preprocess.R`.

## MAQC Notes

The MAQC script uses the public GEO processed matrices:

- `GSE5350-GPL570_series_matrix.txt.gz`: Affymetrix hgu133plus2 expression
  matrix. The reproduction keeps the first-site arrays `MAQC_AFX_1_A1` through
  `MAQC_AFX_1_D5`.
- `GSE5350-GPL4097_TaqMan_series_matrix.txt.gz`: TaqMan qPCR validation
  matrix.

Rossell reports 1296 qPCR validated genes. The public GPL4097 matrix currently
contains 1044 assays, so the qPCR ROC reproduction is a close procedural
replication rather than a bitwise recreation of Figure 2(b).

The MAQC script now also fits the EBarrays Gamma-Gamma model used as the Ga
baseline where possible, writes a PNG copy of the ROC curve, and exports the
top Pattern 4 gene lists for GaGa and MiGaGa2 in
`results/maqc_top_pattern4_genes.csv`.
