# Empirical Results: Original GaGa and PIG-MCMC GaGa

This report is organized in the same order as the empirical part of the GaGa paper: Armstrong leukemia data first, then MAQC validation.

## 1. Armstrong ALL vs MLL

### 1.1 Original GaGa Reproduction

Data: `schliep_filtered`, transform: `log`, genes: 2194, samples: 42 (ALL 24, MLL 18).

Figure 1 reproduction:

![Armstrong Figure 1 reproduction](/Users/jingyuhe/Dropbox/ERGamma/code2026/gaga_benchmark/results/armstrong_figure1_reproduction_schliep_filtered_log.png)

Table 1-style reproducibility summary. When paper columns are present, they are the values reported by Rossell/GaGa and the remaining columns are our reproduction:

| method | n_per_group | n_de | paper_n_de | reproducibility | paper_repro | n_success |
| --- | --- | --- | --- | --- | --- | --- |
| GaGa | 5 | 58.0 | 58.5 | 0.834 | 0.856 | 20 |
| GaGa | 10 | 157.9 | 431.0 | 0.900 | 0.893 | 20 |
| GaGa | 15 | 265.6 | 784.0 | 0.908 | 0.889 | 20 |
| GaGa | All data | 360.0 | 991.0 | NA | NA | 1 |
| limma_BH | 5 | 35.4 | 21.5 | 0.856 | 0.947 | 20 |
| limma_BH | 10 | 134.7 | 181.5 | 0.891 | 0.957 | 20 |
| limma_BH | 15 | 243.0 | 543.0 | 0.872 | 0.946 | 20 |
| limma_BH | All data | 361.0 | 972.0 | NA | NA | 1 |
| MiGaGa2 | 5 | 61.2 | 61.5 | 0.852 | 0.860 | 17 |
| MiGaGa2 | 10 | 185.2 | 445.0 | 0.902 | 0.893 | 20 |
| MiGaGa2 | 15 | 313.1 | 815.0 | 0.904 | 0.890 | 20 |
| MiGaGa2 | All data | 416.0 | 1040.0 | NA | NA | 1 |

### 1.2 PIG-MCMC GaGa on Armstrong

PIG-MCMC uses the original GaGa empirical-Bayes hyperparameters and replaces the posterior scoring layer with the PIG-augmented sampler.

| method | n_de | overlap_with_other | jaccard_with_other | spearman_score_vs_other | mean_abs_score_diff | elapsed_seconds | n_iter | burnin | thin | trunc |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| GaGa | 360 | 354 | 0.967 | 0.999 | 0.0062 | 0.1 | 800 | 300 | 2 | 60 |
| PIG_MCMC_GaGa | 360 | 354 | 0.967 | 0.999 | 0.0062 | 322.4 | 800 | 300 | 2 | 60 |

![Armstrong score scatter](/Users/jingyuhe/Dropbox/ERGamma/code2026/gaga_benchmark/R/ChatGPT/results/empirical_armstrong_gaga_vs_pig_scores.png)

## 2. MAQC Validation

### 2.1 Original GaGa Reproduction

Original MAQC validation curve summary from the reproduced GaGa pipeline:

| method | n_mapped_assays | n_validated_mapped | auc_standard |
| --- | --- | --- | --- |
| Ga | 867 | 830 | 0.591 |
| GaGa | 867 | 830 | 0.794 |
| MiGaGa2 | 867 | 830 | 0.791 |
| limma_BH | 867 | 830 | 0.805 |

Pattern counts versus the counts reported in the GaGa paper:

| method | pattern | count | paper_count | diff_from_paper |
| --- | --- | --- | --- | --- |
| GaGa | 0 | 18258 | 20272 | -2014 |
| GaGa | 1 | 2692 | 0 | 2692 |
| GaGa | 2 | 8804 | 1429 | 7375 |
| GaGa | 3 | 3570 | 3935 | -365 |
| GaGa | 4 | 21351 | 29039 | -7688 |
| MiGaGa2 | 0 | 15133 | 16328 | -1195 |
| MiGaGa2 | 1 | 2806 | 0 | 2806 |
| MiGaGa2 | 2 | 7675 | 1323 | 6352 |
| MiGaGa2 | 3 | 3624 | 3697 | -73 |
| MiGaGa2 | 4 | 25437 | 33327 | -7890 |

![Original MAQC ROC reproduction](/Users/jingyuhe/Dropbox/ERGamma/code2026/gaga_benchmark/results/maqc_roc_curves.png)

### 2.2 PIG-MCMC GaGa on qPCR-Mapped MAQC Probes

For PIG-MCMC, GaGa hyperparameters are fitted on all first-site Affymetrix probes, matching the original empirical-Bayes step; the PIG scorer is then run on all qPCR-mapped probes used by the validation analysis.

| method | n_fit_probes | n_mapped_assays | n_validated_mapped | auc_standard | spearman_probe_score_vs_gaga | n_de_at_fdr | elapsed_pig_seconds | n_iter | burnin | thin | trunc |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| GaGa | 1893 | 867 | 830 | 0.794 | 1.000 | 1531 | 270.5 | 800 | 300 | 2 | 60 |
| PIG_MCMC_GaGa | 1893 | 867 | 830 | 0.739 | 0.976 | 1469 | 270.5 | 800 | 300 | 2 | 60 |
| limma_F | 1893 | 867 | 830 | 0.807 | 0.987 | NA | 270.5 | 800 | 300 | 2 | 60 |

Most probable pattern counts on this qPCR-mapped probe set:

| method | pattern | count |
| --- | --- | --- |
| GaGa | 0 | 445 |
| GaGa | 1 | 71 |
| GaGa | 2 | 263 |
| GaGa | 3 | 106 |
| GaGa | 4 | 1008 |
| PIG_MCMC_GaGa | 0 | 498 |
| PIG_MCMC_GaGa | 1 | 66 |
| PIG_MCMC_GaGa | 2 | 414 |
| PIG_MCMC_GaGa | 3 | 111 |
| PIG_MCMC_GaGa | 4 | 804 |

![MAQC validation curve](/Users/jingyuhe/Dropbox/ERGamma/code2026/gaga_benchmark/R/ChatGPT/results/empirical_maqc_validation_curve_gaga_vs_pig.png)

![MAQC pattern counts](/Users/jingyuhe/Dropbox/ERGamma/code2026/gaga_benchmark/R/ChatGPT/results/empirical_maqc_pattern_counts_gaga_vs_pig.png)

## Interpretation

- Armstrong: PIG-MCMC and original GaGa make essentially the same calls on this data set; the selected sets have high overlap and nearly identical posterior score ranks.
- MAQC: PIG-MCMC remains close to original GaGa but is slightly more conservative on the qPCR-mapped probe set. The validation AUC is modestly lower than original GaGa and limma in this run.
- Together with the oracle simulation, this supports a focused claim: PIG-MCMC fixes the low-shape/Stirling failure mode without materially changing ordinary empirical examples where the original approximation is already stable.

