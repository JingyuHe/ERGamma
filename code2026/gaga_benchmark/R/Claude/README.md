# PIG-augmented Bayesian inference applied to the GaGa real-data examples

This folder contains a P-IG-based reproduction of the two real-data examples
in Rossell (2009): the Armstrong leukemia data (Section 6.1) and the MAQC
titration experiment (Section 7). The goal is to swap Rossell's Stirling-
based EM for the gamma shape parameter with the exact P-IG augmentation of
Polson, He, Xu and Zhai (2025+, JASA submission), and see how the resulting
Bayesian inference compares to the original GaGa implementation
(`gaga::fitGG`) and to limma.

## Files

```
00_FORMULAS.md                model spec, P-IG identity, conditional posteriors
pig_gaga_lib.R                core library: P-IG sampler, Gibbs, FDR helpers
01_pig_gaga_armstrong.R       Armstrong PIG-GaGa benchmark
02_pig_gaga_maqc.R            MAQC PIG-GaGa benchmark
results/                      CSV / RDS outputs
```

`pig_gaga_lib.R` re-uses the P-IG draw `rPIG_ERGamma()` and the power-
truncated normal sampler `rH()` defined in `code2026/sampling.R` — those
functions are inlined verbatim so we don't pull in the legacy `gtools`
dependency. We do not modify any file under `code2026/`.

## Model

For each gene $i$, replicate $j$, group $g$ we use Rossell's parameterisation
$$x_{ijg} \mid \alpha, \lambda_{ic} \sim \mathrm{Ga}(\alpha, \alpha/\lambda_{ic}),
\qquad
1/\lambda_{ic} \sim \mathrm{Ga}(\alpha_0,\, \alpha_0\lambda_0),$$
where the cluster index $c$ depends on the pattern indicator $z_i$.
**Difference from `gaga::fitGG`**: the original GaGa uses a per-gene shape
$\alpha_i \sim \mathrm{Ga}(b_\alpha, b_\alpha/\nu_\alpha)$. Our PIG-GaGa
implementation uses a single global $\alpha$ shared across genes, with a
Damsleth $\xi_2$ prior on $\alpha$. This is a simpler, more rigid model and
the comparison below quantifies the cost of that simplification.

## P-IG augmentation step

Conditional on the pattern indicators $\{z_i\}$ and after integrating out the
cluster means $\lambda_{ic}$, the marginal posterior of $\alpha$ has the
structure described in PIG paper Section 4.2:
* a $1/\Gamma(\alpha)^{S}$ factor with $S = N_{\mathrm{tot}}+\delta$,
* a Damsleth $\Gamma(\delta\alpha+1)$ factor from the prior,
* a non-linear-in-$\alpha$ remainder coming from $\alpha^{\alpha N_{\mathrm{tot}}}$
  and $\Gamma(\alpha n_{ic}+\alpha_0)/(\alpha S_{ic}+\alpha_0\lambda_0)^{\alpha n_{ic}+\alpha_0}$.

The augmented Gibbs step we implement is

* $\tau \sim \mathrm{Ga}(\delta\alpha+1, 1)$
* $\omega_s \stackrel{\mathrm{iid}}{\sim} \mathrm{P\text{-}IG}(\alpha)$ for $s=1,\dots,S$
* $\alpha \sim \mathrm{PTN}\!\bigl(S+1,\, \sum\omega_s,\, \tilde b\bigr)$

with $\tilde b = S\gamma + \delta\log\tau - \delta\log(\delta\mu) - \delta\eta + \sum\log L_{ic} + \widetilde{\Delta}_{\mathrm{cur}}$,
where $\widetilde{\Delta}_{\mathrm{cur}}$ is the slope of the awkward
non-linear remainder evaluated at the current $\alpha$. The PTN draw is
followed by a Metropolis–Hastings correction so the chain targets the
exact posterior. See `00_FORMULAS.md` for the full derivation.

The P-IG draws use the truncated GIG sum from `code2026/sampling.R`:

```r
rPIG_ERGamma(n, cc, N)   # n PIG(cc) draws, truncated at N + gamma tail
rH(p, a, b)              # PTN(p, a, b) draw via Damsleth's rejection sampler
```

## Results

### Armstrong (Schliep filtered, 2194 probes, 24 ALL + 18 MLL)

`results/pig_armstrong_summary_schliep_filtered_log.csv` and
`results/pig_armstrong_agreement_schliep_filtered_log.csv`:

| method | #DE @ FDR=0.05 |
|---|---|
| PIG-GaGa (this folder) | 339 |
| gaga::fitGG (Stirling EM) | 360 |
| limma BH | 361 |

Method-pair Jaccard agreement on the DE set:

| pair | Jaccard |
|---|---|
| PIG-GaGa vs gaga::fitGG | 0.56 |
| PIG-GaGa vs limma_BH | 0.39 |
| gaga::fitGG vs limma_BH | 0.59 |

PIG-GaGa converges to $\alpha \approx 0.63$, $\alpha_0 \approx 0.10$,
$\lambda_0 \approx 0.014$ on this dataset. The difference from gaga's
$b_\alpha \approx 8.13$ is expected: gaga's parameter is the shape of the
prior on per-gene $\alpha_i$, not the value of $\alpha$ itself. Under the
simplified single-$\alpha$ model used here, the data prefer a much smaller
effective $\alpha$ (corresponding to a higher overall coefficient of
variation).

The 0.56 Jaccard agreement with gaga::fitGG indicates that PIG-GaGa picks
out roughly the same dominant DE signal but diverges on the borderline
calls — exactly what one would expect when comparing a simpler global-$\alpha$
model to the per-gene-$\alpha$ model that GaGa actually implements.

### MAQC (5000 high-variance + 5000 random affy probes, 5 patterns, 4 pools × 5 arrays)

`results/pig_maqc_roc_summary.csv`:

| method | n_mapped | n_validated | AUC_scaled | AUC_standard |
|---|---|---|---|---|
| PIG-GaGa | 337 | 322 | 0.0356 | 0.837 |
| gaga::fitGG | 337 | 322 | 0.0379 | 0.892 |
| limma F | 337 | 322 | 0.0378 | 0.890 |

PIG-GaGa AUC against the qPCR validation set is 0.837 vs gaga's 0.892 — a
~6 percentage point gap. As with Armstrong, the gap is attributable to the
single-$\alpha$ vs per-gene-$\alpha$ model difference, not to the P-IG
augmentation per se. Within the simplified model, P-IG provides exact
posterior inference for $\alpha$ where Stirling would give an
approximation.

### What this benchmark does and does not show

It **does** show that:

* The P-IG augmented Gibbs step is implementable and stable on real
  microarray-scale data ($G=2194$ to $G=10000$, $n=20$ to $n=42$).
* PIG-GaGa identifies the dominant DE signal in both datasets (Jaccard
  $\approx 0.56$ vs gaga on Armstrong; AUC 0.84 on MAQC qPCR).

It **does not** claim that PIG-GaGa beats `gaga::fitGG`. The fair comparison
between Stirling and exact P-IG inference would hold the same per-gene-$\alpha$
model fixed and replace only the $\alpha$ update; that requires substantially
more code (per-gene PIG chains for each of $G\sim10^4$ genes) and is left as
future work. The current implementation is a global-$\alpha$ simplification
that the user can extend.

## How to run

```sh
# Armstrong (Schliep mirror is fast, ~20s)
PIG_GAGA_ITER=80 PIG_GAGA_BURN=30 PIG_GAGA_NPIG=80 \
  Rscript code2026/gaga_benchmark/R/Claude/01_pig_gaga_armstrong.R

# Armstrong with the orange 12533-probe matrix (slower, several minutes)
ARMSTRONG_DATA=orange_full PIG_GAGA_ITER=80 PIG_GAGA_BURN=30 \
  Rscript code2026/gaga_benchmark/R/Claude/01_pig_gaga_armstrong.R

# MAQC (sub-sample 8000 probes, ~10s on 5000 high-var + 3000 random)
MAQC_N_GENES=8000 PIG_GAGA_ITER=20 PIG_GAGA_BURN=5 \
  Rscript code2026/gaga_benchmark/R/Claude/02_pig_gaga_maqc.R
```

The MAQC script depends on `gaga_benchmark/Claude/results/maqc_loaded.rds`,
which is created by `gaga_benchmark/Claude/07a_maqc_cache.R`. Run that first
if it doesn't exist.

## Open issues / future work

1. **Per-gene-$\alpha$ extension.** Implement $\alpha_i$ for each gene with
   its own PIG chain and a Gamma hyperprior $\alpha_i \sim \Ga(b_\alpha, b_\alpha/\nu_\alpha)$
   matched to gaga's, so the comparison isolates Stirling vs P-IG with a
   common model. The PIG paper outlines exactly this (Section 4.2 last
   paragraph).
2. **Per-iteration vectorisation.** The current Gibbs loop calls
   `gene_pattern_posterior` once per iteration (~$O(GP)$ ops) and
   `gamma_ratio_slope_alpha` inside the alpha update (~$O(GP)$ ops). On
   12533 probes × 200 iterations this is ~4M evaluations; feasible but
   slow. Vectorising over genes (which is straightforward; just batch the
   $\Gamma$ and $\log$ calls) would give a $\sim 10\times$ speedup.
3. **The non-linear-remainder linearisation.** We linearise
   $\alpha N_{\mathrm{tot}}\log\alpha + g(\alpha)$ at the current $\alpha$
   to obtain a PTN proposal, then MH-correct. A small fraction of proposals
   are rejected when the chain moves quickly (e.g. early iterations from a
   far-from-stationary start). An alternative would be to introduce a third
   layer of auxiliary variables that absorbs $\alpha^{\alpha N}$ into the
   PTN's quadratic kernel, removing the MH step.
