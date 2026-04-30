# Gamma Shape Inference Simulation Results

Generated at: 2026-04-30 19:24:26

## Scope

- Raw result rows: 39,280.
- Successful result rows: 38,680.
- Skipped result rows: 600.
- The numerical grid is used as the one-dimensional ground truth, not as a sampling competitor.
- P-IG Gibbs is run only where the posterior delta is integer under the current construction.

## Main Cross-Regime Simulation

- P-IG Gibbs mean KS over the 18 main cells is 0.0433; median ESS/sec is 81.8.
- Miller-IMH mean KS over the main cells is 0.0137; median ESS/sec is 46916.9.
- Among deterministic approximations, Miller-Gamma has mean KS 0.0024 versus 0.3358 for Stirling-Gamma.
- P-IG Gibbs is the lowest-KS exact sampler in 0 of 18 main cells, and its median ESS/sec ratio to the fastest exact competitor is 0.0017.

## Extreme Low-Shape Stress Test

- P-IG Gibbs: mean KS 0.0141, mean W1 0.0002, coverage 1.000.
- Miller-IMH: mean KS 0.0137, mean W1 0.0002, coverage 1.000.
- Miller-Gamma: mean KS 0.0006, mean W1 2.25e-05, coverage 1.000.
- Stirling-Gamma: mean KS 0.9645, mean W1 0.0301, coverage 0.000.
- Slice-log-alpha: mean KS 0.0131, mean W1 0.0002, coverage 1.000.

## Prior Robustness

- P-IG Gibbs skipped 600 prior rows because the current integer-delta augmentation does not cover noninteger posterior delta.

## Truncation Sensitivity

- With N=200, mean truncation-suite KS is 0.0579; with N=1000 it is 0.0585. Median runtime increases from 3.87s to 16.70s.

## Manuscript-Ready Interpretation

The simulation now separates two questions that were conflated in the earlier draft: numerical correctness of the target distribution and computational efficiency of competing samplers. The grid calculation supplies a stable one-dimensional reference for each generated data set, so KS and Wasserstein errors are interpretable as errors relative to the actual posterior rather than discrepancies among Monte Carlo methods.

The exact P-IG Gibbs sampler is robust in the sense that it targets the grid posterior across small, moderate, and large shapes when the integer posterior-delta condition is satisfied. However, the cost grows directly with posterior delta because each iteration samples a sum of delta P-IG variables. The comparison with slice sampling and Miller-based independence MH should therefore be reported in ESS/sec, not only ESS/iteration.

Miller-Gamma is the strongest deterministic approximation in these experiments and Miller-IMH is the most important exact competitor because it uses that approximation only as a proposal while preserving the target. Stirling-Gamma is useful as a historical baseline but should not be presented as competitive in the low-shape regimes if its KS/Wasserstein errors are visibly larger.

The prior robustness experiment exposes the current limitation of the augmentation: noninteger posterior delta is not covered by the present Gibbs construction. The paper should either state the integer-delta condition explicitly or add a separate noninteger construction before claiming full robustness to arbitrary informative priors.
