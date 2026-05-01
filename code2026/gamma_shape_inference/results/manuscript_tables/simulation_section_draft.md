# Gamma Shape Inference Simulation: Manuscript Draft

## Design

We evaluated posterior inference for the gamma shape parameter under the
posterior density proportional to
`Gamma(delta * alpha + 1) / { Gamma(alpha)^delta (delta * eta / mu)^(delta * alpha) }`.
For each data set, a one-dimensional adaptive log-grid quadrature was used as
the numerical reference distribution. Accuracy was measured by the
Kolmogorov-Smirnov distance and the 1-Wasserstein distance to this numerical
reference; efficiency was measured by ESS per iteration, ESS per second,
lag-1/lag-10 autocorrelation, and runtime.

The main simulation crossed `alpha_true in {0.1, 0.5, 1, 2, 5, 20}` with
`n in {10, 50, 500}`, using 200 independent data sets per cell. The competitors
were P-IG Gibbs with default truncation `N=200`, slice sampling on `log(alpha)`,
Miller's gamma approximation used both directly and as an independence-MH
proposal, and the Stirling/Rossell gamma approximation. We also ran a Damsleth
calibration experiment with 100 independent sampler seeds for each reported
sufficient statistic, an informative-prior experiment with 100 data sets for
each prior setting, a P-IG truncation experiment over
`N in {50, 100, 200, 500, 1000}`, and a low-shape stress test at
`alpha_true in {0.05, 0.10}` with `n=100`.

All simulations are reproducible from a fixed seed. The full benchmark runner
uses seed `20260430` by default and initializes each data set or sampler task
with `seed + 1009 * task_id`. The chain-length exactness experiment also uses
seed `20260430`, with deterministic cell-specific data seeds and independent
chain seeds.

## Main Cross-Regime Results

Across the 18 main cells, P-IG Gibbs produced a mean KS distance of 0.0433
and a median KS distance of 0.0250. The two generic exact baselines were more
accurate at the fixed chain length: Miller-IMH had mean KS 0.0137 and slice
sampling had mean KS 0.0143. This difference is Monte Carlo rather than a
change in target distribution: the P-IG chain has much higher autocorrelation
(mean lag-1 autocorrelation 0.677) than Miller-IMH (0.011) or slice sampling
(0.004).

The efficiency comparison is decisive. Median ESS/sec over the main cells was
81.8 for P-IG Gibbs, 20,431 for slice sampling, and 46,917 for Miller-IMH.
The median ratio of P-IG ESS/sec to the fastest exact competitor was 0.0017,
and the ratio was smallest for large shape and large sample size. For example,
at `alpha_true=20, n=500`, P-IG had median ESS/sec 0.85 versus about 47,906
for Miller-IMH. The computational bottleneck is structural: each P-IG iteration
samples a sum of `delta` P-IG auxiliary variables, so runtime increases with
the posterior information in the shape parameter.

Miller's direct gamma approximation was extremely close to the numerical
reference in this one-dimensional problem, with mean KS 0.0024 across the main
cells. This makes Miller-Gamma the strongest deterministic approximation in the
study. The Stirling approximation was substantially less reliable, with mean
KS 0.3358 and mean empirical coverage 0.7375.

These main results should not be interpreted as a failure of the P-IG identity.
They show that, for a one-dimensional posterior where tailor-made proposals are
available, P-IG Gibbs has larger finite-chain Monte Carlo error at the same
nominal iteration count. Its value is that it gives an exact conditional
augmentation for gamma-function terms, not that it is automatically the fastest
1D sampler.

## Exactness and Chain Length

To separate exactness from finite-chain efficiency, we fixed three
representative posterior targets and ran P-IG Gibbs for `5000`, `20000`, and
`100000` iterations with four independent chain seeds. The three targets were
low shape (`alpha_true=0.1, n=100`), moderate shape (`alpha_true=5, n=50`), and
high shape with large sample size (`alpha_true=20, n=500`). Miller-Gamma and
Stirling-Gamma were evaluated once on the same targets as deterministic
reference lines.

The P-IG error decreased systematically as the chain length increased. In the
low-shape target, P-IG KS dropped from 0.0139 at 5000 iterations to 0.0073 at
20000 and 0.0034 at 100000. The Stirling approximation on the same target had
KS 0.954, which cannot be reduced by running longer. In the moderate-shape
target, P-IG KS dropped from 0.0506 to 0.0286 and then 0.0109, while
Stirling-Gamma stayed at 0.0500. In the high-shape large-`n` target, P-IG KS
dropped from 0.1401 to 0.1062 and then 0.0164, whereas Stirling-Gamma was
0.0545.

This experiment is the cleanest way to state the benefit of exact augmentation.
P-IG error is Monte Carlo error around the exact posterior and can be reduced
with more simulation. Approximation error is algorithmic bias: it can be small,
as for Miller-Gamma in these one-dimensional examples, or large, as for
Stirling-Gamma in low-shape regimes, but it does not vanish with more MCMC
iterations.

## Low-Shape Stress Test

The low-shape stress test isolates the regime where Stirling-type expansions
are expected to be fragile. At `alpha_true=0.05` and `0.10` with `n=100`,
P-IG Gibbs, slice sampling, and Miller-IMH all tracked the grid posterior
closely, with mean KS around 0.013-0.014 and empirical coverage 1.00. The
Miller-Gamma approximation was also very accurate in these settings, with mean
KS 0.0006. In contrast, Stirling-Gamma failed badly: KS was 0.978 at
`alpha_true=0.05` and 0.951 at `alpha_true=0.10`, with zero coverage in both
cells. Thus the clear failure mode demonstrated by the simulations is the
Stirling approximation, not Miller's approximation.

## Prior Robustness

For informative priors, the generic methods and deterministic approximations
can be evaluated for both integer and noninteger posterior `delta`. The current
P-IG Gibbs construction, however, introduces `delta` independent P-IG
variables and is exact only for integer posterior `delta`. Accordingly, the
simulation marks all noninteger posterior-delta P-IG cases as skipped:
600 prior-experiment rows were skipped for prior `delta=0.5` and `2.5`.
For integer posterior-delta settings (`delta=10` and `50` in the prior
experiment), P-IG ran normally and remained close to the grid reference.

This result should be stated explicitly in the paper. The present algorithm is
not a general solution for arbitrary noninteger prior information unless an
additional noninteger-delta augmentation is developed.

## Truncation Sensitivity

The truncation experiment did not show a practically meaningful improvement
from increasing the P-IG truncation level beyond the default `N=200` at the
current chain length. Mean KS was 0.0599, 0.0574, 0.0579, 0.0541, and 0.0585
for `N=50, 100, 200, 500, 1000`, respectively. These differences are dominated
by finite-chain Monte Carlo error rather than a monotone visible truncation
bias. In contrast, runtime increased sharply: median runtime was 1.13s,
2.06s, 3.87s, 9.12s, and 16.70s. Thus `N=200` is a reasonable computational
default, but the simulation does not support a claim that larger `N` improves
posterior-sampling accuracy at the reported chain length.

## Damsleth Calibration

The Damsleth examples are useful as a calibration against the earlier paper but
should not be treated as a repeated-data frequency experiment because the
sufficient statistics are fixed. Over 100 independent sampler seeds, P-IG Gibbs
had KS 0.034, 0.049, and 0.034 for `n=5, 10, 30`, respectively. Miller-IMH and
slice sampling had KS around 0.014, and Miller-Gamma had KS around 0.003-0.005.
The Damsleth results therefore confirm that P-IG targets the correct posterior
but also show its higher Monte Carlo cost.

## Recommended Manuscript Claim

The strengthened simulation supports a narrower and more defensible claim than
the earlier draft. P-IG Gibbs is an exact data-augmentation sampler for the
integer-posterior-delta gamma shape problem and remains accurate across
low-shape, moderate-shape, and high-shape regimes when enough Monte Carlo
effort is used. However, it is not the most efficient method for this
one-dimensional inference problem. Miller-IMH is a very strong exact
competitor, and Miller's direct gamma approximation is nearly indistinguishable
from the numerical posterior in the tested regimes. The principal practical
advantage of P-IG is that it converts gamma-function normalizing constants into
an exact auxiliary-variable update that can be embedded in larger Gibbs
samplers. The simulation section should therefore frame P-IG as an exact
augmentation device that avoids approximation bias, while honestly reporting
its computational cost and its current integer-posterior-delta limitation.
