# Quad-GaGa Formulas and PIG Sampler Status

This folder is isolated from the existing `code2026` scripts. It implements a
shape-integration benchmark for the two real-data examples used in the GaGa
paper: Armstrong leukemia and MAQC. The empirical scripts in this folder use
adaptive numerical quadrature, so the method label is `Quad_GaGa`, not
`PIG_GaGa`.

## GaGa Model

For gene `i`, sample `j`, and biological group `z_j`, GaGa uses

```text
x_ij | alpha_i, lambda_i,b ~ Gamma(alpha_i, rate = alpha_i / lambda_i,b)
lambda_i,b | a0, nu       ~ InvGamma(a0, scale = a0 / nu)
alpha_i | beta, mu        ~ Gamma(beta, rate = beta / mu)
delta_i | pi              ~ Categorical(pi).
```

Here `b` is the block induced by an expression pattern. With `equalcv=TRUE`,
`alpha_i` is shared across the blocks of a gene.

For a candidate pattern `h`, let `B_h` be its blocks. For block `b`,

```text
n_b = number of arrays in block b
S_b = sum_j x_ij over block b
L_b = sum_j log(x_ij) over block b
r   = a0 / nu
```

After integrating out `lambda_i,b`, the collapsed alpha log-kernel is

```text
ell_ih(alpha) =
  (beta - 1) log(alpha) - (beta / mu) alpha
  + sum_b [
      n_b alpha log(alpha)
      - n_b log Gamma(alpha)
      + (alpha - 1) L_b
      + log Gamma(a0 + n_b alpha)
      - (a0 + n_b alpha) log(r + alpha S_b)
    ]
```

Terms not depending on `alpha` are retained in the code when they differ across
patterns:

```text
sum_b [ a0 log(r) - log Gamma(a0) ] + beta log(beta / mu) - log Gamma(beta).
```

The quadrature-based pattern posterior is then

```text
Pr(delta_i = h | x_i, theta)
  proportional to pi_h * integral_0^infinity exp(ell_ih(alpha)) d alpha.
```

In the empirical scripts, `theta = (a0, nu, beta, mu, pi)` is initialized from
the original `gaga::fitGG(..., nclust=1)` estimate. This makes the benchmark a
shape-integration replacement benchmark, not yet a full replacement for GaGa's
C-level EM fitting routine.

## PIG Shape Sampler

The PIG paper gives the reciprocal-gamma identity

```text
1 / Gamma(alpha) =
  alpha exp(gamma alpha) * E_{omega ~ PIG(0)}[ exp(-omega alpha^2) ],
```

and the tilted conditional

```text
omega | alpha ~ PIG(alpha).
```

For the GaGa gamma-shape distribution discussed in the paper,

```text
p(alpha) proportional to prod_s Gamma(alpha + beta_s) / Gamma(alpha)^S
                       * exp(c alpha),
```

the PIG Gibbs updates are

```text
tau_s | alpha ~ Gamma(alpha + beta_s, 1)
omega_s | alpha ~ PIG(alpha)
alpha | tau, omega ~ PTN(S + 1, sum_s omega_s,
                         S * gamma_const + c + sum_s log(tau_s)).
```

The library file implements this sampler as `pig_gas_sampler()`, but the
Armstrong and MAQC benchmark scripts do not call it. They use
`quad_gaga_pp()`, which evaluates the same collapsed alpha integral with
adaptive log-scale quadrature. Therefore these outputs should be interpreted as
GaGa with numerical shape integration, not as a PIG-augmented sampler.

## Empirical Outputs

The scripts write outputs under:

```text
code2026/gaga_benchmark/R/ChatGPT/results/
```

Key files:

```text
armstrong_quad_gaga_full_summary.csv
armstrong_quad_gaga_scores.csv
armstrong_quad_gaga_reproducibility_summary.csv
maqc_quad_gaga_validation_summary.csv
maqc_quad_gaga_validation_curve.csv
maqc_quad_gaga_validation_curve.png
maqc_quad_gaga_probe_scores.csv
```
