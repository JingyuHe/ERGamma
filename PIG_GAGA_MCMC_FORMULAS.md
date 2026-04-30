# PIG-augmented GaGa model: formulas for the implementation

This document specifies the version I will implement for the Armstrong and
MAQC benchmarks.  It is not the quadrature/Laplace version.  Posterior pattern
probabilities are estimated from MCMC draws, and every alpha update uses the
PIG augmentation through `pig_rPIG_original()` and a direct log-alpha slice
update of the augmented conditional.

## 1. GaGa model matched to Rossell/gaga

For gene `i`, sample `j`, and experimental group `g`, pattern `z_i = h`
partitions the groups into clusters `c in C_h`.  Let `c_h(g)` be the cluster
containing group `g`.

```text
x_ijg | alpha_i, lambda_i,c_h(g)
  ~ Gamma(shape = alpha_i, rate = alpha_i / lambda_i,c_h(g))

q_ic = 1 / lambda_ic | a0, nu
  ~ Gamma(shape = a0, rate = a0 / nu)

alpha_i | b_alpha, nu_alpha
  ~ Gamma(shape = b_alpha, rate = b_alpha / nu_alpha)

z_i | pi ~ Categorical(pi_0, ..., pi_{H-1})
```

This uses the same `nu` parameterization as `gaga::simGG`:

```text
lambda_ic = 1 / Gamma(a0, rate = a0 / nu).
```

For a fixed pattern `h` and cluster `c`, define sufficient statistics

```text
n_ic = number of samples in cluster c,
S_ic = sum_{j,g in c} x_ijg,
L_ic = sum_{j,g in c} log x_ijg,
r0   = a0 / nu,
N_i  = sum_c n_ic.
```

## 2. Collapsed cluster likelihood

Integrating out `q_ic = 1/lambda_ic`,

```text
m_ic(alpha_i)
 = int p(x_ic | alpha_i, q_ic) p(q_ic | a0, nu) dq_ic

 = alpha_i^(n_ic alpha_i)
   exp((alpha_i - 1) L_ic)
   / Gamma(alpha_i)^n_ic
   * r0^a0 / Gamma(a0)
   * Gamma(a0 + n_ic alpha_i)
   / (r0 + alpha_i S_ic)^(a0 + n_ic alpha_i).
```

The pattern-conditional likelihood is

```text
m_ih(alpha_i) = prod_{c in C_h} m_ic(alpha_i).
```

For fixed `alpha_i`, the pattern update is exact:

```text
Pr(z_i = h | alpha_i, x_i, theta)
  = pi_h m_ih(alpha_i) / sum_l pi_l m_il(alpha_i).
```

This step uses no numerical integration.

## 3. Conditional target for alpha_i

Given `z_i = h`, the collapsed target for `alpha_i` is proportional to

```text
p(alpha_i | z_i = h, x_i, theta)
  proportional to
  alpha_i^(b_alpha - 1) exp(-(b_alpha / nu_alpha) alpha_i)
  * prod_c m_ic(alpha_i).
```

Dropping constants independent of `alpha = alpha_i`, write

```text
log p(alpha | z_i, x_i, theta)
 = (b_alpha - 1) log alpha - (b_alpha / nu_alpha) alpha
   + sum_c [
       n_c alpha log alpha
       - n_c log Gamma(alpha)
       + (alpha - 1) L_c
       + log Gamma(a0 + n_c alpha)
       - (a0 + n_c alpha) log(r0 + alpha S_c)
     ].
```

The hard factor is the gamma-shape part:

```text
prod_c Gamma(a0 + n_c alpha) / Gamma(alpha)^n_c.
```

## 4. PIG augmentation for alpha_i

For the denominator, use the PIG identity

```text
1 / Gamma(alpha)^N
  proportional to
  alpha^N exp(N gamma_const alpha)
  E_{omega_s}[ exp(-alpha^2 sum_s omega_s) ],

omega_s | alpha ~ PIG(alpha),  s = 1, ..., N,
gamma_const = -digamma(1).
```

For each numerator factor, use the Damsleth gamma augmentation

```text
Gamma(a0 + n_c alpha)
  = int tau_c^(a0 + n_c alpha - 1) exp(-tau_c) d tau_c,

tau_c | alpha ~ Gamma(shape = a0 + n_c alpha, rate = 1).
```

Given current `alpha`, draw

```text
omega_1, ..., omega_N | alpha ~ iid PIG(alpha),
tau_c | alpha ~ Gamma(a0 + n_c alpha, 1),  for c in C_h.
```

With `Omega = sum_s omega_s`, the augmented conditional can be factored as

```text
p(alpha | tau, omega, z_i, x_i, theta)
  proportional to
  alpha^(N + b_alpha - 1)
  exp(-Omega alpha^2 + B alpha)
  * exp(R(alpha)),
```

where the PTN part is

```text
B = N gamma_const
    + sum_c L_c
    - b_alpha / nu_alpha
    + sum_c n_c log tau_c,
```

and the remaining non-PTN factor is

```text
R(alpha)
  = N alpha log alpha
    - sum_c (a0 + n_c alpha) log(r0 + alpha S_c).
```

## 5. Exact log-alpha slice update used in the benchmark

The runnable benchmark samples this one-dimensional conditional on
`y = log(alpha)`.  Including the Jacobian, the augmented log density is

```text
ell(y)
  = (N + b_alpha) y
    - Omega exp(2y)
    + B exp(y)
    + R(exp(y)).
```

The code applies stepping-out and shrinkage slice sampling to `ell(y)`.  This
is an exact MCMC update for the PIG-augmented conditional and avoids the
zero-acceptance behavior of the plain PTN independence proposal on the real
microarray data.

## 6. Optional non-tangent PTN proposal

The implementation uses the same non-tangent proposal as
`GAGA_ORIGINAL_PIG_MCMC_FORMULAS.tex`:

```text
alpha_star | tau, omega
  ~ PTN(p = N + b_alpha,
        a = Omega,
        b = B),

PTN(p, a, b): density proportional to
  x^(p - 1) exp(-a x^2 + b x),  x > 0.
```

Conditional on the sampled auxiliaries, this proposal does not depend on the
current `alpha`.  The Metropolis-Hastings acceptance probability is therefore

```text
A(alpha -> alpha_star)
  = min(1, exp(R(alpha_star) - R(alpha_old))).
```

A tangent PTN proposal can be useful for efficiency, but it is state-dependent.
It would need the full reverse/forward proposal density ratio, including PTN
normalizing constants.  The diagnostic `PIG_MCMC_ALPHA_KERNEL=ptn` uses the
exact non-tangent proposal above, but it is not the benchmark default because
it barely moves on the real data.

## 7. Full per-gene MCMC for posterior pattern probabilities

For each gene `i`, with fixed hyperparameters
`theta = (a0, nu, b_alpha, nu_alpha, pi)`:

```text
Initialize alpha_i > 0 and z_i.

For t = 1, ..., T:
  1. Sample z_i from Pr(z_i = h | alpha_i, x_i, theta).
  2. Given z_i, draw tau_c and omega_s.
  3. Sample y_i = log alpha_i by slice sampling from ell(y_i).
  4. Set alpha_i = exp(y_i).
  5. After burn-in, record z_i and the conditional probability vector
       Pr(z_i = h | alpha_i, x_i, theta).

Posterior pattern probabilities:
  pp_ih = mean_t Pr(z_i = h | alpha_i^(t), x_i, theta).

DE score:
  score_i = 1 - pp_i0.
```

The last line is a Rao-Blackwellized MCMC estimate.  The code also stores hard
`z_i` counts for diagnostics, but the benchmark uses the smoother conditional
probability average because the null pattern can be rare in these data.

This is the first benchmark target, because it holds the GaGa hyperparameters
fixed and replaces only the shape/posterior computation with PIG-augmented
MCMC.  In code, the default hyperparameters will come from `gaga::fitGG` so the
comparison isolates the PIG posterior step from hyperparameter estimation.

## 8. Optional full Bayesian/MCEM extension

If we want to fit all GaGa hyperparameters rather than borrow them from
`gaga::fitGG`, keep `q_ic = 1/lambda_ic` in the chain:

```text
q_ic | alpha_i, z_i, x_i, a0, nu
  ~ Gamma(a0 + n_ic alpha_i, rate = a0 / nu + alpha_i S_ic).
```

Then

```text
pi | z ~ Dirichlet(epsilon + n_0, ..., epsilon + n_{H-1}).
```

For `rho = 1/nu`, using `rho ~ Gamma(a_rho, b_rho)`,

```text
rho | q, a0
  ~ Gamma(a_rho + M a0, b_rho + a0 sum_{i,c} q_ic),
nu = 1 / rho.
```

The shape hyperparameters `a0` and `b_alpha` have gamma-shape conditionals with
`1/Gamma(a0)^M` and `1/Gamma(b_alpha)^G` terms.  They can be updated by the same
PIG/PTN-plus-remainder MH idea, or by a conservative log-scale MH step.  For the
Rossell benchmark I will start with the fixed-hyperparameter version above,
then add hyperparameter sampling only after that is validated against
`gaga::fitGG` on simulated data.
