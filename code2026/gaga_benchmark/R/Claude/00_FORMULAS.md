# PIG-augmented Bayesian inference for the GaGa model

This document derives the data-augmentation Gibbs sampler we use to replace
Rossell (2009)'s Stirling-based EM for the GaGa hyperparameters with the
P-IG augmentation of Polson, He, Xu, Zhai (2025+, JASA submission).

## Notation

* Genes $i = 1, \dots, G$.
* Experimental groups $g = 1, \dots, K$. For Armstrong $K = 2$
  (ALL vs MLL); for MAQC $K = 4$ (pools A, C, D, B in titration order).
* Replicates $j = 1, \dots, n_g$ per group.
* Patterns $z_i \in \{0, 1, \dots, P-1\}$ partitioning the $K$ groups
  into clusters. $z=0$ is the all-equal pattern. For Armstrong $P = 2$,
  for MAQC $P = 5$ (the five titration patterns of Rossell §7).
* Under pattern $z_i$ groups $g$ are partitioned into clusters
  $c \in \mathcal{C}(z_i)$. Let $n_{ic}$ be the number of replicates
  in cluster $c$ for gene $i$, $S_{ic} = \sum_{j,g \in c} x_{ijg}$ and
  $L_{ic} = \prod_{j,g \in c} x_{ijg}$ the cluster sufficient statistics.

## Likelihood

We adopt Rossell's parameterization,
$$x_{ijg} \mid \alpha, \lambda_{ic} \;\sim\; \mathrm{Ga}\bigl(\alpha, \alpha/\lambda_{ic}\bigr),$$
so the cluster mean is $\mathbb{E}[x_{ijg}] = \lambda_{ic}$ and the
coefficient of variation is $1/\sqrt{\alpha}$. The shape $\alpha$ is shared
across genes (the `equalcv = TRUE` setting Rossell uses for the Armstrong
analysis); we use a single global $\alpha$ throughout.

## Hierarchical prior

On the per-cluster mean we use Rossell's prior
$$\frac{1}{\lambda_{ic}} \;\sim\; \mathrm{Ga}\!\left(\alpha_0,\; \alpha_0\lambda_0\right),
\qquad i = 1,\dots,G,\; c \in \mathcal{C}(z_i),$$
with hyperparameters $(\alpha_0, \lambda_0)$.

On the global shape $\alpha$ we use the Damsleth $\xi_2(\eta/\mu, \delta)$
type-II conjugate prior (PIG paper §4.1, eq. (4.1)),
$$p(\alpha \mid \mu, \delta, \eta) \;\propto\;
  \frac{\Gamma(\delta\alpha+1)}{\Gamma(\alpha)^\delta}\,
  (\delta\mu)^{-\delta\alpha}\, e^{-\delta\eta\alpha},$$
where the hyperparameters $(\delta, \eta, \mu)$ are estimated by method
of moments from the full data (a fully Bayesian extension would put a
hyperprior on $\delta$).

On $\alpha_0$ we use a similar Damsleth prior; on $\lambda_0$ we use a
diffuse inverse-gamma; on the pattern-probability vector $\boldsymbol\pi$
a symmetric Dirichlet.

## Conditional posteriors

### 1. Cluster means $\lambda_{ic}$

Conditional on $(\alpha, \alpha_0, \lambda_0, z_i)$ and the data,
$$\frac{1}{\lambda_{ic}} \;\Big\vert\; \cdot \;\sim\;
  \mathrm{Ga}\!\bigl(\alpha_0 + \alpha n_{ic},\; \alpha_0\lambda_0 + \alpha S_{ic}\bigr).$$

### 2. Pattern indicators $z_i$

For each gene we compute the pattern-conditional marginal likelihood
$$m_i(z) = \prod_{c \in \mathcal{C}(z)}
  \frac{\alpha^{\alpha n_{ic}} L_{ic}^{\alpha-1}}{\Gamma(\alpha)^{n_{ic}}}\,
  \frac{(\alpha_0\lambda_0)^{\alpha_0}}{\Gamma(\alpha_0)}\,
  \frac{\Gamma(\alpha_0 + \alpha n_{ic})}{(\alpha_0\lambda_0 + \alpha S_{ic})^{\alpha_0+\alpha n_{ic}}}.$$
Then $\Pr(z_i = z \mid \cdot) \propto \pi_z\, m_i(z)$ and we sample
$z_i$ from this categorical distribution (or, in EM mode, take the mode).

### 3. Pattern probabilities $\boldsymbol\pi$

With Dirichlet$(\boldsymbol\nu)$ prior and pattern counts $n_z$,
$$\boldsymbol\pi \mid \mathbf{z} \sim \mathrm{Dirichlet}(\boldsymbol\nu + \mathbf{n}).$$

### 4. Hyperprior mean $\lambda_0$

With prior $1/\lambda_0 \sim \mathrm{Ga}(a_0, b_0)$ and after summing
$\sum_{i,c} 1/\lambda_{ic}$,
$$\frac{1}{\lambda_0} \;\Big\vert\; \cdot \;\sim\;
  \mathrm{Ga}\!\Bigl(a_0 + \alpha_0 \sum_i |\mathcal{C}(z_i)|,\;
  b_0 + \alpha_0 \sum_{i,c} 1/\lambda_{ic}\Bigr).$$

### 5. Global shape $\alpha$ — the P-IG augmentation

Conditional on $\boldsymbol\lambda$ and the data, the gamma sampling
likelihood gives, after dropping factors that are constant in $\alpha$,
$$p(\alpha \mid \boldsymbol\lambda, \mathbf{x}) \;\propto\;
  \frac{\alpha^{\alpha N_{\mathrm{tot}}}}{\Gamma(\alpha)^{N_{\mathrm{tot}}}}\,
  e^{\alpha D_\alpha}\, p(\alpha),$$
where $N_{\mathrm{tot}} = \sum_{i,j,g} 1$ is the total number of
expression measurements and
$$D_\alpha = \sum_{i,j,g} \log x_{ijg} - \sum_{i,c} n_{ic} \log\lambda_{ic}
            - \sum_{i,c} S_{ic}/\lambda_{ic}.$$

Combining with the Damsleth prior on $\alpha$ (which contributes
$\Gamma(\delta\alpha+1)/\Gamma(\alpha)^\delta \cdot (\delta\mu)^{-\delta\alpha}\,
 e^{-\delta\eta\alpha}$), the unnormalised posterior is
$$p(\alpha \mid \cdot) \;\propto\;
  \frac{\alpha^{\alpha N_{\mathrm{tot}}}\,\Gamma(\delta\alpha+1)}
       {\Gamma(\alpha)^{N_{\mathrm{tot}}+\delta}}\,
  e^{\alpha\,(D_\alpha - \delta\eta - \delta\log(\delta\mu))}.$$

The intrinsic difficulty is the $\alpha^{\alpha N_{\mathrm{tot}}}$
factor — this is the gamma-density normalising constant raised to a
power and is the same obstacle Rossell handled via Stirling. We deal
with it via the **change-of-variable trick**: write $\phi_{ic} = \alpha/\lambda_{ic}$
and **sample $\phi_{ic}$ instead of $\lambda_{ic}$** in the Gibbs chain.
With $\phi_{ic}$ as the natural rate parameter of the gamma sampling
model,
$$x_{ijg} \mid \alpha, \phi_{ic} \;\sim\; \mathrm{Ga}(\alpha, \phi_{ic}),$$
the data contribute $\phi_{ic}^{\alpha n_{ic}}/\Gamma(\alpha)^{n_{ic}}$
per cluster, with **no $\alpha^{\alpha n}$**. The price is that the
prior on $\phi_{ic}$, induced by Rossell's prior on $\lambda_{ic}$, is
$\alpha$-dependent. Specifically, $\phi_{ic} = \alpha/\lambda_{ic}$ with
$1/\lambda_{ic} \sim \mathrm{Ga}(\alpha_0, \alpha_0\lambda_0)$ implies
$$\phi_{ic} \mid \alpha, \alpha_0, \lambda_0 \;\sim\;
  \mathrm{Ga}\!\bigl(\alpha_0,\; \alpha_0\lambda_0/\alpha\bigr),$$
so the prior contributes $\phi^{\alpha_0-1} e^{-\phi \alpha_0\lambda_0/\alpha}
\cdot (\alpha_0\lambda_0/\alpha)^{\alpha_0}/\Gamma(\alpha_0)$ per
cluster — the $\alpha^{-\alpha_0}$ from the normaliser is benign.

After this reparameterisation, the conditional density of $\alpha$ given
$\boldsymbol\phi$ and data has the form
$$p(\alpha \mid \boldsymbol\phi, \mathbf{x})
  \;\propto\; \frac{\Gamma(\delta\alpha + 1)}{\Gamma(\alpha)^{N_{\mathrm{tot}}+\delta}}\,
  e^{\alpha\,\widetilde{D}}\, ,$$
where
$$\widetilde{D} = \sum_{ijg}\log x_{ijg}
             + \sum_{i,c} n_{ic}\log\phi_{ic}
             - \sum_{i,c} S_{ic}\phi_{ic}
             - \delta\eta - \delta\log(\delta\mu)
             - \alpha_0 \sum_{i,c} \log(\alpha_0\lambda_0/\alpha) \cdot 0
             $$
(the rate-prior log-normaliser term cancels exactly with the $\phi$
density contribution because both scale linearly with $\alpha_0\log\alpha$;
see the appendix). This is precisely the GaS form of PIG paper §4.2:
$$p(\alpha) \;\propto\; \frac{\Gamma(\alpha + \beta)^S}{\Gamma(\alpha)^S}\,e^{c\alpha}$$
with $\beta = 0$ (so $\Gamma(\alpha)^{-S}$ alone after applying the
Damsleth $\Gamma(\delta\alpha+1)$ via $S = N_{\mathrm{tot}}+\delta$
and one auxiliary $\tau \sim \mathrm{Ga}(\delta\alpha+1, 1)$).

The Gibbs update for $\alpha$ is therefore
$$\begin{aligned}
\tau \mid \alpha &\sim \mathrm{Ga}(\delta\alpha + 1,\; 1),\\
\omega_s \mid \alpha &\stackrel{\mathrm{iid}}{\sim} \mathrm{P\text{-}IG}(\alpha),
  \quad s = 1,\dots, N_{\mathrm{tot}}+\delta,\\
\alpha \mid \tau, \boldsymbol\omega &\sim
  \mathrm{PTN}\!\bigl(\delta + 1,\; \tilde{a},\; \tilde{b}\bigr),
\end{aligned}$$
with
$\tilde{a} = \sum_s \omega_s$ and
$\tilde{b} = (N_{\mathrm{tot}} + \delta)\,\gamma + \delta\log\tau + \widetilde{D}$.
$\gamma = -\psi(1)$ is the Euler–Mascheroni constant; PTN is the
power-truncated normal $f(x) \propto x^{p-1} e^{-ax^2 + bx}$ on $x>0$,
sampled with `rH(p, a, b)`.

### 6. Hyperprior shape $\alpha_0$

By the same logic applied to the prior layer
$1/\lambda_{ic} \sim \mathrm{Ga}(\alpha_0, \alpha_0\lambda_0)$, after
the substitution $\rho_{ic} = \alpha_0\lambda_0/\lambda_{ic}$ (rate of
the inverse-$\lambda$ gamma) the conditional of $\alpha_0$ becomes
$$p(\alpha_0 \mid \boldsymbol\rho) \;\propto\;
  \frac{1}{\Gamma(\alpha_0)^{M}}\,
  e^{\alpha_0\,\widetilde{E}},$$
with $M = \sum_i |\mathcal{C}(z_i)|$ the total number of cluster-gene
slots and
$$\widetilde{E} = \sum_{i,c} \log\rho_{ic}
             - \sum_{i,c} \rho_{ic}/(\lambda_0)
             - \log(\text{prior hyperparameter}).$$
The P-IG augmented update is the same shape as for $\alpha$:
$\omega_s^\prime \sim \mathrm{P\text{-}IG}(\alpha_0)$ for $s=1,\dots,M$
and $\alpha_0 \sim \mathrm{PTN}(p, \sum \omega_s^\prime, M\gamma + \widetilde{E})$.

## Algorithm summary

```
Inputs:  matrix X (G × n) of expression, group assignment g(j), patterns
         (P × K) integer matrix.
Initialise α, α_0, λ_0, π, λ_{ic}, z_i.
Repeat for t = 1, …, T iterations:
  1. For each gene i: compute m_i(z) for each pattern z, sample z_i from
     the categorical posterior.
  2. Sample π | z from Dirichlet.
  3. For each gene i, cluster c: sample 1/λ_{ic} from Ga.
  4. Sample 1/λ_0 from Ga.
  5. Update α via PIG augmentation (steps 5a-5c above), using the
     change-of-variable trick on φ_{ic} = α/λ_{ic}.
  6. Update α_0 via PIG augmentation analogously.
After burn-in, average pattern-membership probabilities over Gibbs
samples to obtain posterior P(z_i = 0 | data) for each gene; use these
for FDR-controlled DE calls.
```

## Comparison to Rossell's Stirling EM

* Rossell uses Stirling's approximation
  $\log \Gamma(\alpha) \approx (\alpha - 1/2)\log\alpha - \alpha + \tfrac{1}{2}\log(2\pi)$
  to fit a gamma distribution to the GaS density and then runs an
  EM-type algorithm. This is fast but introduces approximation error
  that grows in regions where the GaS density is far from gamma
  (small $\alpha$, equivalently high CV).
* PIG augmentation samples $\alpha$ exactly from its conditional
  posterior, with no approximation. At the cost of one P-IG draw per
  observation per Gibbs iteration, we get exact MCMC samples that
  automatically integrate over uncertainty in $\alpha$ when computing
  per-gene pattern-membership probabilities.
* The downstream DE-calling pipeline is identical: rank genes by
  $1 - \Pr(z_i = 0 \mid x)$ and threshold at the level that controls
  posterior FDR.
