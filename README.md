# Bayesian Stochastic State-Space Epidemiological Models

A Bayesian state-space compartmental modeling framework implemented in **Stan** and **R**, evaluating how biological assumptions, viral genomic signals, vaccination effects, and non-pharmaceutical interventions influence inferred SARS-CoV-2 transmission dynamics across global regions.

**Authors:** Lucas Anderson, Ian McArthur
**Course:** Bio 415 Project Analysis Pipeline

---

## Contributions

- **Luke Anderson:** Model design and biological/statistical formulation, Stan/R implementation, simulation and Bayesian inference pipeline, MCMC diagnostics, data processing, results visualization, drafting of the Results section, and portions of the Methods / Model Ladder section.
- **Ian McArthur:** Manuscript writing for the Abstract, Introduction, Discussion, and Model description sections, and contributions to the project presentation.

Full methodological and results write-up is provided in [`paper/`](./paper).

---

## Project Overview

SIR and SEIR models are widely used to study infectious disease dynamics but often rely on simplifying assumptions such as constant transmission rates and exponentially distributed transition times. This project extends stochastic compartmental models by incorporating:

- Viral genomic signals as time-dependent transmission drivers
- Genomic memory effects through latent forcing states
- Vaccination-driven changes in susceptibility
- Government stringency effects as external intervention forcing
- Erlang-distributed latent and infectious periods
- Partial immunity and reinfection pathways
- Bayesian parameter estimation with uncertainty quantification

The framework evaluates how increasing biological realism affects **model convergence**, **parameter inference**, **transmission dynamics**, and **predictive performance**.

---

## Scientific Question

How do biological assumptions and external forcing mechanisms affect inference of epidemic transmission dynamics?

Specifically, this framework evaluates:

- Whether viral genomic variation improves transmission inference
- Whether genomic memory states better represent evolutionary pressure than fixed lags
- Whether Erlang-distributed dwell times improve epidemic state representation
- How vaccination and policy interventions alter inferred genomic coupling

---

## 24-Model Factorial Design

The model suite follows a **4 (biological mechanism) x 2 (compartment structure) x 3 (intervention tier)** factorial design, yielding 24 fitted models.

| Tier | Included forcing |
|---|---|
| **GENO** | Viral genomic signal only |
| **GENO+VAX** | Viral genomic signal + vaccination effects |
| **FULL** | Viral genomic signal + vaccination + government stringency |

This design isolates the individual contribution of biological and policy mechanisms through structured model ablation.

---

## Biological Model Hierarchy

Each biological complexity tier pairs an SIRS (odd ID) and SEIRS (even ID) variant sharing the same forcing/dwell-time structure:

| Model | Structure | Description |
|---|---|---|
| M1 | SIRS (lagged genomic forcing) | Lagged genomic signal, no exposed compartment |
| M2 | SEIRS (lagged genomic forcing) | Exposed compartment + lagged genomic signal |
| M3 | SIRS + genomic memory | Accumulated, decaying genomic forcing state |
| M4 | SEIRS + genomic memory | Accumulated, decaying genomic forcing state |
| M5 | SIRS + Erlang dwell times | Multi-stage infectious period (k = 4) |
| M6 | SEIRS + Erlang dwell times | Multi-stage exposed and infectious periods (k = 4) |
| M7 | SIRS + Erlang + partial immunity | Erlang dwell times + reinfection via infectious pathway |
| M8 | SEIRS + Erlang + partial immunity | Erlang dwell times + reinfection via exposed pathway |

---

## 24-Model Matrix

| Biological Complexity | Structure | GENO | GENO+VAX | FULL |
|---|---|---|---|---|
| M1/M2 (Lagged) | SEIRS | `GENO_M2_SEIRS.stan` | `NAKED_M2_SEIRS_NoStringency.stan` | `SEIR_nu_stringency.stan` |
| M1/M2 (Lagged) | SIRS | `SIR_nu_stringency_GENO.stan` | `SIR_nu_stringency_GENO_VAX.stan` | `SIR_nu_stringency.stan` |
| M3/M4 (Memory) | SEIRS | `GENO_M4_SEIRS_Memory.stan` | `NAKED_M4_SEIRS_Memory_NoStringency.stan` | `SEIR_nu_smooth.stan` |
| M3/M4 (Memory) | SIRS | `SIR_nu_smooth_GENO.stan` | `SIR_nu_smooth_GENO_VAX.stan` | `SIR_nu_smooth.stan` |
| M5/M6 (Erlang) | SEIRS | `GENO_M6_Erlang_SEIRS.stan` | `NAKED_M6_Erlang_SEIRS_NoStringency.stan` | `SEIR_Erlang.stan` |
| M5/M6 (Erlang) | SIRS | `SIR_Erlang_GENO.stan` | `SIR_Erlang_GENO_VAX.stan` | `SIR_Erlang.stan` |
| M7/M8 (Erlang + Waning) | SEIRS | `GENO_M8_Erlang_SEIRS_Immfrac.stan` | `NAKED_M8_Erlang_SEIRS_Immfrac_NoStringency.stan` | `SEIR_Erlang_waning.stan` |
| M7/M8 (Erlang + Waning) | SIRS | `SIR_Erlang_waning_GENO.stan` | `SIR_Erlang_waning_GENO_VAX.stan` | `SIR_Erlang_waning.stan` |

---

## Model Assumptions

### Base SIRS / SEIRS structure

The population is closed (no births, deaths, or migration).

```text
SIRS:  Susceptible -> Infected -> Recovered -> Susceptible
SEIRS: Susceptible -> Exposed -> Infected -> Recovered -> Susceptible
```

New infections arise through contact between susceptible and infectious individuals. The transmission rate is modulated by external forcing signals; recovered individuals lose immunity and return to susceptibility. Vaccination moves susceptible individuals directly into the recovered compartment.

**Vaccination deduction.** Each week, a bounded flow of susceptibles is diverted to recovered before the epidemic transition step is applied:

$$
\Delta_{vax} = \min(S_t, \nu_{vax,t} \cdot N), \qquad S_t^{*} = S_t - \Delta_{vax}
$$

where $\nu_{vax,t}$ is the weekly incident vaccination proportion and $S_t^{*}$ is the post-vaccination susceptible pool used in the infection step.

**Discrete transition probabilities.** All continuous-time rates are converted to weekly transition probabilities using the exact survival probability $1-e^{-r}$, keeping every flow in $[0,1]$:

$$
p_\lambda = 1-e^{-\lambda_t}, \quad p_\sigma = 1-e^{-\sigma}, \quad p_\gamma = 1-e^{-\gamma}, \quad p_\omega = 1-e^{-\omega}
$$

with per-stage Erlang analogues $p_{\sigma_k}=1-e^{-k\sigma}$ and $p_{\gamma_k}=1-e^{-k\gamma}$ used in M5-M8 (see below).

### Genomic memory (M3, M4)

Genomic pressure accumulates and decays over time rather than depending only on the prior week's signal:

$$
A_t = (1-\delta) A_{t-1} + \nu_{signal,t}, \qquad A_1 = \nu_{signal,1}
$$

where $A_t$ is the accumulated genomic forcing state, $\delta$ is the memory decay rate, and $\nu_{signal,t}$ is the standardized weekly genomic volatility signal. Tier 1 (M1, M2) instead uses a fixed one-week-lagged signal $\nu_{lag,t}$ directly in place of $A_t$.

### Erlang dwell times (M5-M8)

Latent and infectious periods are modeled as Erlang-distributed with $k = 4$ sequential sub-stages instead of a single exponential compartment, each with per-stage rate $k\gamma$ (infectious) or $k\sigma$ (exposed):

```text
I1 -> I2 -> I3 -> I4 -> R
```

M5/M6 apply this to SIRS infectious periods and SEIRS exposed+infectious periods, respectively. M7/M8 extend M5/M6 with partial immunity.

### Partial immunity (M7, M8)

A fraction $\eta = 0.5$ of recovered individuals retain partial immunity and re-enter the infection process directly rather than passing through full susceptibility, while the remaining $(1-\eta)$ fraction wanes normally back to S:

$$
R \rightarrow I_1 \quad \text{(M7, SIRS pathway)}
$$

$$
R \rightarrow E_1 \quad \text{(M8, SEIRS pathway)}
$$

---

## Bayesian Inference Framework

All models are implemented in **Stan** and fitted with **CmdStanR** using Hamiltonian Monte Carlo via the No-U-Turn Sampler (NUTS). The pipeline performs Bayesian posterior estimation, parameter uncertainty quantification, posterior predictive evaluation, structural model comparison, and convergence assessment.

### Convergence criteria

- Rank-normalized split-$\hat R < 1.01$ (models reported in the paper required $\hat R < 1.1$ for acceptance)
- Bulk and tail effective sample size (target > 400 across chains)
- Divergent transition count
- Maximum treedepth saturation
- Energy / E-BFMI diagnostics
- Posterior predictive checks against observed weekly case trajectories

---

## Transmission Model

$$
\beta_t = \beta_0 \exp\big(\min[\alpha_\nu \Theta_t + \alpha_{str} Z_t,\ 10]\big)
$$

$$
\lambda_t = \beta_t \cdot \frac{I_t}{N}
$$

- $\beta_0$: baseline transmission rate, prior $\beta_0 \sim \text{LogNormal}(0, 0.5)$
- $\Theta_t$: genomic forcing signal — $\nu_{lag,t}$ for M1/M2, memory state $A_t$ for M3-M8
- $Z_t$: standardized government stringency index (FULL tier only)
- $\alpha_\nu$: genomic coupling coefficient, prior $\mathcal N(0, 0.3)$
- $\alpha_{str}$: intervention effect coefficient, prior $\mathcal N(0, 0.3)$, FULL tier only
- $\lambda_t$: resulting force of infection applied to the susceptible pool; $I_t=\sum_j I_{t,j}$ for Erlang models

The exponent is clipped at 10 to prevent numerical overflow during sampling.

---

## Observation Model

$$
y_t \sim \text{NegBin2}(\mu_t, \phi^{-1})
$$

The expected observed cases $\mu_t$ are new-infection flow scaled by the region-specific reporting rate $\rho$, with the flow term depending on model structure:

- SIRS backbone (M1, M3, M5): $\mu_t = \Delta_{SI} \cdot \rho$
- M7 (partial immunity): $\mu_t = (\Delta_{SI} + \Delta_{RI}) \cdot \rho$
- SEIRS exponential (M2, M4): $\mu_t = \Delta_{EI} \cdot \rho$
- SEIRS Erlang (M6, M8): $\mu_t = \Delta_{E_4} \cdot \rho$

- $y_t$: observed weekly cases
- $\mu_t$: expected observed cases (defined above)
- $\phi^{-1}$: inverse overdispersion parameter, prior $\phi^{-1}\sim\text{Gamma}(2,0.1)$, accounting for the overdispersed clumping typical of COVID-19 case data

---

## Initial Conditions & Fixed Parameters

| Symbol | Value / Prior | Definition |
|---|---|---|
| $s_0$ | $\text{Beta}(3,3)$ | Initial susceptible fraction; $S_1 = N\cdot s_0$ |
| $i_0$ | $\text{Beta}(1,80)$ | Initial infectious fraction; $I_1 = N\cdot i_0$ |
| $e_0$ | $\text{Beta}(1,80)$, SEIR only | Initial exposed fraction; $E_1 = N\cdot e_0$ |
| $\delta$ | $\text{Beta}(2,5)$, M3-M8 | Genomic memory decay rate |
| $\sigma$ | $1/5\ \text{wk}^{-1}$ (fixed) | Incubation rate; mean latent period 5 days |
| $\gamma$ | $1/7\ \text{wk}^{-1}$ (fixed) | Recovery rate; mean infectious period 7 days |
| $\omega$ | $1/180\ \text{wk}^{-1}$ (fixed) | Waning immunity rate; mean immune duration 180 days |
| $k$ | 4 (fixed) | Erlang shape parameter, M5-M8 |
| $\eta$ | 0.5 (fixed) | Partial immunity fraction, M7/M8 |
| $N$ | fixed | Total regional population |
| $\rho$ | fixed | Region-specific case reporting rate |

---

## Sensitivity (Ablation) Analysis

The three intervention tiers are themselves the ablation design: each successive tier reintroduces one covariate to isolate its effect on the genomic coupling coefficient $\alpha_\nu$.

- **GENO** tier (files prefixed `GENO_...` / `SIR_..._GENO.stan`): both $\alpha_{str}$ and vaccination are removed, using the raw susceptible pool $S_t$ in place of $S_t^{*}$:

$$
\lambda_t^{\text{GENO}} = \beta_0 \exp(\alpha_\nu \Theta_t) \cdot \frac{I_t}{N}
$$

- **GENO+VAX** tier (files prefixed `NAKED_..._NoStringency.stan` / `SIR_..._GENO_VAX.stan`): vaccination is reintroduced (using $S_t^{*}$ from the vaccination deduction step), but $\alpha_{str}$ remains excluded:

$$
\lambda_t^{\text{GENO+VAX}} = \beta_0 \exp(\alpha_\nu \Theta_t) \cdot \frac{I_t}{N}
$$

- **FULL** tier: both vaccination and government stringency are included, matching the transmission model defined above:

$$
\lambda_t^{\text{FULL}} = \beta_0 \exp\big(\min[\alpha_\nu \Theta_t + \alpha_{str} Z_t,\ 10]\big) \cdot \frac{I_t}{N}
$$

Comparing $\alpha_\nu$ across these three fits shows how much of the genomic coupling effect is confounded with policy and vaccination timing versus attributable to variant pressure alone.

---

## Results & Discussion

Across 485 total model runs spanning eight regions (Africa, Asia, Europe, North America, South America, Oceania, and three global replicates), convergence patterns varied systematically with model complexity and intervention tier. Tier 1 (GENO-only) SIRS models converged more reliably than SEIRS counterparts, likely because removing vaccination and stringency covariates forces the genomic forcing term to account for both the rise and fall of case waves, straining identifiability in the additional exposed compartment. At Tier 3 (FULL), the more complex Erlang and memory-based models tended to converge better, consistent with the smoother posterior geometry these structures provide for the NUTS sampler.

The genomic coupling coefficient $\alpha_\nu$ was positive across most regions and tiers, indicating that genomic volatility signals amplify transmission, with notable exceptions in Oceania (near zero, credible interval spanning zero) and shifts to negative values in Asia once vaccination was added. Tier 1 models generally produced smaller $\alpha_\nu$ estimates than higher tiers, since the genomic term must absorb variance later explained by policy and vaccination covariates once those are included.

The genomic memory decay parameter $\delta$ (and its implied half-life) varied widely across regions and tiers, with no universal decay rate. Global and North America replicates showed memory half-lives shrinking to 22-25 days once vaccination was accounted for, while Europe's Tier 3 estimate carried a very wide credible interval (17.8-555.9 days), suggesting possible confounding. Africa and South America showed comparatively stable short-term memory across tiers, while $\delta$ became less identifiable at Tier 3 in several regions when estimated jointly with $\alpha_{str}$.

### Limitations

- Discrete weekly time steps introduce truncation error, particularly near rapidly changing peaks
- Fixed values for $\omega$, $\sigma$, and $\gamma$ (waning, incubation, recovery rates) were taken from early-pandemic estimates and may not hold across all variant periods
- The partial immunity fraction $\eta$ was fixed at 0.5 rather than estimated

### Future Work

- Estimate the partial immunity fraction $\eta$ rather than fixing it
- Allow disease parameters to vary across distinct variant periods
- Incorporate inter-regional mobility using open-source airline passenger flow data to model cross-continental variant spread
- Further examine the relationship between genomic change and the memory kernel structure

---

## Repository Structure

```text
bayesian-stochastic-epi/
├── data/
│   ├── model_data_agg_backup.csv
│   └── modeldataaggweekly.RDS
├── R/
│   └── Lucas Anderson - Ian McArthur - Bio 415 Project Analysis Pipeline.R
├── stan/
│   ├── GENO_M2_SEIRS.stan
│   ├── SEIR_Erlang_waning.stan
│   └── ... (all 24 Stan model implementations)
├── results/
│   ├── figures/
│   ├── posterior_summaries/
│   └── diagnostics/
├── paper/
│   └── (full write-up, PDF/LaTeX source)
└── .gitignore
```

---

## Data Sources

- Weekly SARS-CoV-2 case surveillance data
- Viral genomic volatility signals
- Vaccination rollout data
- Government stringency indices

Models are fit across Africa, Asia, Europe, North America, South America, Oceania, and global replicates.

---

## Future Extensions

- Estimating immunity parameters rather than fixing them
- Allowing disease parameters to vary across variant periods
- Incorporating regional mobility and travel networks
- Extending genomic forcing into non-Markovian memory kernels
- Linking latent epidemic states with genomic surveillance methods

---

## Citation

If you use this framework, please cite the accompanying paper (see [`paper/`](./paper)) and this repository.
