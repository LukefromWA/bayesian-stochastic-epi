# Bayesian Stochastic State-Space Epidemiological Models

A Bayesian state-space compartmental modeling framework implemented in **Stan** and **R**, evaluating how biological assumptions, viral genomic signals, vaccination effects, and non-pharmaceutical interventions influence inferred SARS-CoV-2 transmission dynamics across global regions.

**Authors:** Lucas (Luke) Anderson, Ian McArthur
**Course:** Bio 415 Project Analysis Pipeline

---

## Contributions

- **Luke Anderson:** Model design and biological/statistical formulation, Stan/R implementation, simulation and Bayesian inference pipeline, MCMC diagnostics, data processing, results visualization, drafting of the Results section, and portions of the Methods / Model Ladder section.
- **Ian McArthur:** Manuscript writing for the Abstract, Introduction, Discussion, and Model description sections, and contributions to the project presentation.

Full methodological and results write-up is provided in [`paper/`](./paper) (add PDF/LaTeX source here).

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

| Model | Structure | Description |
|---|---|---|
| M2 | SEIRS (lagged genomic forcing) | Exposed compartment + lagged genomic signal |
| M2 | SIRS (lagged genomic forcing) | Lagged genomic signal, no exposed compartment |
| M4 | SEIRS + genomic memory | Accumulated, decaying genomic forcing state |
| M4 | SIRS + genomic memory | Accumulated, decaying genomic forcing state |
| M6 | SEIRS + Erlang dwell times | Multi-stage exposed and infectious periods (k = 4) |
| M6 | SIRS + Erlang dwell times | Multi-stage infectious period (k = 4) |
| M8 | SEIRS + Erlang + partial immunity | Erlang dwell times + reinfection via exposed pathway |
| M8 | SIRS + Erlang + partial immunity | Erlang dwell times + reinfection via infectious pathway |

---

## 24-Model Matrix

| Biological Complexity | Structure | GENO | GENO+VAX | FULL |
|---|---|---|---|---|
| M2 (Lagged) | SEIRS | `GENO_M2_SEIRS.stan` | `NAKED_M2_SEIRS_NoStringency.stan` | `SEIR_nu_stringency.stan` |
| M2 (Lagged) | SIRS | `SIR_nu_stringency_GENO.stan` | `SIR_nu_stringency_GENO_VAX.stan` | `SIR_nu_stringency.stan` |
| M4 (Memory) | SEIRS | `GENO_M4_SEIRS_Memory.stan` | `NAKED_M4_SEIRS_Memory_NoStringency.stan` | `SEIR_nu_smooth.stan` |
| M4 (Memory) | SIRS | `SIR_nu_smooth_GENO.stan` | `SIR_nu_smooth_GENO_VAX.stan` | `SIR_nu_smooth.stan` |
| M6 (Erlang) | SEIRS | `GENO_M6_Erlang_SEIRS.stan` | `NAKED_M6_Erlang_SEIRS_NoStringency.stan` | `SEIR_Erlang.stan` |
| M6 (Erlang) | SIRS | `SIR_Erlang_GENO.stan` | `SIR_Erlang_GENO_VAX.stan` | `SIR_Erlang.stan` |
| M8 (Erlang + Waning) | SEIRS | `GENO_M8_Erlang_SEIRS_Immfrac.stan` | `NAKED_M8_Erlang_SEIRS_Immfrac_NoStringency.stan` | `SEIR_Erlang_waning.stan` |
| M8 (Erlang + Waning) | SIRS | `SIR_Erlang_waning_GENO.stan` | `SIR_Erlang_waning_GENO_VAX.stan` | `SIR_Erlang_waning.stan` |

---

## Model Assumptions

### Base SIRS / SEIRS structure

The population is closed (no births, deaths, or migration).

```text
SIRS:  Susceptible -> Infected -> Recovered -> Susceptible
SEIRS: Susceptible -> Exposed -> Infected -> Recovered -> Susceptible
```

New infections arise through contact between susceptible and infectious individuals. The transmission rate is modulated by external forcing signals; recovered individuals lose immunity and return to susceptibility. Vaccination moves susceptible individuals directly into the recovered compartment.

### Genomic memory (M4)

Genomic pressure accumulates and decays over time rather than depending only on the prior week's signal:

```math
A_t = (1-\delta) A_{t-1} + \nu_{signal,t}
```

where \(A_t\) is the accumulated genomic forcing state, \(\delta\) is the memory decay rate, and \(\nu_{signal,t}\) is the standardized weekly genomic volatility signal.

### Erlang dwell times (M6, M8)

Latent and infectious periods are modeled as Erlang-distributed with \(k = 4\) sequential sub-stages instead of a single exponential compartment:

```text
I1 -> I2 -> I3 -> I4 -> R
```

M6 applies this to SIRS infectious periods and SEIRS exposed+infectious periods. M8 extends M6 with partial immunity.

### Partial immunity (M8)

A fraction of recovered individuals retain partial immunity and re-enter the infection process directly rather than passing through full susceptibility:

```math
R \rightarrow I_1 \quad \text{(SIRS pathway)}
```
```math
R \rightarrow E_1 \quad \text{(SEIRS pathway)}
```

---

## Bayesian Inference Framework

All models are implemented in **Stan** and fitted with **CmdStanR** using Hamiltonian Monte Carlo via the No-U-Turn Sampler (NUTS). The pipeline performs Bayesian posterior estimation, parameter uncertainty quantification, posterior predictive evaluation, structural model comparison, and convergence assessment.

### Convergence criteria

- Rank-normalized split-\(\hat R < 1.01\)
- Bulk and tail effective sample size (target > 400 across chains)
- Divergent transition count
- Maximum treedepth saturation
- Energy / E-BFMI diagnostics
- Posterior predictive checks against observed weekly case trajectories

---

## Transmission Model

```math
\beta_t = \beta_0 \exp(\alpha_\nu \Theta_t + \alpha_{str} Z_t)
```

- \(\beta_0\): baseline transmission rate
- \(\Theta_t\): genomic forcing signal (lagged or memory-based)
- \(Z_t\): standardized government stringency index
- \(\alpha_\nu\): genomic coupling coefficient
- \(\alpha_{str}\): intervention effect coefficient

Lagged forcing models set \(\Theta_t = \nu_{lag,t}\); memory models set \(\Theta_t = A_t\).

---

## Observation Model

```math
y_t \sim \text{NegBin2}(\mu_t, \phi^{-1})
```

- \(y_t\): observed weekly cases
- \(\mu_t\): expected observed cases, computed as new-infection flow scaled by the region-specific reporting rate \(\rho\) (e.g. \(\Delta_{SI}\cdot\rho\) for the SIRS backbone, \(\Delta_{EI}\cdot\rho\) for exponential SEIRS, \(\Delta_{E_4}\cdot\rho\) for Erlang SEIRS)
- \(\phi^{-1}\): inverse overdispersion parameter, prior \(\phi^{-1}\sim\text{Gamma}(2,0.1)\), accounting for the overdispersed clumping typical of COVID-19 case data

---

## Results & Discussion

Across 485 total model runs spanning eight regions (Africa, Asia, Europe, North America, South America, Oceania, and three global replicates), convergence patterns varied systematically with model complexity and intervention tier. Tier 1 (GENO-only) SIRS models converged more reliably than SEIRS counterparts, likely because removing vaccination and stringency covariates forces the genomic forcing term to account for both the rise and fall of case waves, straining identifiability in the additional exposed compartment. At Tier 3 (FULL), the more complex Erlang and memory-based models tended to converge better, consistent with the smoother posterior geometry these structures provide for the NUTS sampler.

The genomic coupling coefficient \(\alpha_\nu\) was positive across most regions and tiers, indicating that genomic volatility signals amplify transmission, with notable exceptions in Oceania (near zero, credible interval spanning zero) and shifts to negative values in Asia once vaccination was added. Tier 1 models generally produced smaller \(\alpha_\nu\) estimates than higher tiers, since the genomic term must absorb variance later explained by policy and vaccination covariates once those are included.

The genomic memory decay parameter \(\delta\) (and its implied half-life) varied widely across regions and tiers, with no universal decay rate. Global and North America replicates showed memory half-lives shrinking to 22-25 days once vaccination was accounted for, while Europe's Tier 3 estimate carried a very wide credible interval (17.8-555.9 days), suggesting possible confounding. Africa and South America showed comparatively stable short-term memory across tiers, while \(\delta\) became less identifiable at Tier 3 in several regions when estimated jointly with \(\alpha_{str}\).

### Limitations

- Discrete weekly time steps introduce truncation error, particularly near rapidly changing peaks
- Fixed values for \(\omega\), \(\sigma\), and \(\gamma\) (waning, incubation, recovery rates) were taken from early-pandemic estimates and may not hold across all variant periods
- The partial immunity fraction \(\eta\) was fixed at 0.5 rather than estimated

### Future Work

- Estimate the partial immunity fraction \(\eta\) rather than fixing it
- Allow disease parameters to vary across distinct variant periods
- Incorporate inter-regional mobility using open-source airline passenger flow data to model cross-continental variant spread
- Further examine the relationship between genomic change and the memory kernel structure

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