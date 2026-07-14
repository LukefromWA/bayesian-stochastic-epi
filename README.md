# Bayesian Stochastic State-Space Epidemiological Models

A Bayesian state-space compartmental modeling framework implemented in **Stan** and **R**. This repository contains a 24-model ablation pipeline designed to evaluate how different biological assumptions, viral genomic signals, vaccination effects, and non-pharmaceutical interventions influence inferred epidemic dynamics.

The framework models SARS-CoV-2 transmission across multiple regions while increasing biological complexity through different compartment structures, disease progression assumptions, immunity mechanisms, and external forcing mechanisms.

Developed by **Lucas Anderson and Ian McArthur** for the **Bio 415 Project Analysis Pipeline**.

---

## Project Overview

Traditional SIR and SEIR models are powerful tools for studying infectious disease dynamics but often rely on simplifying assumptions such as constant transmission rates and exponentially distributed transition times.

This project extends stochastic compartmental models by incorporating:

- Viral genomic signals as time-dependent transmission drivers
- Genomic memory effects through latent forcing states
- Vaccination-driven changes in susceptibility
- Non-pharmaceutical intervention effects through government stringency
- Erlang-distributed latent and infectious periods
- Partial immunity and reinfection pathways
- Bayesian parameter estimation with uncertainty quantification

The goal of this framework is to evaluate how increasing biological realism changes parameter inference, model convergence, and predictive performance.

---

# 24-Model Factorial Design

The model suite follows an **8 × 3 factorial experimental design**.

Eight biological model structures are evaluated across three covariate intervention tiers:

| Tier | Description |
|---|---|
| **GENO** | Viral genomic forcing only |
| **NAKED** | Viral genomic forcing + vaccination effects |
| **FULL** | Viral genomic forcing + vaccination + government stringency |

This design allows the effect of individual biological and policy mechanisms to be evaluated through structured model ablation.

---

# Biological Model Structures

| Model | Structure | Description |
|---|---|---|
| M2 | Lagged forcing SIRS/SEIRS | Baseline model using pre-lagged genomic forcing |
| M4 | Memory forcing SIRS/SEIRS | Adds a latent genomic memory state |
| M6 | Erlang SIRS/SEIRS | Uses Erlang distributed latent and infectious periods |
| M8 | Erlang + waning immunity | Adds partial immunity and reinfection pathways |

---

# 24-Model Matrix

| Biological Complexity | Compartment Structure | GENO Tier | NAKED Tier | FULL Tier |
|---|---|---|---|---|
| M2 (Lagged) | SEIRS | GENO_M2_SEIRS.stan | NAKED_M2_SEIRS_NoStringency.stan | SEIR_nu_stringency.stan |
| M2 (Lagged) | SIRS | SIR_nu_stringency_GENO.stan | SIR_nu_stringency_GENO_VAX.stan | SIR_nu_stringency.stan |
| M4 (Continuous Memory) | SEIRS | GENO_M4_SEIRS_Memory.stan | NAKED_M4_SEIRS_Memory_NoStringency.stan | SEIR_nu_smooth.stan |
| M4 (Continuous Memory) | SIRS | SIR_nu_smooth_GENO.stan | SIR_nu_smooth_GENO_VAX.stan | SIR_nu_smooth.stan |
| M6 (Erlang Dwell Times) | SEIRS | GENO_M6_Erlang_SEIRS.stan | NAKED_M6_Erlang_SEIRS_NoStringency.stan | SEIR_Erlang.stan |
| M6 (Erlang Dwell Times) | SIRS | SIR_Erlang_GENO.stan | SIR_Erlang_GENO_VAX.stan | SIR_Erlang.stan |
| M8 (Erlang + Waning) | SEIRS | GENO_M8_Erlang_SEIRS_Immfrac.stan | NAKED_M8_Erlang_SEIRS_Immfrac_NoStringency.stan | SEIR_Erlang_waning.stan |
| M8 (Erlang + Waning) | SIRS | SIR_Erlang_waning_GENO.stan | SIR_Erlang_waning_GENO_VAX.stan | SIR_Erlang_waning.stan |

---

# Bayesian Inference Framework

All models are implemented in Stan and fitted using **CmdStanR** with Hamiltonian Monte Carlo through the **No-U-Turn Sampler (NUTS)**.

The framework performs:

- Bayesian posterior estimation
- Parameter uncertainty quantification
- Convergence assessment
- Posterior predictive evaluation
- Structural model comparison

Models were evaluated using:

- $\hat{R} < 1.1$
- Effective sample size
- Posterior diagnostics

---

# Model Components

## Genomic Memory Forcing

Higher complexity models replace a fixed genomic lag with a latent memory state:

$$
A_t=(1-\delta)A_{t-1}+\nu_{signal,t}
$$

where genomic volatility accumulates over time and decays according to the estimated memory parameter.

The transmission rate is modeled as:

$$
\beta_t=\beta_0\exp(\alpha_\nu\Theta_t+\alpha_{str}Z_t)
$$

where:

- $\Theta_t$ represents genomic forcing
- $Z_t$ represents government stringency
- $\alpha_\nu$ estimates genomic coupling strength
- $\alpha_{str}$ estimates intervention effects

---

# Repository Structure

```text
bayesian-stochastic-epi/
│
├── data/
│   ├── model_data_agg_backup.csv
│   └── modeldataaggweekly.RDS
│
├── R/
│   └── Lucas Anderson - Ian McArthur - Bio 415 Project Analysis Pipeline.R
│
├── stan/
│   ├── GENO_M2_SEIRS.stan
│   ├── SEIR_Erlang_waning.stan
│   └── ... (all 24 Bayesian model implementations)
│
├── results/
│   ├── figures/
│   ├── posterior summaries/
│   └── diagnostics/
│
└── .gitignore
```

---

# Future Extensions

Potential extensions include:

- Estimating immunity parameters instead of fixing them
- Allowing disease parameters to vary across variant periods
- Incorporating regional mobility and travel networks
- Extending memory forcing into non-Markovian transmission kernels
- Linking latent epidemic states with genomic surveillance methods

This project provides a framework for studying how hidden biological processes and external interventions influence observed epidemic dynamics through stochastic state-space modeling.