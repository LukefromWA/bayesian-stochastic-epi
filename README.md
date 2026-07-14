# Bayesian Stochastic State-Space Epidemiological Models

A Bayesian state-space compartmental modeling framework implemented in **Stan** and **R**. This repository contains a 24-model ablation pipeline designed to evaluate how different biological assumptions, viral genomic signals, vaccination effects, and non-pharmaceutical interventions influence inferred epidemic dynamics.

The framework models SARS-CoV-2 transmission across multiple regions while increasing biological complexity through different compartment structures, disease progression assumptions, immunity mechanisms, and external forcing mechanisms.

Developed by **Lucas Anderson and Ian McArthur** for the **Bio 415 Project Analysis Pipeline**.

---

# Project Overview

Traditional SIR and SEIR models are widely used for studying infectious disease dynamics but often rely on simplifying assumptions, including constant transmission rates and exponentially distributed transition times.

This project extends stochastic compartmental models by incorporating:

- Viral genomic signals as time-dependent transmission drivers
- Genomic memory effects through latent forcing states
- Vaccination-driven changes in susceptibility
- Non-pharmaceutical intervention effects through government stringency
- Erlang-distributed latent and infectious periods
- Partial immunity and reinfection pathways
- Bayesian parameter estimation with uncertainty quantification

The goal of this framework is to evaluate how increasing biological realism changes:

- Model convergence
- Parameter inference
- Transmission dynamics
- Predictive performance

---

# Scientific Question

How do assumptions about biological realism and external forcing affect inference of epidemic transmission dynamics?

This framework evaluates:

- Does viral genomic variation improve transmission inference?
- Does a latent genomic memory state better represent evolutionary pressure than a fixed lag?
- Do Erlang distributed transition times improve model behavior compared to exponential assumptions?
- How do vaccination and policy interventions change inferred genomic coupling?

---

# 24-Model Factorial Design

The model suite follows an **8 × 3 factorial experimental design**.

Eight biological model structures are evaluated across three intervention tiers:

| Tier | Description |
|---|---|
| **GENO** | Viral genomic forcing only |
| **NAKED** | Viral genomic forcing + vaccination effects |
| **FULL** | Viral genomic forcing + vaccination + government stringency |

This design allows individual biological and policy mechanisms to be evaluated through structured model ablation.

---

# Biological Model Structures

| Model | Structure | Description |
|---|---|---|
| M2 | Lagged forcing SIRS/SEIRS | Baseline model using pre-lagged genomic forcing |
| M4 | Memory forcing SIRS/SEIRS | Adds a latent genomic memory state |
| M6 | Erlang SIRS/SEIRS | Replaces exponential transitions with Erlang distributed dwell times |
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

All models are implemented in **Stan** and fitted using **CmdStanR** with Hamiltonian Monte Carlo through the **No-U-Turn Sampler (NUTS)**.

The framework performs:

- Bayesian posterior estimation
- Parameter uncertainty quantification
- Posterior predictive evaluation
- Structural model comparison
- Convergence assessment

Models were evaluated using:

- $\hat{R} < 1.1$
- Effective sample size
- Posterior diagnostics

---

# Model Components

## Genomic Memory Forcing

Higher complexity models replace a fixed genomic lag with a latent memory state.

The memory process is defined as:

```math
A_t=(1-\delta)A_{t-1}+\nu_{signal,t}
```

where:

- $A_t$ represents accumulated genomic pressure
- $\delta$ controls memory decay
- $\nu_{signal,t}$ represents the standardized genomic volatility signal

The transmission rate is modeled as:

```math
\beta_t=\beta_0\exp(\alpha_\nu\Theta_t+\alpha_{str}Z_t)
```

where:

- $\Theta_t$ represents the active genomic forcing signal
- $Z_t$ represents standardized government stringency
- $\alpha_\nu$ estimates genomic coupling strength
- $\alpha_{str}$ estimates intervention effects

For lower complexity models:

```math
\Theta_t=\nu_{lag,t}
```

For memory models:

```math
\Theta_t=A_t
```

---

# Disease Progression Models

## Erlang Distributed Dwell Times

Higher complexity models replace memoryless exponential transitions with Erlang distributed progression.

The Erlang structure divides compartments into sequential substages:

```math
I_1 \rightarrow I_2 \rightarrow I_3 \rightarrow I_4 \rightarrow R
```

This allows the model to represent more realistic infectious and latent period distributions.

---

# Observation Model

Observed weekly cases are modeled using an overdispersed Negative Binomial observation process:

```math
y_t \sim NegBin2(\mu_t,\phi^{-1})
```

where:

- $y_t$ represents observed weekly cases
- $\mu_t$ represents expected reported infections
- $\phi$ controls overdispersion

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

# Data Sources

The framework integrates:

- Weekly SARS-CoV-2 case surveillance data
- Viral genomic volatility signals
- Vaccination rollout data
- Government stringency indicators

Models are evaluated across:

- Africa
- Asia
- Europe
- North America
- South America
- Oceania
- Global replicates

---

# Future Extensions

Potential extensions include:

- Estimating immunity parameters rather than fixing them
- Allowing disease parameters to vary across variant periods
- Incorporating regional mobility and travel networks
- Extending genomic forcing into non-Markovian memory kernels
- Linking latent epidemic states with genomic surveillance methods

This project provides a framework for studying how hidden biological processes and external interventions influence observed epidemic dynamics through stochastic state-space modeling.