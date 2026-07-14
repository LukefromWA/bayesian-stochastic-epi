# Bayesian Stochastic State-Space Epidemiological Models

A Bayesian state-space compartmental modeling framework implemented in **Stan** and **R**. This repository contains a 24-model ablation pipeline designed to evaluate how biological assumptions, viral genomic signals, vaccination effects, and non-pharmaceutical interventions influence inferred epidemic dynamics.

The framework models SARS-CoV-2 transmission across multiple global regions while increasing biological complexity through changes in compartment structure, disease progression assumptions, immunity mechanisms, and external forcing mechanisms.

Developed by **Lucas Anderson and Ian McArthur** for the **Bio 415 Project Analysis Pipeline**.

---

# Project Overview

SIR and SEIR models are widely used to study infectious disease dynamics but often rely on simplifying assumptions such as constant transmission rates and exponentially distributed transition times.

This project extends stochastic compartmental models by incorporating:

- Viral genomic signals as time-dependent transmission drivers
- Genomic memory effects through latent forcing states
- Vaccination-driven changes in susceptibility
- Government stringency effects as external intervention forcing
- Erlang-distributed latent and infectious periods
- Partial immunity and reinfection pathways
- Bayesian parameter estimation with uncertainty quantification

The goal of this framework is to evaluate how increasing biological realism affects:

- Model convergence
- Parameter inference
- Transmission dynamics
- Predictive performance

---

# Scientific Question

How do biological assumptions and external forcing mechanisms affect inference of epidemic transmission dynamics?

This framework evaluates:

- Whether viral genomic variation improves transmission inference
- Whether genomic memory states better represent evolutionary pressure than fixed lags
- Whether Erlang distributed dwell times improve epidemic state representation
- How vaccination and policy interventions alter inferred genomic coupling

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

# Biological Model Hierarchy

The model ladder progressively increases biological realism.

| Model | Structure | Description |
|---|---|---|
| M1 | SIRS | Baseline susceptible-infected-recovered-susceptible model |
| M2 | SEIRS | Adds an exposed compartment representing latent infection |
| M3 | SIRS + Genomic Memory | Adds accumulated genomic forcing |
| M4 | SEIRS + Genomic Memory | Adds accumulated genomic forcing to SEIRS |
| M5 | Erlang SIRS | Replaces infectious period with Erlang distributed stages |
| M6 | Erlang SEIRS | Adds Erlang distributed exposed and infectious periods |
| M7 | Erlang SIRS + Partial Immunity | Adds reinfection pathway from partially immune recovered individuals |
| M8 | Erlang SEIRS + Partial Immunity | Adds reinfection through exposed pathway |

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

# Model Assumptions

## M1: SIRS Model

The population is assumed to be closed, with no births, deaths, or movement between regions.

Individuals transition between:

```text
Susceptible → Infected → Recovered → Susceptible
```

New infections occur through contact between susceptible and infected individuals.

The transmission rate is affected by external forcing signals. Infected individuals recover at a fixed rate, while recovered individuals lose immunity and return to susceptibility.

Vaccination removes susceptible individuals and transfers them into the recovered compartment.

---

## M2: SEIRS Model

The SIRS model is extended by adding an exposed compartment.

Individuals transition through:

```text
Susceptible → Exposed → Infected → Recovered → Susceptible
```

The exposed compartment represents the latent period between infection and infectiousness.

---

## M3 / M4: Genomic Memory Models

Instead of using only the previous week's genomic volatility signal, genomic pressure is assumed to accumulate and decay over time.

The latent genomic memory state is:

```math
A_t=(1-\delta)A_{t-1}+\nu_{signal,t}
```

where:

- $A_t$ is the accumulated genomic forcing state
- $\delta$ controls memory decay
- $\nu_{signal,t}$ is the standardized genomic volatility signal

---

## M5 / M6: Erlang Dwell-Time Models

The infectious and latent periods are modeled using Erlang distributions with:

```math
k=4
```

Instead of a single exponential compartment, individuals move through sequential substages.

Example infectious progression:

```text
I1 → I2 → I3 → I4 → R
```

This allows the model to represent more realistic variation in disease progression times.

M5 applies Erlang progression to SIRS infectious periods.

M6 applies Erlang progression to both SEIRS exposed and infectious periods.

---

## M7 / M8: Partial Immunity Models

The final model tier introduces partial immunity after recovery.

Instead of all recovered individuals returning directly to susceptibility after immunity wanes, a fraction of individuals retain partial immunity and follow alternative reinfection pathways.

### M7: Erlang SIRS + Partial Immunity

Partially immune individuals bypass full susceptibility:

```math
R \rightarrow I_1
```

### M8: Erlang SEIRS + Partial Immunity

Partially immune individuals return through the exposed pathway:

```math
R \rightarrow E_1
```

These models evaluate how incomplete immunity and reinfection affect epidemic dynamics.

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

# Transmission Model

The transmission rate is modeled as:

```math
\beta_t=\beta_0\exp(\alpha_\nu\Theta_t+\alpha_{str}Z_t)
```

where:

- $\beta_0$ is the baseline transmission rate
- $\Theta_t$ is the genomic forcing signal
- $Z_t$ is the standardized government stringency index
- $\alpha_\nu$ measures genomic coupling
- $\alpha_{str}$ measures intervention effects

For lagged forcing models:

```math
\Theta_t=\nu_{lag,t}
```

For memory models:

```math
\Theta_t=A_t
```

---

# Observation Model

Observed weekly cases are modeled using an overdispersed Negative Binomial observation model:

```math
y_t \sim NegBin2(\mu_t,\phi^{-1})
```

where:

- $y_t$ represents observed cases
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
- Government stringency indices

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