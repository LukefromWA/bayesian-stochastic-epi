# Bayesian Stochastic State-Space Epidemiological Models

A comprehensive, 24-model Bayesian state-space compartmental modeling framework implemented in Stan and R. This pipeline evaluates the explanatory power of viral genomic signals alongside vaccination rollouts and non-pharmaceutical behavioral interventions across multiple layers of biological complexity and policy constraints.

Developed by **Lucas Anderson** and **Ian McArthur** for the **Bio 415 Project Analysis Pipeline**.

---

## Project Overview

This repository uses an **8 x 3 factorial experimental design** to rigorously evaluate how varying biological complexity, transition assumptions, and intervention tiers affect the fit and predictive performance of stochastic epidemic models.

### The 24-Model Factorial Matrix

Our model suite is structured into eight distinct biological configurations evaluated across three covariate intervention tiers:

| Biological Complexity | Compartmental Structure | GENO Tier (Evolution Only) | NAKED Tier (Evolution + Vax) | FULL Tier (Evolution + Vax + Policy) |
| :--- | :--- | :--- | :--- | :--- |
| **M2 (Lagged)** | **SEIRS** | `GENO_M2_SEIRS.stan` | `NAKED_M2_SEIRS_NoStringency.stan` | `SEIR_nu_stringency.stan` |
| | **SIRS** | `SIR_nu_stringency_GENO.stan` | `SIR_nu_stringency_GENO_VAX.stan` | `SIR_nu_stringency.stan` |
| **M4 (Continuous Memory)** | **SEIRS** | `GENO_M4_SEIRS_Memory.stan` | `NAKED_M4_SEIRS_Memory_NoStringency.stan` | `SEIR_nu_smooth.stan` |
| | **SIRS** | `SIR_nu_smooth_GENO.stan` | `SIR_nu_smooth_GENO_VAX.stan` | `SIR_nu_smooth.stan` |
| **M6 (Erlang Dwell Times)** | **SEIRS** | `GENO_M6_Erlang_SEIRS.stan` | `NAKED_M6_Erlang_SEIRS_NoStringency.stan` | `SEIR_Erlang.stan` |
| | **SIRS** | `SIR_Erlang_GENO.stan` | `SIR_Erlang_GENO_VAX.stan` | `SIR_Erlang.stan` |
| **M8 (Erlang + Waning)** | **SEIRS** | `GENO_M8_Erlang_SEIRS_Immfrac.stan` | `NAKED_M8_Erlang_SEIRS_Immfrac_NoStringency.stan` | `SEIR_Erlang_waning.stan` |
| | **SIRS** | `SIR_Erlang_waning_GENO.stan` | `SIR_Erlang_waning_GENO_VAX.stan` | `SIR_Erlang_waning.stan` |

---

## Directory Structure

```text
bayesian-stochastic-epi/
├── data/
│   ├── model_data_agg_backup.csv     # Aggregated dataset backups
│   └── modeldataaggweekly.RDS        # Consolidated weekly surveillance datasets
├── R/
│   └── Lucas Anderson - Ian McArthur - Bio 415 Project Analyisis Pipeline.R # Execution pipeline
├── stan/
│   ├── GENO_M2_SEIRS.stan            # Lagged genomic-only SEIRS
│   ├── SEIR_Erlang_waning.stan       # Erlang waning full-tier model
│   └── ...                           # (All 24 compiled Bayesian scripts)
└── .gitignore                        # Standard protection against CmdStan binary leakages
