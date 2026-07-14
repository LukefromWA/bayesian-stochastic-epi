// NAKED_M8_Erlang_SEIRS_Immfrac_NoStringency.stan
// Identical to full M8 (SEIR_Erlang_waning.stan) except:
//   - stringency_z covariate removed from data and forcing
//   - alpha_str parameter removed
data {
  int<lower=1> n_weeks;
  array[n_weeks] int<lower=0> cases;
  vector[n_weeks] nu_signal;
  vector[n_weeks] vax_prop;
  real<lower=0> pop;
  real<lower=0, upper=1> rho;
  real<lower=0> sigma;
  real<lower=0> gamma;
  real<lower=0> omega;
  real<lower=0, upper=1> immfrac;
  int<lower=1> k_erlang;
}
parameters {
  real<lower=0> beta0;
  real alpha_nu;
  real<lower=0, upper=1> delta;
  real<lower=1e-8, upper=0.1> i0_frac;
  real<lower=1e-8, upper=0.1> e0_frac;
  real<lower=0.1, upper=0.99> s0_frac;
  real<lower=0.01> inv_od;
}
transformed parameters {
  vector[n_weeks] expected_cases;
  vector[n_weeks+1] S;
  matrix[n_weeks+1, k_erlang] E;
  matrix[n_weeks+1, k_erlang] I;
  vector[n_weeks+1] R;
  vector[n_weeks] A;

  real sigma_k = sigma * k_erlang;
  real gamma_k = gamma * k_erlang;

  S[1] = pop * s0_frac;
  for (k in 1:k_erlang) {
    E[1, k] = pop * e0_frac / k_erlang;
    I[1, k] = pop * i0_frac / k_erlang;
  }
  R[1] = pop - S[1] - sum(E[1,]) - sum(I[1,]);
  A[1] = nu_signal[1];

  for (t in 1:n_weeks) {
    if (t > 1)
      A[t] = (1 - delta) * A[t-1] + nu_signal[t];

    real forcing   = fmin(fmax(alpha_nu * A[t], -10), 10);
    real beta_t    = beta0 * exp(forcing);
    real foi       = beta_t * sum(I[t,]) / pop;
    real dVax      = fmin(S[t], vax_prop[t] * pop);
    real S_post    = S[t] - dVax;
    real dSE       = S_post * (1 - exp(-foi));
    real dRS_total = R[t] * (1 - exp(-omega));
    real dRS_to_S  = dRS_total * (1 - immfrac);
    real dRS_to_E  = dRS_total * immfrac;

    vector[k_erlang] dE;
    vector[k_erlang] dI;
    for (k in 1:k_erlang) {
      dE[k] = E[t, k] * (1 - exp(-sigma_k));
      dI[k] = I[t, k] * (1 - exp(-gamma_k));
    }

    S[t+1]    = fmax(S_post - dSE + dRS_to_S, 1.0);
    E[t+1, 1] = fmax(E[t, 1] + dSE + dRS_to_E - dE[1], 1.0);
    for (k in 2:k_erlang)
      E[t+1, k] = fmax(E[t, k] + dE[k-1] - dE[k], 1.0);
    I[t+1, 1] = fmax(I[t, 1] + dE[k_erlang] - dI[1], 1.0);
    for (k in 2:k_erlang)
      I[t+1, k] = fmax(I[t, k] + dI[k-1] - dI[k], 1.0);
    R[t+1] = fmax(R[t] + dI[k_erlang] - dRS_total + dVax, 1.0);

    expected_cases[t] = fmax(dE[k_erlang] * rho, 0.1);
  }
}
model {
  beta0    ~ lognormal(0, 0.5);
  alpha_nu ~ normal(0, 0.3);
  delta    ~ beta(2, 5);
  inv_od   ~ gamma(2, 0.1);
  s0_frac  ~ beta(3, 3);
  i0_frac  ~ beta(1, 80);
  e0_frac  ~ beta(1, 80);
  cases    ~ neg_binomial_2(expected_cases, inv_od);
}
generated quantities {
  vector[n_weeks] log_lik;
  array[n_weeks] int cases_sim;
  for (t in 1:n_weeks) {
    log_lik[t]   = neg_binomial_2_lpmf(cases[t] | expected_cases[t], inv_od);
    cases_sim[t] = (expected_cases[t] < 1e3 && inv_od > 1e-6)
                   ? neg_binomial_2_rng(expected_cases[t], inv_od) : -1;
  }
}
