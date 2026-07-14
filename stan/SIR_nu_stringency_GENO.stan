// SIR_nu_stringency_GENO.stan
// Sensitivity: Genomic forcing only. No stringency. No vaccination.
data {
  int<lower=1> n_weeks;
  array[n_weeks] int<lower=0> cases;
  vector[n_weeks] nu_lag;
  real<lower=0> pop;
  real<lower=0, upper=1> rho;
  real<lower=0> gamma;
  real<lower=0> omega;
}
parameters {
  real<lower=0> beta0;
  real alpha_nu;
  real<lower=1e-8, upper=0.1> i0_frac;
  real<lower=0.1, upper=0.99> s0_frac;
  real<lower=0.01> inv_od;
}
transformed parameters {
  vector[n_weeks] expected_cases;
  vector[n_weeks+1] S;
  vector[n_weeks+1] I;
  vector[n_weeks+1] R;

  S[1] = pop * s0_frac;
  I[1] = pop * i0_frac;
  R[1] = pop - S[1] - I[1];

  for (t in 1:n_weeks) {
    real forcing = fmin(fmax(alpha_nu * nu_lag[t], -10), 10);
    real beta_t  = beta0 * exp(forcing);
    real foi     = beta_t * I[t] / pop;
    real dSI     = S[t] * (1 - exp(-foi));
    real dIR     = I[t] * (1 - exp(-gamma));
    real dRS     = R[t] * (1 - exp(-omega));

    S[t+1] = fmax(S[t] - dSI + dRS, 1.0);
    I[t+1] = fmax(I[t] + dSI - dIR, 1.0);
    R[t+1] = fmax(R[t] + dIR - dRS, 1.0);

    expected_cases[t] = fmax(dSI * rho, 0.1);
  }
}
model {
  beta0    ~ lognormal(0, 0.5);
  alpha_nu ~ normal(0, 0.3);
  inv_od   ~ gamma(2, 0.1);
  s0_frac  ~ beta(3, 3);
  i0_frac  ~ beta(1, 80);
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
