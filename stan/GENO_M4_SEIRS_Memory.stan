// GENO_M4_SEIRS_Memory.stan
// Pure Genomics/Biology with smooth Nu signal. No Stringency, No Vaccination.
data {
  int<lower=1> n_weeks;
  array[n_weeks] int<lower=0> cases;
  vector[n_weeks] nu_signal;
  real<lower=0> pop;
  real<lower=0, upper=1> rho;
  real<lower=0> sigma;
  real<lower=0> gamma;
  real<lower=0> omega;
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
  vector[n_weeks+1] E;
  vector[n_weeks+1] I;
  vector[n_weeks+1] R;
  vector[n_weeks] A;

  S[1] = pop * s0_frac;
  E[1] = pop * e0_frac;
  I[1] = pop * i0_frac;
  R[1] = pop - S[1] - E[1] - I[1];
  A[1] = nu_signal[1];

  for (t in 1:n_weeks) {
    if (t > 1)
      A[t] = (1 - delta) * A[t-1] + nu_signal[t];

    real forcing = fmin(fmax(alpha_nu * A[t], -10), 10);
    real beta_t  = beta0 * exp(forcing);
    real foi     = beta_t * I[t] / pop;
    
    real dSE     = S[t] * (1 - exp(-foi));
    real dEI     = E[t] * (1 - exp(-sigma));
    real dIR     = I[t] * (1 - exp(-gamma));
    real dRS     = R[t] * (1 - exp(-omega));

    S[t+1] = fmax(S[t] - dSE + dRS, 1.0);
    E[t+1] = fmax(E[t] + dSE - dEI, 1.0);
    I[t+1] = fmax(I[t] + dEI - dIR, 1.0);
    R[t+1] = fmax(R[t] + dIR - dRS, 1.0);

    expected_cases[t] = fmax(dEI * rho, 1e-5);
  }
}
model {
  beta0    ~ lognormal(log(0.40), 0.20);
  alpha_nu ~ normal(0, 0.5);
  delta    ~ beta(2, 5);
  i0_frac  ~ beta(1, 99);
  e0_frac  ~ beta(1, 99);
  s0_frac  ~ beta(2, 2);
  inv_od   ~ exponential(0.1);

  cases ~ neg_binomial_2(expected_cases, inv_od);
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
