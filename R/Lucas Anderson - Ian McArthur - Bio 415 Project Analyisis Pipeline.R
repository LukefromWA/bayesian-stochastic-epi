###########################
#### Analysis Pipeline ####
###########################

library(cmdstanr)
library(rstan)
library(loo)
library(dplyr)
library(stringr)
library(tibble)
library(tidyr)

options(mc.cores = 4)

base_dir <- "SET_YOUR_BASE_PATH"
stan_dir <- "SET_YOUR_STAN_FILES_PATH"
output_dir <- "SET_YOUR_OUTPUT_PATH"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

safe_loo <- function(fit) {
  tryCatch({
    log_lik_mat <- loo::extract_log_lik(fit, merge_chains = FALSE)
    if (any(!is.finite(log_lik_mat))) return("FAILED_INSTABILITY")
    r_eff <- loo::relative_eff(exp(log_lik_mat), cores = 1)
    loo::loo(log_lik_mat, r_eff = r_eff, cores = 1)
  }, error = function(e) "LOO_CRASHED")
}

sum_erlang <- function(ext, p_name, s_name) {
  if (p_name %in% names(ext)) arr <- ext[[p_name]]
  else if (s_name %in% names(ext)) arr <- ext[[s_name]]
  else return(NULL)
  apply(arr, c(1, 2), sum)
}

make_init <- function(m_key, emp_i) {
  is_seir <- grepl("M2|M4|M6|M8", m_key)
  is_erlang <- grepl("M3|M4|M5|M6|M7|M8", m_key)
  i_low <- max(1e-8, emp_i * 0.50)
  i_up <- min(0.05, emp_i * 2.00)
  e_low <- max(1e-8, emp_i * 1.50)
  e_up <- min(0.05, emp_i * 3.00)
  function() {
    i <- list(
      beta0 = exp(rnorm(1, mean = log(0.40), sd = 0.20)),
      s0_frac = runif(1, 0.20, 0.45),
      i0_frac = runif(1, i_low, i_up),
      alpha_nu = rnorm(1, 0, 0.03),
      alpha_str = rnorm(1, 0, 0.03),
      inv_od = runif(1, 10.0, 20.0)
    )
    if (is_seir) i$e0_frac <- runif(1, e_low, e_up)
    if (is_erlang) i$delta <- runif(1, 0.15, 0.40)
    i
  }
}

stan_files <- c(
  "M2_NAKED" = "NAKED_M2_SEIRS_NoStringency.stan",
  "M4_NAKED" = "NAKED_M4_SEIRS_Memory_NoStringency.stan",
  "M6_NAKED" = "NAKED_M6_Erlang_SEIRS_NoStringency.stan",
  "M8_NAKED" = "NAKED_M8_Erlang_SEIRS_Immfrac_NoStringency.stan",
  "M2_GENO" = "GENO_M2_SEIRS.stan",
  "M4_GENO" = "GENO_M4_SEIRS_Memory.stan",
  "M6_GENO" = "GENO_M6_Erlang_SEIRS.stan",
  "M8_GENO" = "GENO_M8_Erlang_SEIRS_Immfrac.stan"
)

compiled_models <- setNames(
  lapply(names(stan_files), function(k) {
    cmdstan_model(file.path(stan_dir, stan_files[k]))
  }),
  names(stan_files)
)

jobs <- data.frame(Region = "Appendix", rep_num = 1, m_key = names(stan_files), stringsAsFactors = FALSE)

for (i in seq_len(nrow(jobs))) {
  m_key <- jobs$m_key[i]
  rep_num <- jobs$rep_num[i]
  res_key <- sprintf("Appendix_%s", m_key)
  rds_path <- file.path(output_dir, paste0(res_key, "_APPENDIX_RERUN.rds"))
  
  if (file.exists(rds_path)) next
  
  target_pop <- 8.15e9
  target_rho <- 0.30
  emp_i <- 0.001
  
  stan_dat <- list(
    n_weeks = 156,
    cases = rep(100, 156),
    nu_signal = rep(0, 156),
    nu_lag = rep(0, 156),
    stringency_z = rep(0, 156),
    vax_prop = rep(0, 156),
    pop = target_pop,
    rho = target_rho,
    gamma = 1/7,
    omega = 1/180
  )
  
  if (grepl("M2|M4|M6|M8", m_key)) stan_dat$sigma <- 1/5
  if (grepl("M5|M6|M7|M8", m_key)) stan_dat$k_erlang <- 4L
  if (grepl("M7|M8", m_key)) stan_dat$immfrac <- 0.5
  
  t0 <- proc.time()["elapsed"]
  
  fit_cmd <- try(compiled_models[[m_key]]$sample(
    data = stan_dat,
    seed = 200L + (i * 11L),
    init = make_init(m_key, emp_i),
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 1250,
    iter_sampling = 1250,
    adapt_delta = 0.95,
    max_treedepth = 13
  ))
  
  elapsed <- round((proc.time()["elapsed"] - t0) / 60, 1)
  
  if (!inherits(fit_cmd, "try-error")) {
    fit <- rstan::read_stan_csv(fit_cmd$output_files())
    ext <- rstan::extract(fit)
    sum_df <- as.data.frame(rstan::summary(fit)$summary) %>% rownames_to_column("parameter")
    
    E_post <- if (grepl("M6|M8", m_key)) sum_erlang(ext, "E_sub", "E")
    else if (grepl("M2|M4", m_key)) ext$E
    else NULL
    
    I_post <- if (grepl("M5|M6|M7|M8", m_key)) sum_erlang(ext, "I_sub", "I")
    else ext$I
    
    saveRDS(list(
      summary = sum_df,
      loo_result = safe_loo(fit),
      posteriors = list(S = ext$S, E = E_post, I = I_post, R = ext$R, expected_cases = ext$expected_cases),
      full_draws = ext,
      metadata = list(region = "Appendix", model = m_key, pop = target_pop)
    ), rds_path)
    
    rm(fit_cmd, fit, ext, sum_df, E_post, I_post); gc(); gc()
  } else {
    gc()
  }
}