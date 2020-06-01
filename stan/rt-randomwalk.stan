functions {
  real get_at_time(int n, vector v, int t_min, int t) {
    int i = t - t_min + 1;
    if(i < 1 || i > n) {
      return 0.0;
    }
    return v[i];
  }
  
  real adjust_gamma_param(real x) {
    return 1e-10 + (1.0 - 1e-10) * x;
  }
  
  real beta_binomial_gamma_shape(real n, real p, real rate) {
    return adjust_gamma_param(n * p * rate);
  }
  
  real beta_binomial_gamma_rate(real n, real p, real dispersion) {
    return adjust_gamma_param((
      dispersion + p * (1.0 - p)
    ) / (
      (1.0 - p) * (dispersion + n * p * (1.0 - p))
    ));
  }
}

data {
  // Smallest (most negative) observation delay
  int delay_min;
  
  // Length of observation delay distribution
  int n_delay;
  
  // PMF of delay distribution
  simplex[n_delay] pmf_delay;
  
  // Smallest delay to infection
  int dt_inf_min;
  
  // Length of infection distribution ("infectiousness profile")
  int n_dt_inf;
  
  // PMF of infection distribution
  simplex[n_dt_inf] pmf_inf;
  
  // Initial time of observed time series
  int t_obs_min;
  
  // Length of observed, delayed time series
  int n_obs;
  
  // Delayed observations
  int<lower=0> y_obs[n_obs];
  
  // Initial timepoint for unobserved states
  int t_unobs_min;
  
  // Scale for prior on initial incidence
  real<lower=0> scale_x_unobs_init;
  
  // Scale for prior on initial Rt
  real<lower=0> scale_rt_init;
  
  // Scale for prior on SD for random walk on log(Rt)
  real<lower=0> scale_sd_dlogrt;
  
  // Observation probability
  real<lower=0, upper=1> p_obs;
  
  // Scale for prior on beta-binomial observation dispersion
  real<lower=0> scale_obs_dispersion;
}

transformed data {
  // Largest delay to infection
  int dt_inf_max = dt_inf_min + n_dt_inf - 1;
  
  // Last timepoint: shared between observations and unobserved states
  // (final timepoints will have wide uncertainty)
  int t_max = t_obs_min + n_obs - 1;
  
  // Number of unobserved states
  int n_unobs = t_max - t_unobs_min + 1;
  
  // Largest delay
  int delay_max = delay_min + n_delay - 1;
}

parameters {
  // Initial incidence
  real<lower=0> x_unobs_init;
  
  // Initial rt
  real<lower=0> rt_init;
  
  // SD for the random walk on log(Rt)
  real<lower=0> sd_dlogrt;
  
  // Changes in log(Rt) over time, on normal(0, 1) scale
  vector[n_unobs - 1] dlogrt_unscaled;
  
  // Dispersion parameter for observations
  real<lower=0> obs_dispersion;
}

transformed parameters {
  vector[n_unobs] rt = exp(
    cumulative_sum(append_row(
      log(rt_init),
      sd_dlogrt * dlogrt_unscaled
    ))
  );
  
  vector[n_unobs] x_unobs = append_row(
    x_unobs_init,
    rep_vector(0.0, n_unobs - 1)
  );
  
  vector[n_obs] lambda_obs = rep_vector(1e-10, n_obs);
  
  // Compute x_unobs by directly applying Rt
  for(t_primary in t_unobs_min:t_max) {
    for(i_inf in 1:n_dt_inf) {
      int dt_inf = dt_inf_min + i_inf - 1;
      int t_secondary = t_primary + dt_inf;
      int i_secondary = t_secondary - t_unobs_min + 1;
      if(i_secondary >= 1 && i_secondary <= n_unobs) {
        x_unobs[i_secondary] += x_unobs[t_primary] * rt[t_primary] * pmf_inf[i_inf];
      }
    }
  }
  
  // Compute lambda_obs by directly looping over all unobserved states
  // and all delays
  for(i_unobs in 1:n_unobs) {
    int t_unobs = t_unobs_min + i_unobs - 1;
    for(i_delay in 1:n_delay) {
      int delay = delay_min + i_delay - 1;
      int t_obs = t_unobs + delay;
      int i_obs = t_obs - t_obs_min + 1;
      if(i_obs >= 1 && i_obs <= n_obs) {
        lambda_obs[i_obs] += pmf_delay[i_delay] * x_unobs[i_unobs];
      }
    }
  }
}

model {
  x_unobs_init ~ normal(0, scale_x_unobs_init);
  rt_init ~ normal(0, scale_rt_init);
  dlogrt_unscaled ~ normal(0, 1);
  sd_dlogrt ~ normal(0, scale_sd_dlogrt);
  obs_dispersion ~ normal(0, scale_obs_dispersion);
  
  for(i in 1:n_obs) {
    // Gamma approximation to beta-binomial (Li/Dushoff/Bolker 2018)
    real rate = beta_binomial_gamma_rate(lambda_obs[i], p_obs, obs_dispersion);
    real shape = beta_binomial_gamma_shape(lambda_obs[i], p_obs, rate);
    y_obs[i] ~ gamma(shape, rate);
  }
}
