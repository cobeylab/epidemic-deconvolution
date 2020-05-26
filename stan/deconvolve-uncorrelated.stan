data {
  // Smallest (most negative) delay
  int dt_min;
  
  // Length of delay distribution
  int np;
  
  // PMF of delay distribution
  simplex[np] p;
  
  // Initial time of observed time series
  int ty_min;
  
  // Length of observed, delayed time series
  int ny;
  
  // Delayed observations
  int<lower=0> y[ny];
  
  // Initial time of unobserved, undelayed time series
  int tx_min;
  
  // Length of unobserved, undelayed time series
  int nx;
}

transformed data {
  // Scale of y observations, used to scale priors for x
  real sd_scale_x = max(y);
}

parameters {
  real<lower=0> poisson_floor;
  real<lower=0> scale_x;
  vector<lower=0>[nx] lambda_x;
}

transformed parameters {
  vector[ny] lambda_y;
  
  for(iy in 1:ny) {
    int ty = ty_min + iy - 1;
    lambda_y[iy] = 0.0;
    for(ip in 1:np) {
      int dt = dt_min + ip - 1;
      int t = ty - dt;
      int x_index = t - tx_min + 1;
      if(x_index >= 1 && x_index <= nx) {
        lambda_y[iy] += p[ip] * lambda_x[x_index];
      }
    }
    lambda_y[iy] = poisson_floor + lambda_y[iy];
  }
}

model {
  poisson_floor ~ exponential(1.0);
  scale_x ~ normal(0, sd_scale_x);
  lambda_x ~ normal(0, scale_x);
  y ~ poisson(lambda_y);
}
