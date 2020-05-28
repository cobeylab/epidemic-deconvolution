#' Richard-Lucy-type deconvolution, as per Goldstein et al. 2009 PNAS:
#' https://doi.org/10.1073/pnas.0902958106
#' 
#' Warning: very sensitive to initial conditions, either due to a bug or
#' due to an inherent property of the algorithm.
#' 
#' A fairly direct translation of the methods; totally unoptimized.
#' 
#' Function parameter names have the same meaning as the Stan models.
#' 
#' @param t_obs_min Initial timepoint for observations
#' @param y_obs Delayed observations
#' @param t_unobs_min Initial time of unobserved, undelayed time series
#' @param n_unobs Length of unobserved time series
#' @param n_iterations_max Maximum number of iterations.
#' @param n_iterations Number of iterations.
#' If provided, chi-squared stopping criteria will not be used.
deconvolve_rltype_goldsteinetal <- function(
  t_obs_min, y_obs,
  delay_min, pmf_delay,
  t_unobs_min, n_unobs,
  n_iterations_max = 50,
  n_iterations = NULL,
  x_unobs_init = NULL
) {
  # Check parameters (dynamic typing is a horror; nobody should use R)
  is_int_vector <- function(x) {
    is.numeric(x) && all(x == round(x))
  }
  is_int_scalar <- function(x) {
    is.numeric(x) && length(x) == 1 && x == round(x)
  }
  stopifnot(is_int_scalar(t_obs_min))
  stopifnot(is_int_vector(y_obs) && all(y_obs >= 0))
  stopifnot(is_int_scalar(delay_min))
  stopifnot(
    is.numeric(pmf_delay) && all(pmf_delay >= 0) &&
    abs(sum(pmf_delay) - 1) < 1e-10
  )
  stopifnot(is_int_scalar(t_unobs_min))
  stopifnot(is_int_scalar(n_unobs) && n_unobs >= 1)
  
  # Some indexing variables for convenience
  n_obs <- length(y_obs)
  t_obs_max <- t_obs_min + n_obs - 1
  t_obs_range <- t_obs_min:t_obs_max
  t_unobs_range <- t_unobs_min:(t_unobs_min + n_unobs - 1)
  n_delay <- length(pmf_delay)
  delay_range <- delay_min:(delay_min + n_delay - 1)
  
  # Function to index by time, returning 0 when out of bounds
  # (super-inefficient, non-vectorized)
  get_at_time <- function(v, t, t_min) {
    i <- t - t_min + 1
    if(i < 1 || i > length(v)) 0 else v[i]
  }
  
  # Compute q: the probability that a death in the observed series will
  # be observed at any time in the observations, used to adjust for
  # right-censoring. Sum rewritten in terms of possible delays.
  q <- sapply(t_unobs_range, function(t) {
    sum(sapply(delay_range, function(delay) {
      if(t + delay >= t_obs_min && t + delay <= t_obs_max) {
        get_at_time(pmf_delay, delay, delay_min)
      }
      else {
        0
      }
    }))
  })
  
  # Function to compute y_expected from x_unobs.
  # Sum rewritten in terms of possible delays.
  compute_y_expected <- function(x_unobs) {
    sapply(t_obs_range, function(t_obs) {
      sum(sapply(delay_range, function(delay) {
        get_at_time(
          pmf_delay, delay, delay_min
        ) * get_at_time(
          x_unobs, t_obs - delay, t_unobs_min
        )
      }))
    })
  }
  
  # Function to compute new x_unobs from previous, y_obs_expected
  # Equation 2 in methods from Goldstein et al.;
  # Sum rewritten in terms of possible delays.
  # Inner loop could easily be vectorized.
  update_x_unobs <- function(x_unobs, y_expected) {
    y_obs_over_expected <- y_obs / y_expected
    x_unobs / q * sapply(t_unobs_range, function(t) {
      sum(sapply(delay_range, function(delay) {
        get_at_time(
            pmf_delay, delay, delay_min
        ) * get_at_time(
          y_obs_over_expected, t + delay, t_obs_min
        )
      }))
    })
  }
  
  # Function to compute chi-squared statistic
  compute_chi_squared <- function(y_expected) {
    squared_error <- (y_expected - y_obs)^2
    1 / n_obs * sum(
      ifelse(
        squared_error == 0 & y_expected == 0, 0,
        squared_error / y_expected
      )
    )
  }
  
  # Initial guess: shift the observed time series by the mode
  # of the delay distribution, with a placeholder for zero values
  delay_at_pmf_mode <- delay_range[which.max(pmf_delay)]
  x_unobs <- if(is.null(x_unobs_init)) {
    x <- sapply(t_unobs_range, function(t) {
      get_at_time(y_obs, t + delay_at_pmf_mode, t_obs_min)
    })
    ifelse(x == 0, y_obs[length(y_obs)], x)
  } else {
    x_unobs_init
  }
  
  y_expected <- compute_y_expected(x_unobs)
  
  # Main loop
  iteration <- 0
  done <- FALSE
  while(!done) {
    iteration <- iteration + 1
    
    x_unobs <- update_x_unobs(x_unobs, y_expected)
    y_expected <- compute_y_expected(x_unobs)
    chi_squared <- compute_chi_squared(y_expected)
    
    # Termination condition:
    # chi-squared statistic or specified number of iterations
    if(is.null(n_iterations)) {
      if(chi_squared < 1) {
        done <- TRUE
      }
    }
    else {
      if(iteration == n_iterations) {
        done <- TRUE
      }
    }
    
    # Extra termination condition: maximum number of iterations
    if(iteration == n_iterations_max) {
      done <- TRUE
    }
  }
  
  # Return everything
  list(
    q = q,
    x_unobs = x_unobs,
    y_expected = y_expected,
    chi_squared = chi_squared,
    n_iterations = iteration
  )
}
