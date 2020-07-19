#' Simulate a SEIR model, deterministic or stochastic.
#' 
#' written by Ed Baskerville, reviewed/edited by Timothy M Pollington
#' 18 July 2020 (v2)
#' 
#' No age structure.
#' 
#' @param n_t Number of units of time to simulate.
#' @param n_steps_per_t Number of discrete simulation steps to take per
#' unit of simulation time.
#' @param arnaught Basic reproduction number (ratio), squiggly-R-0.
#' Average number of new infections produced by an infection in a
#' susceptible population. A scalar or a vector of length `n_t + 1`,
#' which specifies R0 at the start of each timestep. R0 is linearly
#' interpolated between timesteps.
#' @param t_E Mean latent period. If set to 0, the model reduces to an
#' SIR.
#' @param t_I Mean duration of infectiousness.
#' 
#' @return A dataframe containing `time`, all compartments
#' (`S`, `E`, `I`, `R`), and transition counts (`dS`, `dEI`, `dIR`).
simulate_seir <- function(
  arnaught, t_E, t_I,
  N, S_init, E_init, I_init,
  n_t, n_steps_per_t = 1,
  method = 'stochastic'
) {
  func <- if(method == 'stochastic') simulate_seir_stochastic else simulate_seir_ode
  func(
    arnaught, t_E, t_I,
    N, S_init, E_init, I_init,
    n_t, n_steps_per_t
  )
}

simulate_seir_stochastic <- function(
  arnaught, 
  t_E, 
  t_I,
  N, 
  S_init, 
  E_init, 
  I_init,
  n_t, 
  n_steps_per_t
) {
  check_args(
    arnaught, t_E, t_I, N, S_init, E_init, I_init, n_t, n_steps_per_t
  ) 
  
  delta_t <- 1 / n_steps_per_t
  
  # Draws a binomial based on a rate
  draw <- function(n, rate) {
    stopifnot(rate >= 0)
    p <- 1 - exp(-rate * delta_t)
    rbinom(1, n, p)
  }
  
  # Function to compute beta at a particular time
  beta <- construct_beta(arnaught, t_I, n_t)
  
  # Step forward from t to t + delta_t
  step <- function(t, S_prev, E_prev, I_prev) {
    dS <- draw(S_prev, beta(t) * I_prev / N)
    dIR <- draw(I_prev, 1 / t_I)
    
    if(t_E > 0) {
      # SEIR model
      dEI <- draw(E_prev, 1 / t_E)
      list(
        S = S_prev - dS,
        E = E_prev + dS - dEI,
        I = I_prev + dEI - dIR,
        dS = dS,
        dEI = dEI,
        dIR = dIR
      )
    }
    else {
      # SIR model
      list(
        S = S_prev - dS,
        E = 0,
        I = I_prev + dS - dIR,
        dS = dS,
        dEI = 0,
        dIR = dIR
      )
    }
  }
  
  # Set up state vectors over time
  S <- numeric(n_t + 1)
  S[1] <- S_init
  E <- numeric(n_t + 1)
  I <- numeric(n_t + 1)
  if(t_E > 0) {
    # SEIR model
    E[1] <- E_init
    I[1] <- I_init
  }
  else {
    # SIR model: all initial E get dumped in I
    E[1] <- 0
    I[1] <- E_init + I_init
  }
  
  # Track transitions over time
  dS <- rep(NA, n_t + 1)
  dEI <- rep(NA, n_t + 1)
  dIR <- rep(NA, n_t + 1)
  
  # Simulate
  for(i in 1:n_t) {
    S_prev <- S[i]
    E_prev <- E[i]
    I_prev <- I[i]
    
    dS[i+1] <- 0
    dEI[i+1] <- 0
    dIR[i+1] <- 0
    for(j in 1:n_steps_per_t) {
      state_next <- step(i + delta_t * (j - 1), S_prev, E_prev, I_prev)
      S_prev <- state_next$S
      E_prev <- state_next$E
      I_prev <- state_next$I
      dS[i+1] <- dS[i+1] + state_next$dS
      dEI[i+1] <- dEI[i+1] + state_next$dEI
      dIR[i+1] <- dIR[i+1] + state_next$dIR
    }
    
    S[i+1] <- S_prev
    E[i+1] <- E_prev
    I[i+1] <- I_prev
  }
  
  # Return each compartment over time
  data.frame(
    time = 0:n_t,
    S = S,
    E = E,
    I = I,
    R = N - S - E - I,
    dS = dS,
    dEI = dEI,
    dIR = dIR
  )
}

simulate_seir_ode <- function(
  arnaught, t_E, t_I,
  N, S_init, E_init, I_init,
  n_t,
  n_steps_per_t = 1 # Ignored; included so the function signature matches stochastic version
) {
  check_args(
    arnaught, t_E, t_I, N, S_init, E_init, I_init, n_t, n_steps_per_t
  )
  
  beta <- construct_beta(arnaught, t_I, n_t)
  d_dt <- function(t, y, params) {
    dS <- y['S'] * beta(t) * y['I'] / N
    dIR <- y['I'] / t_I
    
    if(t_E > 0) {
      # SEIR model
      dEI <- y['E'] / t_E
      list(c(
        S = -dS,
        E = dS - dEI,
        I = dEI - dIR,
        R = dIR,
        cum_dS = dS,
        cum_dEI = dEI
      ), NULL)
    }
    else {
      # SIR model
      list(c(
        S = -dS,
        E = 0,
        I = dS - dIR,
        R = dIR,
        cum_dS = dS,
        cum_dEI = dS
      ), NULL)
    }
  }
  
  y_init <- c(
    S = S_init,
    E = if(t_E > 0) E_init else 0,
    I = if(t_E > 0) I_init else E_init + I_init,
    R = 0,
    cum_dS = 0,
    cum_dEI = 0
  )
  as.data.frame(deSolve::ode(y_init, 0:n_t, d_dt, NULL)) %>%
    mutate(dS = cum_dS - lag(cum_dS, 1)) %>%
    mutate(dEI = cum_dEI - lag(cum_dEI, 1)) %>%
    mutate(dIR = R - lag(R, 1))
}

#' Helper to construct a beta(t) function (constant or time-varying)
construct_beta <- function(arnaught, t_I, n_t) {
  beta_t_all <- arnaught / t_I
  if(length(arnaught) == 1) {
    function(t) beta_t_all
  } else {
    approxfun(0:n_t, beta_t_all)
  }
}

#' Helper to validate function arguments
check_args <- function(
  arnaught, t_E, t_I, N, S_init, E_init, I_init, n_t, n_steps_per_t
) {
  # Check all input parameters
  stopifnot(length(S_init) == 1 && S_init > 0 && round(S_init) == S_init)
  stopifnot(length(E_init) == 1 && E_init >= 0 && round(E_init) == E_init)
  stopifnot(length(I_init) == 1 && I_init >= 0 && round(I_init) == I_init)
  stopifnot(length(N) == 1 && N >= 1 && round(N) == N)
  stopifnot(length(n_t) == 1 && n_t >= 1 && round(n_t) == n_t)
  stopifnot(
    length(n_steps_per_t) == 1 && n_steps_per_t >= 1 &&
    round(n_steps_per_t) == n_steps_per_t
  )
  stopifnot(
    is.numeric(arnaught) && arnaught > 0 &&
      (length(arnaught) == 1 || length(arnaught) == n_t + 1)
  )
  stopifnot(
    is.numeric(t_E) && length(t_E) == 1 && t_E >= 0
  )
  stopifnot(
    is.numeric(t_I) && length(t_I) == 1 && t_I > 0
  )
}
