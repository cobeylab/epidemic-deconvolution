specify_arnaught <- function(max_R0, min_R0, end_max_time, start_min_time, n_t, CONTINUOUS = FALSE){
  arnaught <- if(max_R0==min_R0){
    #constant arnaught
    max_R0
  } else if(end_max_time==start_min_time){
    #step function arnaught 
    swap_time <- end_max_time + 1
    c(rep(max_R0, swap_time), rep(min_R0, n_t-swap_time)) 
    
  } else if (CONTINUOUS == FALSE){ 
    #linearly decreasing intervention arnaught
    c(rep(max_R0, end_max_time), seq(max_R0, min_R0, length.out=(start_min_time-end_max_time+1)), rep(min_R0, n_t - start_min_time))
    
  } else { #"continuous" arnaught decrease (arctan)
    mid_time <- floor((end_max_time-start_min_time)/2)+end_max_time
    transform_fn <- function(point){
      (1/pi)*(max_R0-min_R0)*(-atan(point-mid_time)+pi/2) + min_R0
    }
    sapply(seq(0, n_t, 1), transform_fn)
  }
}

simulate_sir_example <- function(
  arnaught, t_I, N, I_init, n_t, n_steps_per_t = 10,
  method = 'stochastic'
) {
  simulate_seir(
    arnaught = arnaught,
    t_E = 0,
    t_I = t_I,
    N = N,
    S_init = N - I_init,
    E_init = 0,
    I_init = I_init,
    n_t = n_t, n_steps_per_t = n_steps_per_t,
    method = method
  )
}

simulate_seir_example <- function(
  arnaught, t_E, t_I, N, E_init, I_init, n_t, n_steps_per_t = 10,
  method = 'stochastic'
) {
  simulate_seir(
    arnaught = arnaught,
    t_E = t_E,
    t_I = t_I,
    N = N,
    S_init = N - E_init - I_init,
    E_init = E_init,
    I_init = I_init,
    n_t = n_t, n_steps_per_t = n_steps_per_t,
    method = method
  )
}



load_sims_for_one_R0 <-  function(arnaught, model_type = 'seir', method = 'stochastic'){
  data.frame(fns = list.files(path = sprintf('R0-%.1f', arnaught))) %>%
    filter(grepl(fns, pattern = method)) %>%
    pull(fns) %>%
    lapply(FUN = function(xx){
      readRDS(paste0(sprintf('R0-%.1f/', arnaught), xx)) -> tmp
      tmp$sim_df %>% 
        mutate(int_time = as.character(tmp$intervention_time), 
               dec_dur = as.character(tmp$decrease_duration))
    }) %>%
    bind_rows()
}
