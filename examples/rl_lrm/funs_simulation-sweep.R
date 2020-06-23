#!/usr/bin/env Rscript
source('simulation.R')
source('util.R')

##Author: Katie Gostic

## Define functions to simulate and plot data using SEIR 

sim_sweep <- function(PARAMS){
  
  dirname <- sprintf("R0-%.1f", PARAMS$pre_intervention_R0)
  if(!dir.exists(dirname)){
    sprintf('creating new output directory, %s', dirname)
    dir.create(dirname)
  }
  setwd(dirname)
  saveRDS(PARAMS, paste0('PARAMS_df.Rds'))
  for(intervention_time_1 in PARAMS$intervention_time_1){
    for(decrease_duration in PARAMS$days_intervention_to_min){
      for(intervention_time_2 in PARAMS$intervention_time_2){
        for(increase_duration in PARAMS$days_to_Rt_rise){
          arnaught <- c(
            rep(PARAMS$pre_intervention_R0, intervention_time_1),
            seq(PARAMS$pre_intervention_R0, PARAMS$intervention_R0, length.out = decrease_duration+1),
            rep(PARAMS$intervention_R0, intervention_time_2-decrease_duration-intervention_time_1-1),
            seq(PARAMS$intervention_R0, PARAMS$partially_lifeted_R0, length.out = increase_duration+1),
            rep(PARAMS$partially_lifeted_R0, length.out = PARAMS$n_t - increase_duration -intervention_time_2)
          )
          #plot(arnaught); abline(v = intervention_time+c(0, decrease_duration), col = 'red'); 
          do_one_arnaught(arnaught, PARAMS$intervention_time_1, PARAMS$days_intervention_to_min, PARAMS)
        }
      }
    }
  }
  setwd('..')
}

do_one_arnaught <- function(arnaught, intervention_time, decrease_duration, PARAMS) {
  for(method in PARAMS$methods) {
    for(model_type in PARAMS$model_types) {
      sim_df <- if(model_type == 'sir') {
        simulate_sir_example(
          arnaught = arnaught,
          t_I = PARAMS$t_I,
          N = PARAMS$N, I_init = PARAMS$I_init,
          n_t = PARAMS$n_t,
          method = method
        ) %>%
          mutate(
            incidence = round(dS),
            obs_cases = NA
          )
      } else {
        simulate_seir_example(
          arnaught = arnaught,
          t_E = PARAMS$t_E, t_I = PARAMS$t_I,
          N = PARAMS$N, E_init = PARAMS$E_init, I_init = PARAMS$I_init,
          n_t = PARAMS$n_t,
          method = method
        ) %>%
          mutate(
            incidence = round(dS),
            obs_cases = round(dEI)
          )
      }
      
      saveRDS(
        list(
          sim_df = sim_df %>% 
            mutate(true_r0 = arnaught,
                   true_rt = arnaught*S/(S+E+I+R)),
          arnaught = arnaught,
          intervention_time = intervention_time,
          decrease_duration = decrease_duration,
          params = PARAMS,
          method = method,
          model_type = model_type
        ),
        sprintf('%s_%s_dec%.0f-%.0f_sim.rds', model_type, method, intervention_time, decrease_duration)
      )
    }
  }
}
  
 
#sim_sweep(PARAMS)
#END UP WITH one folder for each R0, which contains: SEIR stoc {decreases}; SEIR ode {decreases}, parameters



## Check outputs
testplots <- function(PARAMS) {
  for(arnaught in PARAMS$pre_intervention_R0) {
    for(model_type in c('seir')) {
      for(method in PARAMS$methods) {
        # Load results from all interventions applied to a given R0, and bind into a single data frame
       sim_results <- data.frame(fns = list.files(path = sprintf('R0-%.1f', arnaught)), stringsAsFactors = FALSE) %>%
         filter(grepl(fns, pattern = method)) %>%
         filter(grepl(fns, pattern = '_sim')) %>%
         pull(fns) %>%
         lapply(FUN = function(xx){
           readRDS(paste0(sprintf('R0-%.1f/', arnaught), xx)) -> tmp
           tmp$sim_df %>% 
             mutate(int_time = as.character(tmp$intervention_time), 
                    dec_dur = as.character(tmp$decrease_duration))
         }) %>%
         bind_rows()
        
            cat(sprintf('Plotting %s, %s', model_type, method))
            # sim_result <- readRDS(sprintf('R0-%.1f/%s_%s_dec%.0f-%.0f_sim.rds', 
            #                               arnaught, model_type, method, intervention_time, decrease_duration))
            
            ## Plot prevalence in each compartment
            sim_results %>%
              select(time:R, int_time, dec_dur) %>%
              pivot_longer(S:R, names_to = 'Compartment', values_to = 'Prevalence') %>%
              ggplot() +
              geom_line(aes(x = time, y = Prevalence, color = int_time, lty = dec_dur)) +
              facet_wrap(.~Compartment, scales = 'free_y') +
              ggtitle(sprintf('%s - %s - R0=%.1f', model_type, method, arnaught))
            if(!dir.exists('simplots/')){ dir.create('simplots/')}
            ggsave(filename = sprintf('simplots/prevalence_R0-%.1f_%s_%s.png', arnaught, model_type, method), width = 11, height = 8.5)
            
            ## Plot incidence in each compartment
            sim_results %>%
              select(time, int_time, dec_dur, dEI, dIR, incidence) %>%
              pivot_longer(dEI:incidence, names_to = 'transition', values_to = 'incidence') %>%
              ggplot() +
              geom_line(aes(x = time, y = incidence, color = int_time, linetype = dec_dur), size = 1, alpha = .5) +
              theme_bw()+
              facet_wrap(.~transition) +
              ggtitle(sprintf('%s - %s - R0=%.1f', model_type, method, arnaught))
            ggsave(filename = sprintf('simplots/incidence_R0-%.1f_%s_%s.png', arnaught, model_type, method), width = 11, height = 8.5)
          }
        }
      }
    }

#testplots(PARAMS)
