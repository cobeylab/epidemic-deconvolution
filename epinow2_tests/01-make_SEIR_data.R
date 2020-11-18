## Generate synthetic data
## Later, load outputs using get_sim_df()
## Plot the outptus in figs/SEIR_sim.png


source('util.R')
source('simulation.R')

## Set simulation parameters
parlist <- {
  list(
    N = 2e6, #total population size
    E_init = 0, # initial in E class
    I_init = 60, # initial in I class
    t_E = 4, # mean time in E (latent period)
    t_I = 4, # mean time in I (duration of infectiousness)
    n_t = 300, # total timesteps
    R0_vals = c(2, 0.7, 1.1), ## Steady-state R0 values held at different points in the simulation
    change_starts = c(50, 90), ## Times at which R0 starts to decrease from its previous steady state value
    change_ends = c(57, 97), ## Times at which R0 reaches its new steady state
    model_types = c('seir'), # Can also choose sir
    methods = c('ode', 'stochastic') # could also choose ode
  )
}

## Derive the mean and variance of the serial interval from the input parameters
parlist$true_mean_GI = (parlist$t_E+parlist$t_I)
parlist$true_var_GI = 2*(parlist$true_mean_GI/2)^2



## Simulate SEIR data using a stochastic (ode) model
## Results are saved to subdirectory rds/
arnaught <- with(parlist, specify_arnaught(R0_vals, change_starts, change_ends, n_t))
sim_wrapper(arnaught, parlist)
write_rds(parlist, path = 'true_pars.rds')

## Visualize synthetic data
get_sim_df() %>%
  filter(time < 300) %>%
  ggplot() +
  geom_line(aes(x = time, y = infections))+
  geom_vline(aes(xintercept = parlist$change_starts[1]), lty = 2)+ ## Dahsed line where Rt starts to decrease
  geom_vline(aes(xintercept = parlist$change_starts[2]), lty = 2)+
  ggtitle('Epidemic curve') -> inc

get_sim_df() %>% 
  filter(time < 300) %>%
  ggplot()+
  geom_line(aes(x = time, y = true_rt)) +
  geom_hline(aes(yintercept = 1), lty = 2)+
  ylab(expression(paste(R[t])))+
  ggtitle(expression(paste('Underlying ', R[t], ' values'))) -> R0
cowplot::plot_grid(R0, inc,  nrow = 2)
dir_check('figs')
ggsave('figs/SEIR_sim.png', width = 5, height = 5, units = 'in', dpi = 300)