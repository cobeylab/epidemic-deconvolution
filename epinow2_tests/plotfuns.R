library(ggplot2)
theme_set(theme_bw())


plot_pr_vs_post <- function(
  post_samples,
  incubation_period,
  generation_time,
  delay,
  delay_type = 'case',
  testpars,
  parlist
){
  gi_prior <- data_frame(
    xx = seq(0, generation_time$max, by = 0.01),
    prior = dlnorm(xx, generation_time$mean, generation_time$sd),
    true = dgamma(xx, get_shape(parlist$true_mean_GI, parlist$true_var_GI), get_rate(parlist$true_mean_GI, parlist$true_var_GI))
  )
  inc_prior <- data_frame(
    xx = seq(0, incubation_period$max, by = 0.01),
    prior = dlnorm(xx, incubation_period$mean, incubation_period$sd),
    true = dlnorm(xx, log(testpars$true_mean_inc), log(testpars$true_sd_inc))
  )
  del_prior <- data_frame(
    xx = seq(0, case_delay$max, by = 0.01),
    prior = dlnorm(xx, delay$mean, delay$sd)
  )
  if(delay_type == 'case'){
    del_prior$true = dlnorm(del_prior$xx, log(testpars$true_mean_case_delay), log(testpars$true_sd_case_delay))
  }else{
    del_prior$true = dlnorm(del_prior$xx, log(testpars$true_mean_death_delay), log(testpars$true_sd_death_delay))
  }
  
  ## Write a function to evaluate the posterior distribution at each timestep
  get_post <- function(ests, truepars, vname = 'gen_interval'){
    post_samples %>%
      as.tbl() %>%
      filter(grepl(pattern = 'gt', variable)) %>%
      group_by(variable, sample) %>%
      summarise(value = unique(value)) %>%
      pivot_wider(names_from = variable, values_from = value) %>%
      apply(MARGIN = 1, FUN = function(ss) dlnorm(seq(0, truepars$max, by = 0.01), ss[2], ss[3])) %>%
      apply(MARGIN = 1, function(rr){
        c(lower = quantile(rr, .025), 
          median = quantile(rr, .5),
          upper = quantile(rr, .975),
          mean = mean(rr))
      }) %>%
      t() %>%
      as.data.frame() %>%
      mutate(xx = seq(0, truepars$max, by = 0.01),
             variable = vname)
  }
  
  gi_post <- get_post(ests, truepars = generation_time, vname = 'gi_post')
  inc_post <- get_post(ests, truepars = incubation_period, vname = 'inc_post')
  delay_post <- get_post(ests, truepars = delay, vname = 'delay_post')
  
  
  gg <- merge(gi_prior, gi_post, by = 'xx') %>%
    pivot_longer(cols = c('prior', 'true', 'mean')) %>%
    mutate(name = ifelse(name == 'mean', 'posterior', name)) %>%
    ggplot()+
    geom_line(aes(x = xx, y = value, color = name, lty = name), alpha = .5)+
    geom_ribbon(aes(x = xx, ymin = `lower.2.5%`, ymax = `upper.97.5%`), fill = 'yellow', alpha = .5)+
    scale_color_manual("", values = c('goldenrod', 'red', 'black'))+
    scale_linetype_manual("", values = c(1, 2, 1))+
    ylab('dens')+xlab('value')+
    ggtitle('Generation interval')+
    theme(legend.position = c(.8, .8))
  
  ii <- merge(inc_prior, inc_post, by = 'xx') %>%
    pivot_longer(cols = c('prior', 'true', 'mean')) %>%
    mutate(name = ifelse(name == 'mean', 'posterior', name)) %>%
    ggplot()+
    geom_line(aes(x = xx, y = value, color = name, lty = name), alpha = .5)+
    geom_ribbon(aes(x = xx, ymin = `lower.2.5%`, ymax = `upper.97.5%`), fill = 'yellow', alpha = .5)+
    scale_color_manual("", values = c('goldenrod', 'red', 'black'))+
    scale_linetype_manual("", values = c(1, 2, 1))+
    ylab('dens')+xlab('value')+
    ggtitle('Incubation period')+
    theme(legend.position = c(.8, .8))
  
  dd <- merge(del_prior, delay_post, by = 'xx') %>%
    pivot_longer(cols = c('prior', 'true', 'mean')) %>%
    mutate(name = ifelse(name == 'mean', 'posterior', name)) %>%
    ggplot()+
    geom_line(aes(x = xx, y = value, color = name, lty = name), alpha = .5)+
    geom_ribbon(aes(x = xx, ymin = `lower.2.5%`, ymax = `upper.97.5%`), fill = 'yellow', alpha = .5)+
    scale_color_manual("", values = c('goldenrod', 'red', 'black'))+
    scale_linetype_manual("", values = c(1, 2, 1))+
    ylab('dens')+xlab('value')+
    ggtitle('Reporting delay')+
    theme(legend.position = c(.8, .8))
  
  cowplot::plot_grid(gg, ii, dd, ncol = 3)
  ggsave2(filename = sprintf('figs/priors_vs_post-%s-%s.png', path, delay_type), width = 8, height = 4, units = 'in', dpi = 300)
}


plot_rt <- function(post_summary){
  merge(
    ## True rt values from simulation
    get_sim_df() %>%
      mutate(date = as.Date(min(post_summary$date, na.rm = T))+time) %>%
      select(date, true_rt),
    ## Posterior rt estimates
    post_summary %>% as.tbl() %>% filter(variable == 'R'),
    by = 'date'
  ) %>%
    select(date, true_rt, lower, upper, mean) %>%
    rename(posterior = mean) %>%
    pivot_longer(c(posterior, true_rt)) %>%
    ggplot(aes(x = date))+
    geom_line(aes(y = value, color = name)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = 'yellow', alpha = .7) +
    scale_color_manual("", values = c('goldenrod', 'black'))
  ggsave2(filename = sprintf('figs/rt-%s-%s.png', path, delay_type), width = 5, height = 4, units = 'in', dpi = 300)
}




plot_epi_curve <- function(post_summary, obs_cases){
  merge(
    ## True rt values from simulation
    get_sim_df() %>%
      mutate(date = as.Date(min(post_summary$date, na.rm = T))+time) %>%
      select(date, incidence),
    ## Posterior rt estimates
    post_summary %>% as.tbl() %>% filter(variable == 'infections'),
    by = 'date'
  ) %>% 
  merge(obs_cases, by = 'date')%>%
  select(date, incidence, lower, upper, mean, confirm) %>%
  setnames(c('date', 'true infections', 'lower', 'upper', 'estimated infections', 'true observations')) %>%
  pivot_longer(c(`true infections`, `estimated infections`, `true observations`)) %>%
    ggplot(aes(x = date))+
    geom_line(aes(y = value, color = name)) + 
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = 'yellow', alpha = .5)+
    scale_color_manual("", values = c('goldenrod', 'black', 'blue3'))
  ggsave2(filename = sprintf('figs/infections-%s-%s.png', path, delay_type), width = 5, height = 4, units = 'in', dpi = 300)
}




make_plots <- function(path,
                       delay_type = 'cases'){
  
  
  gen_int <- readRDS(sprintf('%s/gen_interval.rds', path))
  inc_pd <- readRDS(sprintf('%s/incubation_pd.rds', path, delay_type))
  case_delay <- readRDS(sprintf('%s/case_delay.rds', path, delay_type))
  death_delay <- readRDS(sprintf('%s/death_delay.rds', path, delay_type))
  testpars <- readRDS(sprintf('%s/testpars.rds', path, delay_type))
  parlist <- readRDS('true_pars.rds')
  post_samples <- readRDS(sprintf('%s/%s/latest/estimate_samples.rds', path, delay_type))
  post_summary <- readRDS(sprintf('%s/%s/latest/summarised_estimates.rds', path, delay_type))
  case_summary <- readRDS(sprintf('%s/%s/latest/summarised_estimated_reported_cases.rds', path, delay_type))
  obs_cases <- readRDS(sprintf('%s/%s/latest/reported_cases.rds', path, delay_type))
  
  ## Plot posterior vs. prior vs. true distributions
  ##    Posterior is yellow
  ##    Prior is dashed orange
  ##    Truth is black line
  plot_pr_vs_post(post_samples, inc_pd, gen_int, case_delay, delay_type = 'case', testpars, parlist)
  
  ## Plot posterior Rt (yellow) vs true value
  plot_rt(post_summary)
  
  ## Plot inferred incidence vs. truth
  plot_epi_curve(post_summary, obs_cases)
}

