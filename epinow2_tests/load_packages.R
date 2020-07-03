## Load packages
## Install packages
rm(list = ls())
require(rstan, quietly = TRUE)
rstan_options(auto_write = TRUE)
require(dplyr, quietly = TRUE)
require(readr, quietly = TRUE)
require(tidyr, quietly = TRUE)
require(ggplot2, quietly = TRUE)
require(cowplot, quietly = TRUE)
require(EpiNow2, quietly = TRUE)
require(NCoVUtils, quietly = TRUE)
require(deSolve, quietly = TRUE)
require(data.table, quietly = TRUE) 
require(future, quietly = TRUE)
#require(forecastHybrid, quietly = TRUE)

