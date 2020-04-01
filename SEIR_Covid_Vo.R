################################################################################
# Code to fit SEIR-type compartmental model to two observations of Covid-19    #
# symptomatic and asymptomatic prevalence in Vo' Euganeo, Italy in February    #
# and March 2020                                                               #
#                                                                              #
# Cite as:                                                                     #
# Lavezzo et. al., 2020, Suppression of COVID-19 outbreak in the municipality  #
# of Vo’, Italy                                                                #
#                                                                              #
# Input data:                                                                  #
# - observed number of negative, symptomatic and asymptomatic individuals      #
#   after two screening runs of large parts of the population of Vo' Euganeo   #
#                                                                              #
# Description of parameters:                                                   #
# tSeed : time first case has been infected                                    #
# time1 : time of first sampling                                               #
# time2 : time of second sampling                                              #
# tQ    : time quarantine started                                              #
# N     : population resident in Vo' Euganeo                                   #
# R0_1  : basic reproduction number before implementation of quarantine        #
# 1 - w : proportional reduction of the reproduction number due to the         #
#         implementation of quarantine                                         #
# seed  : number of infected individuals at time tSeed that started the        #
#         epidemic in Vo' Euganeo                                              #
# p     : probability of being asymptomatic upon the onset of infectiousness   #
# 1 / delta : incubation period                                                #
# 1 / gamma : infectious period                                                #
# 1 / sigma : viral load detection period                                      #
#                                                                              #
################################################################################


# Preamble --------------------------------------------------------------------#

# Set working directory if not using R-project
# setwd("yourpath")

library(deSolve)
library(reshape2)
library(dplyr)
library(ggplot2)

source("R/model.R")
source("R/clean.R")
source("R/figures.R")

# Output folders
dir_output  <- "output"
dir_figures <- "figures"
dir.create(dir_output,  showWarnings = FALSE)
dir.create(dir_figures, showWarnings = FALSE)


# Define parameters -----------------------------------------------------------#

# MCMC parameters
mcmc_iterations <- 1000#100000 # number MCMC iterations
sample_spacing  <- 10#100    # every how many iterations do we save the MCMC output

# Parameters to clean results
nr_burnin <- 1#100 # number stored parameter combinations do we condiser to be burnin
nr_sample <- 9#100 # number stored parameter combinations do we sample for plotting

# Observed data at first and second sampling
# weighted average time of first sampling = 24 Feb 2020
# weighted average time of second sampling = 07 Mar 2020
data <- data.frame(Tested       = c(2812, 2343),
                   Asymptomatic = c(  30,   13),
                   Symptomatic  = c(  43,   16))

# Parameters for which we test different values
all_looped_parameters <- data.frame(expand.grid(
  R0_1  = c(2.1, 2.4, 2.7),
  sigma = 1 / c(2, 4, 6, 8, 10, 12)
))

# Parameters whose value is fixed in the model
fixed_parameters <- data.frame(
  tSeed   = 0,
  time1   = 22,
  time2   = 32,
  tQ      = 20,
  N       = 3275,
  delta   = 1/5,
  gamma   = 1/2
)

# Parameters to fit
limits <- data.frame(seed = c(1,    100),
                     p    = c(0.01,   1),
                     w    = c(0.001,  1))
random_walk_rate <- c(seed = 0.05,
                      p    = 0.01,
                      w    = 0.05)


# Run models (loop) -----------------------------------------------------------#
side_effects <- lapply(seq_len(nrow(all_looped_parameters)),
                       wrapper_model, # main model function
                       solve_seir            = solve_seir,
                       data                  = data,
                       fixed_parameters      = fixed_parameters,
                       all_looped_parameters = all_looped_parameters,
                       limits                = limits,
                       random_walk_rate      = random_walk_rate,
                       mcmc_iterations       = mcmc_iterations,
                       sample_spacing        = sample_spacing,
                       dir_output            = dir_output)


# Process results -------------------------------------------------------------#

# Clean results
clean_mcmc_results(dir_output, nr_burnin, nr_sample)

# Make figures
fig_chains(dir_output, dir_figures)
fig_acceptance_rates(dir_output, dir_figures)
fig_SEIR(dir_output, dir_figures, data)

# Make tables
table_fitted(dir_output)
table_TN(dir_output, nr_sample)
