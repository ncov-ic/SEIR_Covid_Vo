################################################################################
# Code to fit SEIR-type compartmental model to two observations of Covid-19    #
# symptomatic and asymptomatic prevalence in Vo' Euganeo, Italy in February    #
# and March 2020                                                               #
#                                                                              #
# Cite as:                                                                     #
# Lavezzo et. al., 2020, Suppression of COVID-19 outbreak in the municipality  #
# of Vo, Italy, medRxiv, doi: 10.1101/2020.04.17.20053157                                                                #
#                                                                              #
# Input data:                                                                  #
# - observed number of negative, symptomatic and asymptomatic individuals      #
#   after two screening runs of large parts of the population of Vo' Euganeo   #
#                                                                              #
# Descriptions of classes/compartments:                                        #
# S    : susceptible individuals                                               #
# E    : individuals that are incubating the virus: they are infected,         #
#        not infectious and tests can not yet detect their infection           #
# TPp  : individuals with a higher viral load that can be detected with test   #
# I_S  : infectious individuals that show symptoms                             #
# I_A  : infectious individuals that do not show symptoms                      #
# TP_S : symptomatic individuals that are no longer infectious but have a      #
#        detectable viral load                                                 #
# TP_A : asymptomatic individuals that are no longer infectious but have a     #
#        detectable viral load                                                 #
# TN   : individuals that test negative                                        #
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
# q     : relative infectiouness of class TPp w.r.t classes I_A and I_S        #
# 1 / delta : incubation period                                                #
# 1 / gamma : infectious period                                                #
# 1 / sigma : viral load detection period                                      #
#                                                                              #
################################################################################

# Preamble --------------------------------------------------------------------#

# Set working directory if not using R-project
# setwd("yourpath")
cat("Please set your working directory or create an R-project before running code!")

library(odin)
library(reshape2)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(gridExtra)

source("R/model.R")
source("R/clean.R")
source("R/figures.R")
source("R/tables.R")

# Output folders
dir_output  <- "output"
dir_figures <- "figures"
dir.create(dir_output,  showWarnings = FALSE)
dir.create(dir_figures, showWarnings = FALSE)


# Define code parameters ------------------------------------------------------#

cat("Set the number of iterations.
As a default, each parameter combo will run three chains in sequence, each of
which will run for 200,000 iterations. This should take about an hour to run.
N.B. code can easily be parallelised using mclapply instead of lapply below")
# MCMC parameters
mcmc_iterations <- 200000 # number MCMC iterations
sample_spacing  <- 100    # every how many iterations do we save the MCMC output
id_chain    <- 1:3 # chains to run

# Parameters to clean results
cat("Make sure nr_sample < mcmc_iterations/sample_spacing - nr_burnin")
nr_burnin <- 200 # number stored parameter combinations we condiser to be burnin
nr_sample <- 100 # number stored parameter combinations we sample for plotting

# Observed data at first and second sampling
# weighted average time of first sampling = 24 Feb 2020
# weighted average time of second sampling = 07 Mar 2020
data <- data.frame(Tested         = c(2812, 2343),
                   Asymptomatic   = c(  29,   13),
                   Symptomatic    = c(  34,   15),
                   Presymptomatic = c(  10,    1))


# Define model parameters -----------------------------------------------------#

# Parameters for which we test different values
cat("Decide what parameter combinations you want to loop through.")
all_looped_parameters <- data.frame(expand.grid(
  R0_1  = c(2.1, 2.4, 2.7, 3),
  sigma = 1 / c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20)
))

# Parameters whose value is fixed in the model
fixed_parameters <- data.frame(
  tSeed   = 0,
  time1   = 22,
  time2   = 32,
  tQ      = 20,
  N       = 3275,
  generation_time = 7
)

# Parameters to fit (MCMC options)
limits <- data.frame(seed      = c(1,    100),
                     p         = c(0.01,   1),
                     w         = c(0.001,  1),
                     inv_nu    = c(0.01, fixed_parameters$generation_time),
                     inv_delta = c(0.01, fixed_parameters$generation_time),
                     q         = c(0.01,   1))
random_walk_rate <- c(seed      = 0.05,
                      p         = 0.01,
                      w         = 0.05,
                      inv_nu    = 0.01,
                      inv_delta = 0.01,
                      q         = 0.05)


# Run models (loop) -----------------------------------------------------------#
cat("N.B. code can easily be parallelised using mclapply instead of lapply")
side_effects <- lapply(seq_len(nrow(all_looped_parameters)),
                       wrapper_model, # main model function
                       data                  = data,
                       fixed_parameters      = fixed_parameters,
                       all_looped_parameters = all_looped_parameters,
                       limits                = limits,
                       random_walk_rate      = random_walk_rate,
                       mcmc_iterations       = mcmc_iterations,
                       sample_spacing        = sample_spacing,
                       dir_output            = dir_output,
                       id_chain              = id_chain)

# Process results -------------------------------------------------------------#

# Clean results
clean_mcmc_results(dir_output, nr_burnin, nr_sample)

# Make test figures
fig_chains(dir_output, dir_figures)
fig_acceptance_rates(dir_output, dir_figures)

# Generate tables for publication
table_fitted(dir_output, dir_figures)
table_TN(dir_output, dir_figures, nr_sample)
rmarkdown::render("tables.Rmd",
                  output_dir = dir_figures,
                  params = list(dir_output = dir_output))

# Generate figures for publication
## SI figure
p <- fig_prevalence(dir_output, dir_figures, data, do.best = FALSE)
ggsave(filename = file.path(dir_figures, "FigSX.tiff"),
       plot = p, device = "tiff",
       width = 183, height = 90,
       units = "mm", dpi = 300, limitsize = TRUE)

## main figure 3
pA <- fig_prevalence(dir_output, dir_figures, data, do.best = TRUE)
pB <- fig_incidence(dir_output)
pC <- fig_final_size(dir_output)
p <- arrangeGrob(pA, pB, pC,
                 layout_matrix = matrix(c(1,1,2,3), byrow = TRUE, ncol = 2))
ggsave(filename = file.path(dir_figures, "Fig3.tiff"),
       plot = p, device = "tiff",
       width = 89, height = 80,
       units = "mm", dpi = 300, limitsize = TRUE)

## main figure 1c
pC <- fig_timeline()
ggsave(filename = file.path(dir_figures, "Fig1c.tiff"),
       plot = pC, device = "tiff",
       width = 89, height = 50,
       units = "mm", dpi = 300, limitsize = TRUE)
