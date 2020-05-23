################################################################################
# Code to fit SEIR-type compartmental model to two observations of Covid-19    #
# symptomatic and asymptomatic prevalence in Vo' Euganeo, Italy in February    #
# and March 2020                                                               #
#                                                                              #
# Cite as:                                                                     #
# Lavezzo et. al., 2020, Suppression of COVID-19 outbreak in the municipality  #
# of Vo, Italy, medRxiv, doi: 10.1101/2020.04.17.20053157                      #
#                                                                              #
# Description:                                                                 #
# See README.md                                                                #
################################################################################

# Preamble --------------------------------------------------------------------#

# Set working directory if not using R-project
# setwd("yourpath")
cat("Please set your working directory to this folder or create an R-project\
 before running code!")

# load packages
library(odin)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggrepel)
library(gridExtra)

# source scripts
source("R/model.R")
source("R/clean.R")
source("R/figures.R")
source("R/tables.R")

# create output folders
dir_output  <- "output"
dir_clean   <- "clean"
dir_figures <- "figures"
dir.create(dir_output,  showWarnings = FALSE)
dir.create(dir_clean,   showWarnings = FALSE)
dir.create(dir_figures, showWarnings = FALSE)

# Define analysis parameters --------------------------------------------------#

# MCMC parameters
cat("Set the number of iterations.
As a default, each parameter combo will run three chains in sequence, each of
which will run for 200,000 iterations. This should take about an hour to run.
N.B. code can easily be parallelised using mclapply instead of lapply below")
mcmc_iterations <- 200000 # number MCMC iterations
sample_spacing  <- 100    # every how many iterations do we save the MCMC output
id_chain <- 1:3 # how many chains to run for each job - will be run in sequence

# Parameters to clean results
cat("Make sure nr_sample < mcmc_iterations/sample_spacing - nr_burnin")
nr_burnin <- 200 # number stored parameter combinations we condiser to be burnin
nr_sample <- 100 # number stored parameter combinations we sample for plotting


# Observed data ---------------------------------------------------------------#

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
looped_parameters <- data.frame(expand.grid(
  R0_1  = seq(2.1, 2.7, by = 0.3),
  sigma = 1 / seq(2, 12, by = 2)
))

# Parameters whose value is fixed in the model
fixed_parameters <- data.frame(
  tSeed   = 0,
  time1   = 22,
  time2   = 32,
  tQ      = 20,
  N       = 3275,
  q_TPp   = 1,
  q_A     = 1,
  generation_time = 7
)

# Parameters to fit
lims <- data.frame(seed      = c(1, 3), # initialisation (uniform sample)
                   p         = c(0.35, 0.45),
                   w         = c(0.05, 0.15),
                   inv_nu    = c(1.5,  2.5),
                   inv_delta = c(1,    2))
limits <- data.frame(seed      = c(1,    100), # MCMC interval boundaries
                     p         = c(0.01,   1),
                     w         = c(0.001,  1),
                     inv_nu    = c(0.01, fixed_parameters$generation_time),
                     inv_delta = c(0.01, fixed_parameters$generation_time))
random_walk_rate <- data.frame(seed      = 0.05,
                               p         = 0.01,
                               w         = 0.05,
                               inv_nu    = 0.01,
                               inv_delta = 0.01)


# Run model -------------------------------------------------------------------#

cat("N.B. code can be parallelised using parallel::mclapply instead of lapply")
side_effects <- lapply(seq_len(nrow(looped_parameters)),
                       wrapper_model, # main model function
                       data              = data,
                       fixed_parameters  = fixed_parameters,
                       looped_parameters = looped_parameters,
                       lims              = lims,
                       limits            = limits,
                       random_walk_rate  = random_walk_rate,
                       mcmc_iterations   = mcmc_iterations,
                       sample_spacing    = sample_spacing,
                       dir_output        = dir_output,
                       id_chain          = id_chain)


# Clean results ---------------------------------------------------------------#

clean_posteriors(dir_output, dir_clean, nr_burnin, nr_sample)


# Create figures and tables ---------------------------------------------------#

# Make test figures
fig_chains(dir_clean, dir_figures)
fig_acceptance_rates(dir_clean, dir_figures)

# Generate tables for publication
table_fitted(dir_clean, dir_figures, data)
table_TN(dir_clean, dir_figures, nr_sample)
rmarkdown::render("tables.Rmd",
                  output_dir = dir_figures,
                  params = list(dir_clean = dir_clean))

# Generate figures for publication

## SI figure
p <- fig_prevalence(dir_clean, dir_figures, data, do = "all")
ggsave(filename = file.path(dir_figures, paste0("FigSX.tiff")),
       plot = p, device = "tiff",
       width = 183, height = 90,
       units = "mm", dpi = 300,
       type = "cairo", compression = "lzw")

## Main figure 3
pA <- fig_prevalence(dir_clean, dir_figures, data, do = "paper_estimate")
pB <- fig_incidence(dir_clean, do = "paper_estimate")
pC <- fig_final_size(dir_clean, do = "paper_estimate")
p <- arrangeGrob(pA, pB, pC,
                 layout_matrix = matrix(c(1,1,2,3), byrow = TRUE, ncol = 2))
ggsave(filename = file.path(dir_figures, "Fig3paper_estimate.tiff"),
       plot = p, device = "tiff",
       width = 89, height = 80,
       units = "mm", dpi = 300,
       type = "cairo", compression = "lzw")

## Main figure 1c
pC <- fig_timeline()
ggsave(filename = file.path(dir_figures, "Fig1c.tiff"),
       plot = pC, device = "tiff",
       width = 89, height = 50,
       units = "mm", dpi = 300, limitsize = TRUE)
