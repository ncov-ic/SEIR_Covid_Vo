Please cite as:
Lavezzo et. al., 2020, Suppression of COVID-19 outbreak in the municipality of Vo, Italy, medRxiv, doi: 10.1101/2020.04.17.20053157

# Analysis structure
Each subfolder is a separate piece of analysis. Please set your working directory to the subfolder of the chosen analysis.

## Logistic regression
Logistic regression on linelist data.
Testing association between positivity and age/sex.

##


## SEIR model
Calibrate a SEIR-like compartmental model to prevalence data from Vo', Italy on the number of symptomatic and asymptomatic infections
 during two screenings in February and March 2020.
Computations might take a long time, decrease the number of implemented models and/or of MCMC iterations.

## Serial interval, effective teproduction number and clusters
Estimating the serial interval for the whole study period and before and after lockdown, along with the 95% confidence interval using bootstrapping.
Estimating the effective reproduction number before and after lockdown and the 95% confidence intervals. 
Code may take a long time to run -- reduce the number of iterations.
Also contains code for plotting clusters from linelist and contact data. 
