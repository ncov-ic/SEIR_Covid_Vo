Please cite as:
Lavezzo et. al., 2020, Suppression of COVID-19 outbreak in the municipality of Vo, Italy, medRxiv, doi: 10.1101/2020.04.17.20053157

# Repository structure
This repository contains 
- the data collected in Vo to assess COVID-19 prevalence at two subsequent surveys conducted in February and March 2020;
- the code used to fit a compartmental SEIR-like model to the observed prevalence data; 
- the code used to estimate the serial interval and reproduction number from a reconstruction of the transmission chains; 
- the code used to assess the association between COVID-19 infection and age or sex via logistic regression. 

## Data
The data folder contains the file `anonymised_data_public_final.xlsx`, which has four sheets:
- anonymised data (sheet 1) contains a line list with contact data, dates and results of swab testing, symptoms and hospitalisation data;
- legend (sheet 2) contains a description of the information contained in sheet 1;
- RT_PCR_comparison (sheet 3) contains the viral load data (both Ct values and genome equivalents) and a statistical analysis of these data; 
- RT_PCR_DATA (sheet 4) contains the viral load data in sheet 3 linked to the subjects IDs of sheet 1. 

## Scripts
The scripts folder contains the main scripts for each separate analysis
- `SEIR_Covid_Vo.R` is the main script for fitting a SEIR-like model to prevalence data
- `italy_analysis_script.R` is the main script for estimating the serial interval and the reproduction number. This script requires script `italy_cleaning_data_script.R` to be run beforehand. 
- `Logistic_regression.R` is the main script to assess the association between swab positivity and age or sex 

## R 
The R folder contains the functions used by the main scripts. 

`SEIR_Covid_Vo.R` uses the following scripts:
- model.R
- clean.R
- figures.R 
- tables.R

`italy_analysis_script.R` uses the following scripts:
- checking_functions.R
- serial_interval.functions.R 
- reproduction_number_functions.R 
- plot_clusters.R 
- plot_tree.R


##

## Notes on the fit of the compartmental model to prevalence data
The code calibrates a compartmental model of SARS-CoV-2 transmission to the prevalence data observed in Vo using the Metropolis-Hastings Markov Chain Monte Carlo (MCMC) method. The computational time depends on the number of MCMC iterations and chains. Convergence was assessed from three MCMC chains starting at different initial points: for all parameter combinations we run 200,000 MCMC iterations, thinned the chains taking every 100th iteration and removed the first 200 observations (burnin). 

### Input data
Number of subjects tested and observed number of pre-symptomatic, symptomatic and asymptomatic study participants testing positive to SARS-CoV-2 at two surveys conducted in the municipality of Vo, Italy in February and March 2020. 

### Descriptions of model compartments
- `S`    : susceptible individuals
- `E`    : individuals that are incubating the virus: infected but not yet infectious and with undetectable viral load
- `TP`  :  individuals with detectable viral load before the onset of symptoms, we assume they are infectious
- `I_S`  : infectious individuals with detectable viral load who show symptoms
- `I_A`  : infectious individuals with detectable viral load who do not show symptoms
- `TP_S` : symptomatic individuals that are no longer infectious but have a detectable viral load
- `TP_A` : asymptomatic individuals that are no longer infectious but have a detectable viral load
- `TN`   : individuals who test negative (undetectable viral load)

### Description of parameters:
- `tSeed` : time infection is seeded in Vo
- `time1` : time of first survey
- `time2` : time of second survey
- `tQ`    : time lockdown started
- `N`     : Vo resident population
- `R0_1`  : basic reproduction number before the implementation of lockdown
- `1 - w` : proportional reduction of the basic reproduction number due to the lockdown 
- `seed`  : number of infectious individuals at time tSeed, needed to trigger the epidemic
- `p`     : proportion of asymptomatic infections (not developing symptoms throughout the whole infection)
- `1 / nu`    : average time from infection to virus detectability
- `1 / delta` : average time from virus detectability to symptom onset
- `1 / gamma` : average duration of infectiousness from time of symptom onset 
- `1 / sigma` : average duration of virus detectability beyond the infectious period 

##

## Notes on the serial interval and effective reproduction number estimation
The code reconstructs transmission chains following the algorithms described in Supplementary Text S1 and S2. The serial interval and the effective reproduction number are estimated from the reconstructed transmission chains. The 95% confidence interval around the central estimates is calculated by bootstrapping. The results presented in the paper were obtained with 10,000 iterations. The computational time depends on the number of iterations used. 

### Data preparation 
`italy_cleaning_data_script.R` uses the data file to generate the following outputs: 
- `all_contacts.rds` lists all contacts between individuals (direct, indirect, imputed household)
- `all_contacts.tsv` - same as `all_contacts.rds` but in tsv format
- `direct_contacts.rds` lists all contacts between individuals (direct and imputed household contacts only)
- `direct_contacts.tsv` - same as `direct_contacts.rds` but in tsv format
- `linelist.rds` lists all individuals who tested positive for SARS-CoV-2, their ids, household ids, dates of first positive and negative test, last positive and negative test, and date of onset (where available)
- `linelist.tsv` - same as `linelist.rds` but in tsv format

These output files are used in `italy_analysis_script.R` to execute the algorithms describes in Supplementary Text S1 and S2 for estimating the central estimates of the serial interval and effective reproduction number and their 95% confidence intervals.

### Clusters plot
`plot_clusters.R` plots the observed transmission clusters as shown in Extended Data Figure 4b. 

##

## Notes on logistic regression
The code in `Logistic_regression.R` does not depend on any external function. It runs logistic regresison on the line list data to test the association between positivity and age and sex.

