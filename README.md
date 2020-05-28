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

Main script: `SEIR_Covid_Vo.R`

### Input data
Observed number of negative, pre-symptomatic, symptomatic and asymptomatic individuals after two screening runs of large parts of the population of the municipality of Vo, Italy in February and March 2020

### Descriptions of classes/compartments
- `S`    : susceptible individuals
- `E`    : individuals that are incubating the virus: they are infected, not infectious and tests can not yet detect their infection
- `TPp`  : individuals with a higher viral load that can be detected with test
- `I_S`  : infectious individuals that show symptoms
- `I_A`  : infectious individuals that do not show symptoms
- `TP_S` : symptomatic individuals that are no longer infectious but have a detectable viral load
- `TP_A` : asymptomatic individuals that are no longer infectious but have a detectable viral load
- `TN`   : individuals that test negative

### Description of parameters:
- `tSeed` : time first case has been infected
- `time1` : time of first sampling
- `time2` : time of second sampling
- `tQ`    : time quarantine started
- `N`     : population resident in Vo' Euganeo
- `R0_1`  : basic reproduction number before implementation of quarantine
- `1 - w` : proportional reduction of the reproduction number due to the implementation of quarantine
- `seed`  : number of infected individuals at time tSeed that started the epidemic in Vo' Euganeo
- `p`     : probability of being asymptomatic upon the onset of infectiousness
- `q`     : relative infectiouness of class TPp w.r.t classes I_A and I_S
- `1 / nu`    : average time from infection to virus detectability
- `1 / delta` : average time from virus detectability to symptoms onset
- `1 / gamma` : average duration of symptoms
- `1 / sigma` : average time from virus detectability to recovery

## Serial interval, effective reproduction number and clusters
Estimating the serial interval for the whole study period and before and after lockdown, along with the 95% confidence interval using bootstrapping.
Estimating the effective reproduction number before and after lockdown and the 95% confidence intervals. 
Code may take a long time to run -- reduce the number of iterations.
Also contains code for plotting clusters from linelist and contact data. 

### Data cleaning code
final.data.xlsx - dataset

italy_cleaning_data.Rmd - R markdown file that converts the xlsx into the datasets needed for the analysis script. the outputs are the data files listed below

### Data files
all_contacts.rds - list of all contacts between individuals (direct, indirect, imputed household)

all_contacts.tsv - same as above but tsv format

direct_contacts.rds - list of all contacts between individuals (direct and imputed household contacts only)

direct_contacts.tsv - same as above but tsv format

linelist.rds - list of all individuals who tested positive for SARS-CoV-2, their ids, household ids, dates of first positive and negative test, last positive and negative test, and date of onset (where available)

linelist.tsv - same as above but tsv format

### Source scripts for cluster analysis
checking_functions.R - performs necessary checks needed in the other two scripts

plot_clusters.R - plots clusters of transmission as seem in extended data figure 4b

plot_tree.R - plots transmission trees - we didn't use this hear due to missing onset dates but have included here 

### Analysis code
italy_analysis.Rmd - code that executes the algorithms describes in the supplementary information for estimating the serial interval, effective reproduction number and their respective confidence intervals.

Code may run slowly due to the number of iterations (N = 10,000). 
