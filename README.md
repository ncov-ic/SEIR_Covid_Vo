# Cite as:
Lavezzo et. al., 2020, Suppression of COVID-19 outbreak in the municipality
of Vo, Italy, medRxiv, doi: 10.1101/2020.04.17.20053157

# SEIR model
Calibrate a SEIR-like compartmental model to prevalence data from Vo', Italy on the number of symptomatic and asymptomatic infections during two screenings in February and March 2020
Main script: `SEIR_Covid_Vo.R`

## Inpu data
Observed number of negative, pre-symptomatic, symptomatic and asymptomatic individuals after two screening runs of large parts of the population of the municipality of Vo, Italy in February and March 2020

## Descriptions of classes/compartments
`S`    : susceptible individuals
`E`    : individuals that are incubating the virus: they are infected, not infectious and tests can not yet detect their infection
`TPp`  : individuals with a higher viral load that can be detected with test
`I_S`  : infectious individuals that show symptoms
`I_A`  : infectious individuals that do not show symptoms
`TP_S` : symptomatic individuals that are no longer infectious but have a detectable viral load
`TP_A` : asymptomatic individuals that are no longer infectious but have a detectable viral load
`TN`   : individuals that test negative

## Description of parameters:
`tSeed` : time first case has been infected
`time1` : time of first sampling
`time2` : time of second sampling
`tQ`    : time quarantine started
`N`     : population resident in Vo' Euganeo
`R0_1`  : basic reproduction number before implementation of quarantine
`1 - w` : proportional reduction of the reproduction number due to the implementation of quarantine
`seed`  : number of infected individuals at time tSeed that started the epidemic in Vo' Euganeo
`p`     : probability of being asymptomatic upon the onset of infectiousness
`q`     : relative infectiouness of class TPp w.r.t classes I_A and I_S
`1 / delta` : incubation period
`1 / gamma` : infectious period
`1 / sigma` : viral load detection period
