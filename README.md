# SEIR_Covid_Vo
Calibrate a SEIR-like compartmental model to prevalence data from Vo', Italy on the number of symptomatic and asymptomatic infections during two screenings in February and March 2020

## Cite as:                                                                    
Lavezzo et. al., 2020, Suppression of COVID-19 outbreak in the municipality of Vo’, Italy                                                               
                                                                            
## Input data:                                                                 
observed number of negative, symptomatic and asymptomatic individuals during two screening runs of large parts of the population of Vo', Italy 
                                                                            
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
`1 / delta` : incubation period                                               
`1 / gamma` : infectious period                                               
`1 / sigma` : viral load detection period                                     
                                                                            
