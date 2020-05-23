## Data cleaning code
final.data.xlsx - dataset

italy_cleaning_data.Rmd - R markdown file that converts the xlsx into the datasets needed for the analysis script. the outputs are the data files listed below

## Data files
all_contacts.rds - list of all contacts between individuals (direct, indirect, imputed household)

all_contacts.tsv - same as above but tsv format

direct_contacts.rds - list of all contacts between individuals (direct and imputed household contacts only)

direct_contacts.tsv - same as above but tsv format

linelist.rds - list of all individuals who tested positive for SARS-CoV-2, their ids, household ids, dates of first positive and negative test, last positive and negative test, and date of onset (where available)

linelist.tsv - same as above but tsv format

## Source scripts for cluster analysis
checking_functions.R - performs necessary checks needed in the other two scripts

plot_clusters.R - plots clusters of transmission as seem in extended data figure 4b

plot_tree.R - plots transmission trees - we didn't use this hear due to missing onset dates but have included here 

## Analysis code
italy_analysis.Rmd - code that executes the algorithms describes in the supplementary information for estimating the serial interval, effective reproduction number and their respective confidence intervals.

Code may run slowly due to the number of iterations (N = 10,000). 
