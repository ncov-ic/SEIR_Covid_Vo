## packages needed
list_of_packages <- c("tidyr",
                      "dplyr",
                      "readr",
                      "janitor",
                      "lubridate",
                      "magrittr",
                      "ggplot2",
                      "plotly",
                      "htmlwidgets",
                      "epicontacts",
                      "igraph",
                      "GGally",
                      "epitrix",
                      "ggpubr",
                      "flextable",
                      "network",
                      "sna",
                      "scales",
                      "intergraph",
                      "cowplot")

## install missing packages
new_packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if(length(new_packages) > 0) {
  install.packages(new.packages)
}

## load packages 
library(tidyr)
library(dplyr)
library(readr)
library(janitor)
library(lubridate)
library(magrittr)
library(ggplot2)
library(plotly)
library(htmlwidgets)
library(epicontacts)
library(igraph)
library(GGally)
library(epitrix)
library(ggpubr)
library(flextable)
library(network)
library(sna)
library(scales)
library(intergraph)
library(cowplot)

## load source codes
source("R/serial_interval_functions.R")
source("R/reproduction_number_functions.R")
source("R/checking_functions.R")
source("R/plot_clusters.R")
source("R/plot_tree.R")

## load data - these are the outputs of italy_cleaning_data_script.R
full_data <- readRDS("full_data.rds")
linelist <- readRDS("linelist.rds")
all_contacts <- readRDS("all_contact.rds")
direct_contacts <- readRDS("direct_contact.rds")


## set number of iterations
iteration_value <- 10000


contacts <- all_contacts[,c("id", "contact_id")]
names(contacts) <- c("to", "from")

## find individuals for whom we know the delay
known_delay <- linelist %>%
  dplyr::filter(!is.na(onset) == TRUE) %>%
  dplyr::mutate(onset_to_confirmation = as.numeric(difftime(first_positive, onset, units = "days")))

onset_to_confirmation <- known_delay$onset_to_confirmation

## test sampling from this set of values
# sample(onset_to_confirmation, size = 50, replace = TRUE)

## just the cases of villagers
linelist_villager <- linelist %>%
  dplyr::filter(!is.na(first_positive) == TRUE)

## I know there should be 80 entries
if(nrow(linelist_villager) != 80) {
  stop("This dataset does not include all inhabitants of Vo' +ve for SARS-CoV-2")
}

## remove any duplicated entries 
df_sort <- t(apply(contacts, 1, sort))
contacts <- contacts[!duplicated(df_sort),]

## contacts between villagers
village_contacts <- contacts %>%
  dplyr::filter((to %in% linelist_villager$id) == TRUE) %>%
  dplyr::filter((from %in% linelist_villager$id) == TRUE)

### serial interval - whole study period 
central_SI <- central_serial_interval(linelist = linelist_villager, 
                                      contacts = village_contacts, 
                                      onset_delay = onset_to_confirmation,
                                      iterations = iteration_value)

## central SI values
central_mean <- mean(central_SI$mean)
central_variance <- mean(central_SI$variance)

## derive shape and scale parameters from their mean and variance 
central_shape <- (central_mean)^2/central_variance
central_scale <- central_variance/central_mean

# 95% confidence interval on the central estimate for the serial interval:
confidence_interval <- confidence_serial_interval(linelist = linelist_villager, 
                                                  contacts = village_contacts, 
                                                  onset_delay = onset_to_confirmation,
                                                  iterations = iteration_value)

lower_ci <- quantile(confidence_interval$mean, c(.025, .975))[1]
upper_ci <- quantile(confidence_interval$mean, c(.025, .975))[2]

## Sensitivity Analysis - Confidence Interval 
# central_sensitivity_analysis <- central_sensitivity(iterations = 10000)
central_sensitivity_analysis <- central_sensitivity(linelist = linelist_villager, 
                                                    contacts = village_contacts, 
                                                    onset_delay = onset_to_confirmation, 
                                                    iterations = iteration_value)

pre_lockdown_si <- central_sensitivity_analysis %>%
  dplyr::filter(status == "pre-lockdown") 

post_lockdown_si <- central_sensitivity_analysis %>%
  dplyr::filter(status == "post-lockdown") 

pre_mean <- mean(pre_lockdown_si$mean)
post_mean <- mean(post_lockdown_si$mean)

pre_variance <- mean(pre_lockdown_si$variance)
post_variance <- mean(post_lockdown_si$variance)

## parameters for the pre and post lockdown distributions 
pre_parm <- data.frame(mean = pre_mean,
                       shape = (pre_mean^2)/pre_variance,
                       scale = pre_variance/pre_mean)

post_parm <- data.frame(mean = post_mean,
                       shape = (post_mean^2)/post_variance,
                       scale = post_variance/post_mean)

## 95% CI on the central estimates pre and post lockdown
confidence_sensitivity_analysis <- confidence_sensitivity(linelist = linelist_villager, 
                                                          contacts = village_contacts, 
                                                          onset_delay = onset_to_confirmation, 
                                                          iterations = iteration_value)

pre_lockdown_confidence <- central_sensitivity_analysis %>%
  dplyr::filter(status == "pre-lockdown") %>% 
  dplyr::summarise(lower_si = quantile(mean, 0.025),
                   upper_si = quantile(mean, 0.975))

post_lockdown_confidence <- confidence_sensitivity_analysis %>%
  dplyr::filter(status == "post-lockdown") %>% 
  dplyr::summarise(lower_si = quantile(mean, 0.025),
                   upper_si = quantile(mean, 0.975))

## plot the serial interval
si_df <- data.frame(days = seq(0, 30, length.out = 10000),
                    overall = dgamma(x = seq(0, 30, length.out = 10000), 
                                     shape = central_shape, 
                                     scale = central_scale),
                    pre_lockdown = dgamma(x = seq(0, 30, length.out = 10000), 
                                          shape = pre_parm$shape, 
                                          scale = pre_parm$scale),
                    post_lockdown = dgamma(x = seq(0, 30, length.out = 10000), 
                                           shape = post_parm$shape, 
                                           scale = post_parm$scale))

si_df <- tidyr::pivot_longer(si_df, overall:post_lockdown,
                             names_to = "duration", values_to = "mle")

si_df$duration <- factor(si_df$duration, levels = c("overall", 
                                                    "pre_lockdown",
                                                    "post_lockdown"))

p1 <- ggplot(si_df, aes(x = days, y = mle, col = duration)) + 
  geom_line() + 
  theme_bw() + 
  xlim(c(0, 15)) +
  labs(x = "serial interval (days)", y = "probability")+
  theme(text=element_text(size=7, family="Sans")) + 
  theme(legend.title = element_blank()) + 
  theme(legend.position = "bottom")

ggsave("serial_interval.png", 
       plot = p1, 
       width = 8.9, 
       height = 5, 
       units = "cm", 
       dpi = 500)


## Estimating Rt

#Central estimate:

effective_reproductive_number <- central_rep_num(linelist = linelist_villager, 
                                                 contacts = village_contacts, 
                                                 onset_delay = onset_to_confirmation, 
                                                 iterations = iteration_value, 
                                                 pre_shape = pre_parm$shape, 
                                                 pre_scale = pre_parm$scale,
                                                 post_shape = post_parm$shape, 
                                                 post_scale = post_parm$scale)

cohort_1 <- effective_reproductive_number %>% dplyr::filter(cohort == "1")
cohort_2 <- effective_reproductive_number %>% dplyr::filter(cohort == "2")

mean(cohort_1$R)
mean(cohort_2$R)

#95% confidence interval:

## final long run
bootstrapped_R <- confidence_rep_num(linelist = linelist_villager, 
                                     contacts = village_contacts, 
                                     onset_delay = onset_to_confirmation, 
                                     iterations = iteration_value, 
                                     pre_shape = pre_parm$shape, 
                                     pre_scale = pre_parm$scale,
                                     post_shape = post_parm$shape, 
                                     post_scale = post_parm$scale)

cohort_1_boot <- bootstrapped_R %>%
  dplyr::filter(cohort == "1")
cohort_2_boot <- bootstrapped_R %>%
  dplyr::filter(cohort == "2")

mean(cohort_1$R)
lower_ci_R_1 <- quantile(cohort_1_boot$R, c(0.025)) 
upper_ci_R_1 <- quantile(cohort_1_boot$R, c(0.975)) 

mean(cohort_2$R)
lower_ci_R_2 <- quantile(cohort_2_boot$R, c(0.025)) 
upper_ci_R_2 <- quantile(cohort_2_boot$R, c(0.975)) 


## Cluster Analysis

linelist_village <- cluster_add_func(linelist_villager, village_contacts)
plot_clusters(linelist_village, village_contacts)

p2 <- plot_clusters(linelist_village, village_contacts)

ggsave("clusters.png", dpi = 500, width = 8.9, height = 5)

##  multifigure plot
plot_theme <- theme(legend.position = "none",
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank(),
                    panel.spacing.x = unit(0.2, "lines"),
                    panel.spacing.y = unit(0.2, "lines"),
                    plot.margin = unit(c(0, 0, 0, 0), "cm"),
                    legend.title = element_blank(),
                    text=element_text(size=7, family="Sans"))
## serial interval figure
p1 <- ggplot(si_df, aes(x = days, y = mle, col = duration)) + 
  geom_line() + 
  theme_bw() + 
  xlim(c(0, 15)) +
  labs(x = "serial interval (days)", y = "probability") + 
  plot_theme

legend <- get_legend(
  p1 +
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom",
          legend.direction = "horizontal",
          legend.key=element_blank(),
          legend.margin=margin(0,0,0,0),
          legend.text = element_text(size = 5, family = "Sans"),
          legend.box.spacing = unit(0, "cm"),
          legend.key.size = grid::unit(0.5, "lines"))
)

p1_leg <- plot_grid(p1, legend, labels = "", nrow = 2, 
                    rel_heights = c(1, 0.2))

tiff("si_clusters.tiff", width=8.9, height=5, units="cm", res=500)
plot_grid(p1_leg,p2, labels = c('a', 'b'), ncol=2, label_size = 8, label_fontface = "bold")
dev.off()
