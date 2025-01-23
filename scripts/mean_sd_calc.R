library(readr)
library(tidyverse)
library(TwoSampleMR)
library(ieugwasr)
library(ggforestplot)
library(janitor)

setwd("C:/Users/ss23664/OneDrive - University of Bristol/Thesis/IL6_bioavailability/Biobank")

# Where to save outputs
output_folder <- "C:/Users/ss23664/OneDrive - University of Bristol/Thesis/IL6_bioavailability/Biobank"
#

p_val_thresh = 5e-8

options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')


# -----------------------------------------------------#
### IL6 data ###
# -----------------------------------------------------#

il6_data <- read.csv("il6_levels_datasets.csv") %>% 
  select(-c(il6)) %>%
  filter(!is.na(il6_mean), 
         female_perc != 0 & female_perc != 100,
         il6_mean < 30,
         il6_sd < 10,
         age >= 40 & age < 70)

il6_data %>% colnames()
il6_data %>% head()

il6_data <- il6_data %>% 
  select(name, n, age, female_perc, il6_mean, il6_sd)

il6_data %>% View()

grand.mean <- function(m, n) {weighted.mean(m, n)}

grand.sd   <- function(S, M, N) {sqrt(weighted.mean(S^2 + M^2, N) -
                                        weighted.mean(M, N)^2)}

grand.mean(il6_data$il6_mean, il6_data$n)
grand.sd(il6_data$il6_sd, il6_data$il6_mean, il6_data$n)
