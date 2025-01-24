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

#########################################
## IL6 activity creation example

# binding affinities
k1 <- 0.5
k2 <- 0.05

il6 <- c(0.3,
         1.02,
         0.66,
         1.62)

il6r <- c(27.9887467638307,
          28.8646122214421,
          31.461509682922,
          32.05117)

il6st <- c(325.623471882641,
           346.01466992665,
           714.37125748503,
           496.963364512936)

data <- cbind(il6, il6r, il6st) %>% as.data.frame()

data <- data %>% 
  mutate(il6_nm = il6/(1000*23.7),
         il6r_nm = il6r/50,
         il6st_nm = il6st/100)

data <- data %>%
  mutate(
         bin_1 = (0.5*il6r_nm)+(0.5*il6_nm)+(0.5*k1),
         bin_2 = il6r_nm**2 + il6_nm**2 + 2*(il6_nm*k1) + k1**2) %>%
  mutate(bin_3 = 0.5*(bin_2**0.5)) %>%
  mutate(bin_complex = bin_1 - bin_3)

data <- data %>%
  mutate(tern_1 = (0.5*bin_complex)+(0.5*il6st_nm)+(0.5*k2),
         tern_2 = bin_complex**2 + il6st_nm**2 + 2*(bin_complex*k2) + k2**2) %>%
  mutate(tern_3 = 0.5*(tern_2**0.5)) %>%
  mutate(tern_complex = tern_1 - tern_3) %>%
  mutate(il6_activity = bin_complex/tern_complex)


