library(readr)
library(tidyverse)
library(TwoSampleMR)
library(ieugwasr)
library(ggforestplot)
library(janitor)
library(MRInstruments)

setwd("C:/Users/ss23664/OneDrive - University of Bristol/Thesis/IL6_bioavailability/Modelling_snps")

# Where to save outputs
output_folder <- "C:/Users/ss23664/OneDrive - University of Bristol/Thesis/IL6_bioavailability/Modelling_snps"
#

p_val_thresh = 5e-8


# -----------------------------------------------------#
### Dataset variable names ###
# -----------------------------------------------------#

# Biobank column names
ukb_vars <- read.table("ukb_clean_vars.txt", 
                       header = T, 
                       sep="") %>% 
  colnames()

# eqtlgen column names
eqtlgen_vars <- read.csv("eqtl_vars.csv", 
                         header = T) %>%
  colnames()


# -----------------------------------------------------#
### IEU GWAS database ###
# -----------------------------------------------------#

options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')

ao <- available_outcomes()
ao %>% colnames(0)
ao %>% head()

ao %>% select(trait) %>% unique() %>% nrow()
ao %>% select(trait) %>% unique() %>% as.data.frame() %>% head(50)

# Patterns to filter
patterns <- c("IL6", "il6", "il-6", "IL-6", 
              "Interleukin-6", "gp130", "ADAM17", 
              "ENSG00000136244", #IL6
              "ENSG00000160712", #IL6R
              "ENSG00000134352", #IL6ST
              "ENSG00000151694" #ADM17
              )

ao %>% select(trait) %>% unique() %>% as.data.frame() %>% 
  filter(grepl("IL6", trait)) %>% nrow()

filtered_ao <- ao %>% filter(grepl(paste(patterns, collapse='|'), trait))

# -----------------------------------------------------#
### CRP ###
# -----------------------------------------------------#

crp_snps <- extract_instruments("ebi-a-GCST90029070")

crp_snps %>% glimpse()
crp_snps %>% head()

crp_snps %>% summarise(min_pval = min(pval.exposure), max_pval = max(pval.exposure))

# Clump SNPs
crp_snps <- clump_data(crp_snps)

# -----------------------------------------------------#
### IL6 ###
# -----------------------------------------------------#

outcome_il6_eqtlgen <- extract_outcome_data(crp_snps$SNP, 
                                            "eqtl-a-ENSG00000136244") %>%
  mutate(outcome = "IL6 eqtlgen")


# -----------------------------------------------------#
### IL6R ###
# -----------------------------------------------------#

outcome_il6r_eqtlgen <- extract_outcome_data(crp_snps$SNP, 
                                             "eqtl-a-ENSG00000160712") %>%
  mutate(outcome = "IL6R eqtlgen")


# -----------------------------------------------------#
### IL6ST ###
# -----------------------------------------------------#

outcome_il6st_eqtlgen <- extract_outcome_data(crp_snps$SNP, 
                                              "eqtl-a-ENSG00000134352") %>%
  mutate(outcome = "IL6ST eqtlgen")


# -----------------------------------------------------#
### ADAM17 ###
# -----------------------------------------------------#

outcome_adam17_eqtlgen <- extract_outcome_data(crp_snps$SNP, 
                                               "eqtl-a-ENSG00000151694") %>%
  mutate(outcome = "ADAM17 eqtlgen")


# -----------------------------------------------------#
### All outcomes ###
# -----------------------------------------------------#

all_outcome_data <- rbind(outcome_il6_eqtlgen,
                          outcome_il6r_eqtlgen,
                          outcome_il6st_eqtlgen,
                          outcome_adam17_eqtlgen)


#Harmonise
crp_dat <- harmonise_data(crp_snps, all_outcome_data)
crp_dat <- crp_dat %>% mutate(exposure = "C-reactive protein levels")

# 
#Results
crp_res <- mr(crp_dat, method_list = c("mr_wald_ratio", "mr_egger_regression", 
                                       "mr_weighted_median", "mr_ivw"))


# -----------------------------------------------------#
### Plots ###
# -----------------------------------------------------#

#save
save_func <- function(file_name, plot_name)
{
  png(file_name, res=450, height=5500, width=3500)
  print(plot_name)
  dev.off()
  
}

# Scatter plot
p1 <- mr_scatter_plot(crp_res, crp_dat)

length(p1)

p1[[1]]


# -----------------------------------------------------#
### ggforestplot ###
# -----------------------------------------------------#
# 

# Forest plot for each exposure

crp_forest_plot <- forestplot(crp_res %>%
                                 #filter(method == "Inverse variance weighted") %>% 
                                 arrange(exposure),
                               name = outcome,
                               estimate = b,
                               se = se,
                               pvalue = pval,
                               colour = method,
                              shape = method,
                               xlab = "",
                               #ylab = "Outcome",
                               title = "MR of CRP levels on proteins involved \n in IL6 signalling",
                               logodds = FALSE
) + theme(axis.text = element_text(size = 10),
          axis.title.x = element_text(size = 10))

crp_forest_plot

save_func(paste0(output_folder,"/MR_CRP/CRP_exposure_forest_plot.png"), crp_forest_plot)


