library(readr)
library(tidyverse)
library(TwoSampleMR)
library(ieugwasr)
library(ggforestplot)
library(janitor)

setwd("C:/Users/ss23664/OneDrive - University of Bristol/Thesis/IL6_bioavailability/Modelling_snps")

# Where to save outputs
output_folder <- "C:/Users/ss23664/OneDrive - University of Bristol/Thesis/IL6_bioavailability/Modelling_snps"
#

p_val_thresh = 5e-8

options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')


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
### IL6 ###
# -----------------------------------------------------#



#Uk Biobank
il6_ukb <- read.table("IL6/IL6_signif.txt", 
                      header = F, 
                      sep="")

il6_ukb <- data.frame(setNames(il6_ukb, ukb_vars)) %>%
  mutate(olink_target_fullname = "IL6_UKB")


# Format data
il6_ukb <- format_data(
  il6_ukb,
  type = "exposure",
  snps = NULL,
  header = TRUE,
  phenotype_col = "olink_target_fullname",
  snp_col = "rsid",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "P",
  samplesize_col = "N",
  id_col = "HGNC.symbol",
  chr_col = "CHROM",
  pos_col = "POS38",
  gene_col = "ensembl_id",
  eaf_col= "A1FREQ"
)

#il6_ukb <- clump_data(il6_ukb)

#il6_ukb$rsid<- il6_ukb$SNP
#il6_ukb$pval<- il6_ukb$pval.exposure



#GTEx - Blood
il6_gtex_blood <- read.table("IL6/IL6_gtex_blood.txt", 
                             header = T, 
                             sep="")

il6_gtex_blood %>% colnames()
il6_gtex_blood %>% head()

il6_gtex_blood %>%
  filter(p < (0.05/333)) %>%
  nrow()
# 0 

il6_gtex_blood <- il6_gtex_blood %>%
  #filter(p < (0.05/333)) %>% 
  mutate(phenotype = "IL6_gtex_blood") %>%
  mutate(N = NA_real_)


# Format data
il6_gtex_blood <- format_data(
  il6_gtex_blood,
  type = "exposure",
  snps = NULL,
  header = TRUE,
  phenotype_col = "phenotype",
  snp_col = "SNP",
  beta_col = "b",
  se_col = "SE",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "p",
  samplesize_col = "N",
  id_col = "phenotype",
  chr_col = "Chr",
  pos_col = "BP",
  gene_col = "Gene",
  eaf_col= "Freq"
)



#GTEx - Liver
il6_gtex_liver <- read.table("IL6/IL6_liver.txt", 
                             header = T, 
                             sep="")

il6_gtex_liver %>% colnames()
il6_gtex_liver %>% head()

il6_gtex_liver %>%
  filter(p < (0.05/333)) %>%
  nrow()
# 0

il6_gtex_liver <- il6_gtex_liver %>%
  #filter(p < (0.05/333)) %>%
  mutate(phenotype = "IL6_gtex_liver") %>% 
  mutate(N = NA_real_)

# Format data
il6_gtex_liver <- format_data(
  il6_gtex_liver,
  type = "exposure",
  snps = NULL,
  header = TRUE,
  phenotype_col = "phenotype",
  snp_col = "SNP",
  beta_col = "b",
  se_col = "SE",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "p",
  samplesize_col = "N",
  id_col = "phenotype",
  chr_col = "Chr",
  pos_col = "BP",
  gene_col = "Gene",
  eaf_col= "Freq"
)

#il6st_gtex_liver <- clump_data(il6st_gtex_liver)


#GWAS - CHARGE
il6_charge <- read.csv("IL6/charge_gwas_il6.csv") %>% 
  janitor::clean_names() %>% mutate(phenotype = "IL6_GWAS")


# Keep SNPs with p-value below threshold
il6_charge <- il6_charge %>% 
  filter(p_value < p_val_thresh)

il6_charge %>% nrow()

# Format data
il6_charge <- format_data(
  il6_charge,
  type = "exposure",
  snps = NULL,
  header = TRUE,
  phenotype_col = "phenotype",
  snp_col = "snp_rs_id",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "p_value",
  samplesize_col = "n_discovery",
  #id_col = "HGNC.symbol",
  chr_col = "chromosome",
  pos_col = "position_bp_gr_ch37",
  gene_col = "gene_nearest_gene_s",
  eaf_col= "effect_allele_frequency"
)

# Clump data
#il6_charge <- clump_data(il6_charge)



# -----------------------------------------------------#
### Combine data ###
# -----------------------------------------------------#

compare_df_cols(il6_ukb, il6_gtex_blood, il6_gtex_liver, il6_charge) %>% View()


il6_exposure_data <- rbind(il6_ukb, il6_gtex_blood, il6_gtex_liver, il6_charge)


# -----------------------------------------------------#
### F-statistic cal ###
# -----------------------------------------------------#

exp_info <- il6_exposure_data %>% select(SNP, ends_with("exposure")) %>% 
  distinct() %>% 
  mutate(f = ((beta.exposure/se.exposure)^2),
         r2 = 2 * (beta.exposure^2) * eaf.exposure * (1-eaf.exposure)
  ) %>%
  group_by(id.exposure) %>%
  reframe(
    target      = unique(gene.exposure),
    exposure    = unique(exposure),
    id_exposure = unique(id.exposure),
    Nsnps       = n(),
    N           = round(median(samplesize.exposure), 0),
    R2          = round(sum(r2), 3),
    f_mean           = round(mean(f), 0)
  )



# Get instruments for IL6 in eqtlgen

il6_eqtlgen <- extract_outcome_data(il6_exposure_data$SNP, "eqtl-a-ENSG00000136244")

il6_dat <- harmonise_data(il6_exposure_data, il6_eqtlgen, action = 1)

#Results
il6_res <- mr(il6_dat, method_list = c("mr_wald_ratio", 
                                         "mr_ivw"))


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


# -----------------------------------------------------#
### ggforestplot ###
# -----------------------------------------------------#
# 

# Forest plot for each exposure

il6_forest_plot <- forestplot(il6_res %>%
                                 #filter(method == "Inverse variance weighted") %>% 
                                 arrange(exposure),
                               name = exposure,
                               estimate = b,
                               se = se,
                               pvalue = pval,
                               colour = method,
                               shape = method,
                               xlab = "",
                               #ylab = "Outcome",
                               title = "MR of IL6 from UKB and GTEx \n on IL6ST in Eqtlgen",
                               logodds = FALSE
) + theme(axis.text = element_text(size = 10),
          axis.title.x = element_text(size = 10))

il6_forest_plot


# -----------------------------------------------------#
### IL6R ###
# -----------------------------------------------------#


il6r_ukb <- read.table("IL6R/IL6R_signif.txt", 
                       header = F,
                       sep = "\t")

il6r_ukb <- data.frame(setNames(il6r_ukb, ukb_vars)) %>%
  mutate(olink_target_fullname = "IL6R_UKB")

# Format data
il6r_ukb <- format_data(
  il6r_ukb,
  type = "exposure",
  snps = NULL,
  header = TRUE,
  phenotype_col = "olink_target_fullname",
  snp_col = "rsid",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "P",
  samplesize_col = "N",
  id_col = "HGNC.symbol",
  chr_col = "CHROM",
  pos_col = "POS38",
  gene_col = "ensembl_id",
  eaf_col= "A1FREQ"
)

#il6r_ukb <- clump_data(il6r_ukb)



#GTEx - Blood
il6r_gtex_blood <- read.table("IL6R/IL6R_gtex_blood.txt", 
                              header = T, 
                              sep="")

il6r_gtex_blood %>% colnames()
il6r_gtex_blood %>% head()

il6r_gtex_blood %>%
  filter(p < (0.05/333)) %>%
  nrow()
# 143 

il6r_gtex_blood <- il6r_gtex_blood %>%
  filter(p < (0.05/333)) %>% 
  mutate(phenotype = "IL6R_gtex_blood") %>%
  mutate(N = NA_real_)


# Format data
il6r_gtex_blood <- format_data(
  il6r_gtex_blood,
  type = "exposure",
  snps = NULL,
  header = TRUE,
  phenotype_col = "phenotype",
  snp_col = "SNP",
  beta_col = "b",
  se_col = "SE",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "p",
  samplesize_col = "N",
  id_col = "phenotype",
  chr_col = "Chr",
  pos_col = "BP",
  gene_col = "Gene",
  eaf_col= "Freq"
)

#il6r_gtex_blood <- clump_data(il6r_gtex_blood)

il6r_exposures <- rbind(il6r_ukb, il6r_gtex_blood)


# Outcome data
outcome_il6r_eqtlgen <- extract_outcome_data(il6r_exposures$SNP, 
                                             "eqtl-a-ENSG00000160712") %>%
  mutate(outcome = "IL6R eqtlgen")

# Harmonise

il6r_dat <- harmonise_data(il6r_exposures, outcome_il6r_eqtlgen, action = 1)

#Results
il6r_res <- mr(il6r_dat, method_list = c("mr_wald_ratio", 
                                         "mr_ivw"))


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


# -----------------------------------------------------#
### ggforestplot ###
# -----------------------------------------------------#
# 

# Forest plot for each exposure

il6r_forest_plot <- forestplot(il6r_res %>%
                                 #filter(method == "Inverse variance weighted") %>% 
                                 arrange(exposure),
                               name = exposure,
                               estimate = b,
                               se = se,
                               pvalue = pval,
                               colour = method,
                               shape = method,
                               xlab = "",
                               #ylab = "Outcome",
                               title = "MR of IL6R from UKB and GTEx \n on IL6ST in Eqtlgen",
                               logodds = FALSE
) + theme(axis.text = element_text(size = 10),
          axis.title.x = element_text(size = 10))

il6r_forest_plot

#save_func(paste0(output_folder,"/MR_CRP/all_forest_plots_clumped.png"), all_forest_plots)



# -----------------------------------------------------#
### IL6ST ###
# -----------------------------------------------------#

il6st_ukb <- read.table("IL6ST/IL6ST_signif.txt", 
                        header = F, 
                        sep="\t")

il6st_ukb <- data.frame(setNames(il6st_ukb, ukb_vars)) %>%
  mutate(olink_target_fullname = "IL6ST_UKB")

# Format data
il6st_ukb <- format_data(
  il6st_ukb,
  type = "exposure",
  snps = NULL,
  header = TRUE,
  phenotype_col = "olink_target_fullname",
  snp_col = "rsid",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "P",
  samplesize_col = "N",
  id_col = "HGNC.symbol",
  chr_col = "CHROM",
  pos_col = "POS38",
  gene_col = "ensembl_id",
  eaf_col= "A1FREQ"
)

#il6st_ukb <- clump_data(il6st_ukb)



#GTEx - Blood
il6st_gtex_blood <- read.table("IL6ST/IL6ST_gtex_blood.txt", 
                               header = T, 
                               sep="")

il6st_gtex_blood %>% colnames()
il6st_gtex_blood %>% head()

il6st_gtex_blood %>%
  filter(p < (0.05/333)) %>%
  nrow()
# 3 

il6st_gtex_blood <- il6st_gtex_blood %>%
  filter(p < (0.05/333)) %>% 
  mutate(phenotype = "IL6ST_gtex_blood") %>% 
  mutate(N = NA_real_)


# Format data
il6st_gtex_blood <- format_data(
  il6st_gtex_blood,
  type = "exposure",
  snps = NULL,
  header = TRUE,
  phenotype_col = "phenotype",
  snp_col = "SNP",
  beta_col = "b",
  se_col = "SE",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "p",
  samplesize_col = "N",
  id_col = "phenotype",
  chr_col = "Chr",
  pos_col = "BP",
  gene_col = "Gene",
  eaf_col= "Freq"
)

#il6st_gtex_blood <- clump_data(il6st_gtex_blood)


#GTEx - Liver
il6st_gtex_liver <- read.table("IL6ST/IL6ST_liver.txt", 
                               header = T, 
                               sep="")

il6st_gtex_liver %>% colnames()
il6st_gtex_liver %>% head()

il6st_gtex_liver %>%
  filter(p < (0.05/333)) %>%
  nrow()
# 2

il6st_gtex_liver <- il6st_gtex_liver %>%
  filter(p < (0.05/333)) %>%
  mutate(phenotype = "IL6ST_gtex_liver") %>% 
  mutate(N = NA_real_)

# Format data
il6st_gtex_liver <- format_data(
  il6st_gtex_liver,
  type = "exposure",
  snps = NULL,
  header = TRUE,
  phenotype_col = "phenotype",
  snp_col = "SNP",
  beta_col = "b",
  se_col = "SE",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "p",
  samplesize_col = "N",
  id_col = "phenotype",
  chr_col = "Chr",
  pos_col = "BP",
  gene_col = "Gene",
  eaf_col= "Freq"
)

#il6st_gtex_liver <- clump_data(il6st_gtex_liver)


# Combine exposures
il6st_exposures <- rbind(il6st_ukb, il6st_gtex_blood, il6st_gtex_liver)


# Outcome data
outcome_il6st_eqtlgen <- extract_outcome_data(il6st_exposures$SNP, 
                                              "eqtl-a-ENSG00000134352") %>%
  mutate(outcome = "IL6ST eqtlgen")

# Harmonise

il6st_dat <- harmonise_data(il6st_exposures, outcome_il6st_eqtlgen, action = 1)

#Results
il6st_res <- mr(il6st_dat, method_list = c("mr_wald_ratio", "mr_ivw"))


# -----------------------------------------------------#
### Plots ###
# -----------------------------------------------------#

# Scatter plot
p1 <- mr_scatter_plot(il6st_res, il6st_dat)

length(p1)

p1[[1]]


# -----------------------------------------------------#
### ggforestplot ###
# -----------------------------------------------------#
# 

# Forest plot for each exposure

il6st_forest_plot <- forestplot(il6st_res %>%
                                 #filter(method == "Inverse variance weighted") %>% 
                                 arrange(exposure),
                               name = exposure,
                               estimate = b,
                               se = se,
                               pvalue = pval,
                               colour = method,
                               shape = method,
                               xlab = "",
                               #ylab = "Outcome",
                               title = "MR of IL6ST from UKB and GTEx \n on IL6ST in Eqtlgen",
                               logodds = FALSE
) + theme(axis.text = element_text(size = 10),
          axis.title.x = element_text(size = 10))

il6st_forest_plot

#save_func(paste0(output_folder,"/MR_CRP/all_forest_plots_clumped.png"), all_forest_plots)


