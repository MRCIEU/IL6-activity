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

il6_ukb <- clump_data(il6_ukb)

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
il6_charge <- clump_data(il6_charge)



# -----------------------------------------------------#
### Combine data ###
# -----------------------------------------------------#

compare_df_cols(il6_ukb, il6_charge) %>% View()


il6_exposure_data <- rbind(il6_ukb, il6_charge)


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

dat <- harmonise_data(il6_exposure_data, il6_eqtlgen)
res <- mr(dat)
