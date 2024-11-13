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

# P-value threshold
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

il6_eqtlgen <- read.csv("IL6/IL6_eqtlgen.csv") %>% 
  mutate(EXPOSURE = paste0(EXPOSURE, "_eqtlgen"))

# Format data
il6_eqtlgen <- format_data(
  il6_eqtlgen,
  type = "exposure",
  snps = NULL,
  header = TRUE,
  phenotype_col = "EXPOSURE",
  snp_col = "RSID",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "EA",
  other_allele_col = "OA",
  pval_col = "Pvalue",
  samplesize_col = "N",
  id_col = "EXPOSURE",
  chr_col = "GeneChr",
  pos_col = "BP",
  gene_col = "Gene",
  eaf_col= "EAF"
)


## -- Clump exposure data ####
options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')

il6_eqtlgen <- clump_data(il6_eqtlgen)


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


#Clumping
# 
# exposure_clumped<- ieugwasr::ld_clump(
#   dat = il6_ukb_exposure,
#   clump_kb = 10000,
#   clump_r2 = 0.01,
#   clump_p = 0.99,
#   pop = "EUR",
#   bfile = "/path/REF_PANEL/EUR",
#   plink_bin ="/hpcpath/plink-1.90/plink")
# 
# 


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
### IL6R ###
# -----------------------------------------------------#

il6r_eqtlgen <- read.csv("IL6R/IL6R_eqtlgen.csv") %>%
  mutate(EXPOSURE = paste0(EXPOSURE, "_eqtlgen"))

# Format data
il6r_eqtlgen <- format_data(
  il6r_eqtlgen,
  type = "exposure",
  snps = NULL,
  header = TRUE,
  phenotype_col = "EXPOSURE",
  snp_col = "RSID",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "EA",
  other_allele_col = "OA",
  pval_col = "Pvalue",
  samplesize_col = "N",
  id_col = "EXPOSURE",
  chr_col = "GeneChr",
  pos_col = "BP",
  gene_col = "Gene",
  eaf_col= "EAF"
)

il6r_eqtlgen <- clump_data(il6r_eqtlgen)


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

il6r_ukb <- clump_data(il6r_ukb)



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
  mutate(N = NA)


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

il6r_gtex_blood <- clump_data(il6r_gtex_blood)


#GTEx - Liver
il6r_gtex_liver <- read.table("IL6R/IL6R_liver.txt", 
                             header = T, 
                             sep="")

il6r_gtex_liver %>% colnames()
il6r_gtex_liver %>% head()

il6r_gtex_liver %>%
  filter(p < (0.05/333)) %>%
  nrow()

# 0



# -----------------------------------------------------#
### IL6ST ###
# -----------------------------------------------------#

il6st_eqtlgen <- read.csv("IL6ST/IL6ST_eqtlgen.csv") %>% 
  mutate(EXPOSURE = paste0(EXPOSURE, "_eqtlgen"))

# Format data
il6st_eqtlgen <- format_data(
  il6st_eqtlgen,
  type = "exposure",
  snps = NULL,
  header = TRUE,
  phenotype_col = "EXPOSURE",
  snp_col = "RSID",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "EA",
  other_allele_col = "OA",
  pval_col = "Pvalue",
  samplesize_col = "N",
  id_col = "EXPOSURE",
  chr_col = "GeneChr",
  pos_col = "BP",
  gene_col = "Gene",
  eaf_col= "EAF"
)

il6st_eqtlgen <- clump_data(il6st_eqtlgen)


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

il6st_ukb <- clump_data(il6st_ukb)



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
  mutate(N = NA)


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

il6st_gtex_blood <- clump_data(il6st_gtex_blood)


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
  mutate(N = NA)

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

il6st_gtex_liver <- clump_data(il6st_gtex_liver)

# -----------------------------------------------------#
### ADAM17 ###
# -----------------------------------------------------#

adam17_eqtlgen <- read.csv("ADAM17/ADAM17_eqtlgen.csv", header = F)
adam17_eqtlgen <- data.frame(setNames(adam17_eqtlgen, eqtlgen_vars)) %>% 
  mutate(EXPOSURE = paste0(EXPOSURE, "_eqtlgen"))

# Format data
adam17_eqtlgen <- format_data(
  adam17_eqtlgen,
  type = "exposure",
  snps = NULL,
  header = TRUE,
  phenotype_col = "EXPOSURE",
  snp_col = "RSID",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "EA",
  other_allele_col = "OA",
  pval_col = "Pvalue",
  samplesize_col = "N",
  id_col = "EXPOSURE",
  chr_col = "GeneChr",
  pos_col = "BP",
  gene_col = "Gene",
  eaf_col= "EAF"
)

adam17_eqtlgen <- clump_data(adam17_eqtlgen)


# GWAS
adam17_gwas <- read.csv("ADAM17/adam17_gwas.csv")

adam17_gwas %>% colnames()
adam17_gwas %>% View()

adam17_gwas <- adam17_gwas %>%
  filter(Exposure == "extracellular") %>%
  mutate(phenotype = "ADAM17_gwas",
         gene = NA)

# Format data
adam17_gwas <- format_data(
  adam17_gwas,
  type = "exposure",
  snps = NULL,
  header = TRUE,
  phenotype_col = "phenotype",
  snp_col = "SNP",
  beta_col = "Coefficient",
  se_col = "SE",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "P",
  samplesize_col = "N",
  id_col = "phenotype",
  chr_col = "Chr",
  pos_col = "Pos",
  gene_col = "gene",
  eaf_col= "EAF"
)

adam17_gwas <- clump_data(adam17_gwas)




#GTEx - Blood
adam17_gtex_blood <- read.table("ADAM17/ADAM17_gtex_blood.txt", 
                             header = T, 
                             sep="")

adam17_gtex_blood %>% colnames()
adam17_gtex_blood %>% head()

adam17_gtex_blood %>%
  filter(p < (0.05/333)) %>%
  nrow()
# 0 


#GTEx - Liver
adam17_gtex_liver <- read.table("ADAM17/ADAM17_liver.txt", 
                             header = T, 
                             sep="")

adam17_gtex_liver %>% colnames()
adam17_gtex_liver %>% head()

adam17_gtex_liver %>%
  filter(p < (0.05/333)) %>%
  nrow()
# 1

adam17_gtex_liver <- adam17_gtex_liver %>%
  filter(p < (0.05/333)) %>% 
  mutate(phenotype = "adam17_gtex_liver") %>% 
  mutate(N = NA)


# Format data
adam17_gtex_liver <- format_data(
  adam17_gtex_liver,
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



# -----------------------------------------------------#
### Combine data ###
# -----------------------------------------------------#

compare_df_cols(il6_eqtlgen, il6_ukb, il6_charge,
                il6r_eqtlgen, il6r_ukb, il6r_gtex_blood,
                il6st_eqtlgen, il6st_ukb, il6st_gtex_blood, il6st_gtex_liver,
                adam17_eqtlgen, adam17_gwas, adam17_gtex_liver) %>% View()


all_exposure_data <- rbind(il6_eqtlgen, il6_ukb, il6_charge,
                           il6r_eqtlgen, il6r_ukb, il6r_gtex_blood,
                           il6st_eqtlgen, il6st_ukb, il6st_gtex_blood, il6st_gtex_liver,
                           adam17_eqtlgen, adam17_gwas, adam17_gtex_liver)


# -----------------------------------------------------#
### F-statistic cal ###
# -----------------------------------------------------#

exp_info <- all_exposure_data %>% select(SNP, ends_with("exposure")) %>% 
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

write.csv(exp_info, paste0(output_folder,"/MR_CRP/exposure_info.csv"))

# -----------------------------------------------------#
### CRP ###
# -----------------------------------------------------#
# crp_snps <- read_tsv('C:/Users/ss23664/Downloads/35459240-GCST90029070-EFO_0004458-Build37.tsv')
# 
# crp_snps %>% colnames()
# crp_snps %>% nrow()
# # 11106737
# 
# crp_snps %>% head() %>% View()
# 
# crp_snps %>% 
#   summarise(min(odds_ratio), max(odds_ratio), 
#             min(ci_lower), max(ci_lower), 
#             min(effect_allele_frequency), max(effect_allele_frequency)) %>% 
#   View()
# 
# crp_snps <- crp_snps %>% 
#   select(-c(odds_ratio, ci_lower, ci_upper)) %>%
#   filter(variant_id %in% all_exposure_data$SNP)
# 
# 
# # Re-format outcome dataset
# outcome_dat_crp <- format_data(crp_snps,
#                            type = "outcome",
#                            snp_col = "variant_id",
#                            beta_col = "beta",
#                            se_col = "standard_error",
#                            effect_allele_col = "effect_allele",
#                            other_allele_col = "other_allele",
#                            eaf_col = "effect_allele_frequency",
#                            pval_col = "p_value",
#                            pos_col = "base_pair_location"
#                            
# )

# Simplified extraction method

outcome_dat_crp <- extract_outcome_data(all_exposure_data$SNP, "ebi-a-GCST90029070")


# -----------------------------------------------------#
### Check to see which exposure SNPs are in outcome ###
# -----------------------------------------------------#

all_exposure_data$SNP %in% outcome_dat_crp$SNP
all_exposure_data$SNP[!(all_exposure_data$SNP %in% outcome_dat_crp$SNP)]
# c("rs115153905", "rs2894379", "rs139564096", "rs9272460", 
# "rs6457457", "rs444921", "rs56086522")

all_exposure_data %>% 
  filter(SNP %in% c("rs115153905", "rs2894379", "rs139564096", 
                    "rs9272460", "rs6457457", "rs444921", "rs56086522")) %>% 
  select(SNP, exposure) %>%
  View()




# -----------------------------------------------------#
### MR ###
# -----------------------------------------------------#
dat <- harmonise_data(
  exposure_dat = all_exposure_data,
  outcome_dat = outcome_dat_crp
)

write.csv(dat, paste0(output_folder,"/MR_CRP/harmonised_data.csv"))

# -----------------------------------------------------#
### Read harmonised data ###
# avoids running previous code
# -----------------------------------------------------#

dat <- read.csv(paste0(output_folder,"/MR_CRP/harmonised_data.csv"))

dat %>% colnames()


# Look at list of exposures
unique(dat$exposure)

# Check all SNPs have same alleles in exposure and outcome datasets
dat %>% 
  mutate(check_effect_allele = 
           ifelse(effect_allele.exposure == effect_allele.outcome, 0, 1),
         check_other_allele = 
           ifelse(other_allele.exposure == other_allele.outcome, 0, 1)) %>%
  summarise(sum(check_effect_allele), sum(check_other_allele))

# All MR methods
#res <- mr(dat)

mr_method_list() %>% select(obj, Description) %>% View()


# Select only MR Egger, weighted median, and Inverse Variance Weighted (IVW)
res <- mr(dat, method_list = c("mr_wald_ratio", "mr_egger_regression", 
                               "mr_weighted_median", "mr_ivw"))

# Heterogeneity test
mr_heterogeneity(dat, method_list = c("mr_ivw")) %>% 
  select(outcome, Q, Q_df, Q_pval)

# Egger intercept term
mr_pleiotropy_test(dat) %>% collect() %>% print()


# -----------------------------------------------------#
### Table of results ###
# -----------------------------------------------------#

outcome_table_all_methods <- res %>% 
  mutate(outcome = "CRP") %>%
  #filter(method == "Inverse variance weighted") %>% 
  select(exposure, outcome, method, nsnp, b, se, pval) %>% 
  arrange(exposure, outcome, method)


write.csv(outcome_table_all_methods, paste0(output_folder,"/MR_CRP/results_all_MR_methods_clumped.csv"))


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


p1 <- mr_scatter_plot(res, dat)

length(p1)

p1[[1]]


for (i in 1:length(p1)){
  save_func(paste0(output_folder,"/MR_CRP/", 
                   unique(res$exposure)[i], "scatter_plot.png"), p1[[i]])
}

# Single SNP MR using the Wald ratio
res_single <- mr_singlesnp(dat)

# -----------------------------------------------------#
### ggforestplot ###
# -----------------------------------------------------#
# 

# Forest plot for each exposure

all_forest_plots <- forestplot(res %>%
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
                          title = "MR of proteins involved in IL6 signalling \n on CRP levels",
                          logodds = FALSE
) + theme(axis.text = element_text(size = 10),
          axis.title.x = element_text(size = 10))

all_forest_plots

save_func(paste0(output_folder,"/MR_CRP/all_forest_plots_clumped.png"), all_forest_plots)


# -----------------------------------------------------#
### IL6 ###
# -----------------------------------------------------#

il6_snps <- read.csv("IL6/charge_gwas_il6.csv") %>% 
  janitor::clean_names()

il6_snps$snp_rs_id %in% crp_snps$variant_id
il6_snps$snp_rs_id[!(il6_snps$snp_rs_id %in% crp_snps$variant_id)]
il6_snps$snp_rs_id[il6_snps$snp_rs_id %in% crp_snps$variant_id]


crp_filtered_il6 <- crp_snps %>%
  filter(variant_id %in% il6_snps$snp_rs_id) %>%
  mutate(gene = "IL6")


# -----------------------------------------------------#
### IL6R ###
# -----------------------------------------------------#

il6r_snps <- read.csv("IL6R/il6r_snps.csv") %>% 
  janitor::clean_names()

il6r_snps %>% colnames()

il6r_snps$snp %in% crp_snps$variant_id
il6r_snps$snp[!(il6r_snps$snp %in% crp_snps$variant_id)]
il6r_snps$snp[il6r_snps$snp %in% crp_snps$variant_id]


crp_filtered_il6r <- crp_snps %>%
  filter(variant_id %in% il6r_snps$snp) %>%
  mutate(gene = "IL6R")

