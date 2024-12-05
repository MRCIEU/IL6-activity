# IL6-activity
Modelling the effect of IL-6 bioavailability/activity using genetic variants.

These scripts look at initial work to determine the viability of deriving a measure of IL-6 activity using genetic variants. 


To do this, we look at four proteins involved in IL-6 signalling - IL-6, IL-6R, IL6ST, ADAM17 - and measure their causal effect on CRP levels using Mendelian Randomisation (MR). Since CRP is downstream of these proteins in the IL-6 signalling process, we expect to see some causal associations. Bidirectional MR is performed with the hypothesis that CRP should not have a direct effect on the levels of these proteins. This is done with numerous data sources per protein where data is available.

#Contents:

# mr_crp.R
This takes, separately, summary level genetic data on IL6, IL6R, IL6ST, and ADAM17 as the exposures and a GWAS from Said et al. (2022) for genetic variants associated with CRP levels as the outcome. For the exposures, we look at data from the UK Biobank, GTEx (whole blood and liver), Eqtlgen, and appropriate GWAS.

# crp_as_exposure.R
This script performs MR of genes associated with CRP levels from Said et al. on to Eqtlgen datasets for IL6, IL6R, IL6ST, ADAM17.

We require the following files:

1. UK Biobank - summary data of pro