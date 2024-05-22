## Mendelian randomization analysis using MRSample
# EVS 5/7/2024

#library(remotes)
library(tidyverse)
#install_github("MRCIEU/TwoSampleMR")
#remotes::install_github('MRCIEU/ieugwasr')
library("TwoSampleMR")
library("ieugwasr")

## the built-in server access is broken; this must be fed manually for each command
# details: https://mrcieu.github.io/ieugwasr/articles/guide.html#authentication
mytoken <- "/mytokencharactersfromAPIthingy"

## ---- outcomes ----

# get outcome data - runs for a while
ao <- TwoSampleMR::available_outcomes(opengwas_jwt = mytoken)

## get outcomes related to ischemia and coronary heart disease
myoutcomes <- c("Major coronary heart disease event", "Emergency coronary revascularization (for ACS) (no controls excluded)", "Coronary artery disease (SPA correction)", "Coronary artery disease (Firth correction)", "Major coronary heart disease event excluding revascularizations", "Coronary heart disease", "Coronary atherosclerosis (no controls excluded)", "Coronary artery disease", "Ischemic Heart Disease", "Ischemic heart diseases")

# subset
cvd <- ao %>% 
  filter(trait %in% myoutcomes)

## ---- test just the heart snps ----


# get Kazachstania
kaz <- read.table("gwas-out/updateITS_results_allsnp_allvariant/allsnp_allvariants.Genus_Kazachstania.glm.linear", sep = "\t", header = TRUE, comment.char = "")

# adjust significance
kaz <- kaz %>% 
  mutate(pfdr = p.adjust(P, method = "fdr"))

# format for package
instr1 <- TwoSampleMR::format_data(kaz %>% filter(ID %in% c("rs12149890")), 
                                   snp_col = "ID",
                                   beta_col = "BETA",
                                   se_col = "SE",
                                   eaf_col = "A1_FREQ",
                                   effect_allele_col = "ALT",
                                   other_allele_col = "REF",
                                   pval_col = "P"
) %>% mutate(exposure = "Kazachstania")

# get outcomes
outc1 <- extract_outcome_data(
  snps = instr1$SNP,
  outcomes = unique(cvd$id),
  opengwas_jwt = mytoken
)

mydat1 <- outc1 %>% filter(outcome %in% tokeep)

# harmonize data
hdat1 <- harmonise_data(exposure_dat = instr1,
                        outcome_dat = mydat1)

# run MR 
res1 <- mr(hdat1, method = "mr_wald_ratio")
generate_odds_ratios(res1)
# get sigs
sig1 <- res1 %>% filter(pval < 0.05) # 6 are significant; check similarity of outcomes
finn.r2 <- get_r_from_lor(lor = 0.97, af = 0.276, ncase = 23363, ncontrol = 195429, prevalence = 23363/(23363+195429))
finn.r2e <- get_r_from_bsen(b=1.115, se = 0.21, n=125)
hdat1$r.outcome = finn.r2
hdat1$r.exposure = finn.r2e
directionality_test(hdat1 %>% filter(id.outcome == "finn-b-I9_CORATHER_EXNONE"))
mr_steiger(p_exp = 5.21e-7,
           p_out = 6.0e-6,
           n_exp = 125,
           n_out = 122733+424528,
           r_exp = finn.r2e,
           r_out = finn.r2)
