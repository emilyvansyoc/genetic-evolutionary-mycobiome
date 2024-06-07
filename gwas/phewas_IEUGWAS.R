### PheWAS using IEU OpenGWAS
# EVS 5/1/2024

library(tidyverse)
library(ieugwasr)
library("TwoSampleMR")

### website to get authentication token (login through Github):
# https://api.opengwas.io/profile/

# janky; has a hard time reading the Renviron file
mytoken <- "/mytokencharactersfromAPIthingy"


# get SNP list
myids <- read.table("data/unique_mycAVS_SNPandStructural.txt", sep = "\t", header = TRUE) 
snplist <- myids$SNPID


# test validation of heart snps in other coronary disease cohorts
phewas(variants = c("rs12149890", "rs12929586"), batch = "bbj-a-159", opengwas_jwt = mytoken)
phewas(variants = c("rs12149890", "rs12929586"), batch = "ebi-a-GCST003116", opengwas_jwt = mytoken)

# get associations from a particular study (fast)
associations(variants = c("rs12149890", "rs12929586"),
             id = c("bbj-a-159", "ebi-a-GCST003116"))

# what does 1000 Genomes say about these variants?
ieugwasr::afl2_rsid(rsid = c("rs12149890", "rs12929586"), opengwas_jwt = mytoken)


