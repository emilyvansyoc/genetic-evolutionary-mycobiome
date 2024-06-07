### summary statistics and linear model for Kazachstania
####NOTE: heavily retracted since the demographic/clinical data is restricted on dbGaP

# EVS 3/2024

library(microViz)
library(phyloseq)
library(tidyverse)
library(ggpubr)
library(rstatix)

# get the input abundance values for GWAS
abun <- read.table("private/pheno_ITS.txt", sep = "\t", header = TRUE)

## get filtered and rarefied phyloseq
load("data/phylo_ITS_resolvedNA.RData")
psf <- psnewname

# get Kazachstania
ka <- psf %>% tax_select("Kazachstania", "Genus")  # 6 taxa
# these are in the family Saccharomycetaceae; of 6 taxa, 4 have species 
otus <- paste0("OTU_", taxa_names(ka))

# were any OTUs tested in the GWAS?
abun %>% dplyr::select(any_of(otus)) %>% head() # only OTU 266 had high enough prevalence; unknown species

### investigate its prevalence and abundance
psg <- psf %>% tax_fix() %>% tax_glom("Genus")
taxa_names(psg) <- psg@tax_table[,"Genus"]
tax_top(psg) # not in the top 10
tax_top(psg, n = 20) # 14th ranked in abundance
# get relative abundance
kar <- psf %>% tax_glom("Genus") %>% tax_transform("compositional") %>% tax_select("Kazachstania", "Genus") %>% 
  ps_melt()
kar %>% get_summary_stats(Abundance)
# remove zeros
kar %>% filter(Abundance > 0) %>% get_summary_stats(Abundance)

kaz <- psg %>% ps_melt() %>% 
  filter(Genus == "Kazachstania")
hist(kaz$Abundance) # most in very low abundance
summary(kaz$Abundance) # one person has very high abundance 
summary(kaz$Abundance[kaz$Abundance < 1000])
length(which(kaz$Abundance == 0)) # absent in 32 indivduals of 125 ~ 25%
length(which(kaz$Abundance > 0))

# get CLR abundances
kc <- psg %>% tax_transform("clr", zero_replace = "halfmin") %>% ps_melt() %>% filter(Genus == "Kazachstania")
kaz <- kaz %>% 
  left_join(kc %>% dplyr::select(OTU, RANDSID, CLR = Abundance))


### ----- kazachstania and demographic/clinical vars ----

## this data and wrangling is retracted; restricted use on dbGaP ##

## test: age, gender, BMI, tobacco use, ethnicity, blood pressure
mod <- lm(kaz_CLR ~ GENDER + AGEENR + sr.race.num + DTPSYSTL + DTPBMI + DVDTOBC + DHXCARD + DHXGI, data = df)
car::Anova(mod, type = "II")
#plot(mod)
summary(mod) # AGE is the only significance; gender is close 

# plot
ggscatter(df, x = "AGEENR", y = "kaz_CLR", add = "reg.line")
