## MWAS in a linear model
# EVS 4/2024

library(tidyverse)
library(rstatix)
library(phyloseq)
library(microViz)
library(vegan)
library(foreach)
library(parallel)
library(doParallel)

## get fungal data 
load("data/updated/phylo-ITS-filteredrarefied_nosubset.RData")

# get pathways
load("wgs-analyses/wgs-data/stool_pathways.RData")


# for now; get just community level
# remove pathways with zero abundnace
# and require 10% prevalence
pathpa <- pathf %>% 
  filter(!str_detect(Pathway, "\\|")) %>% 
  mutate(pa = if_else(avgAbund > 0, 1, 0)) %>% 
  group_by(Pathway) %>% 
  summarize(tot = sum(pa)) %>% 
  filter(tot > (83*0.1))
pathfilt <- pathf %>% 
  filter(!str_detect(Pathway, "\\|")) %>% 
  filter(Pathway %in% pathpa$Pathway)
length(unique(pathfilt$Pathway)) # 248

# filter individuals to match shotgun data
psfilt <- psf %>% 
  ps_filter(sample_names(psf) %in% pathfilt$RANDSID)

# make fungi vertical and get genus level abundance data
psclr <- psfilt %>% 
  
  tax_fix() %>% 
  tax_glom("Genus") %>% 
  # filter for 10% prevalence
  tax_filter(min_prevalence = (83*0.1)) %>% # 26 genera
  tax_transform("clr", zero_replace = "halfmin")
taxa_names(psclr) <- psclr@tax_table[,"Genus"]
fungv <- psclr %>% 
  ps_melt() %>% 
  dplyr::select(Fungi = OTU, Sample, Abundance) %>% 
  mutate(RANDSID = as.numeric(Sample)) # 26 fungal genera

# join
moddf1 <- fungv %>% 
  full_join(pathfilt, relationship = "many-to-many")
# remove pathways with zero total abundance
compa <- moddf1 %>% 
  mutate(pa = if_else(avgAbund > 0, 1, 0)) %>% 
  group_by(Pathway) %>% 
  summarize(tot = sum(pa)) %>% 
  filter(tot > 0) %>% ungroup() # 361
moddf1 <- moddf1 %>% 
  filter(Pathway %in% compa$Pathway)

### build loop
fungs <- unique(moddf1$Fungi)
npath <- unique(moddf1$Pathway)
outstat1 <- data.frame()

for(i in 1:length(npath)) { # open i loop
  
  # print progress every 10th pathway
  if(i %% 10==0) {cat("\n working on pathway ", i, "of ", length(npath))}
  
  for(j in 1:length(fungs)) { # open j loop
    
    # model
    mod <- lm(Abundance ~ log10(avgAbund + 0.1), data = moddf1 %>% filter(Fungi == fungs[j] & Pathway == npath[i]))
    
    # get output
    outstat1 <- rbind(outstat1, data.frame(
      pathway = npath[i],
      fungi = fungs[j],
      estimate = broom::tidy(summary(mod))$estimate[2],
      stderror = broom::tidy(summary(mod))$std.error[2],
      stat = broom::tidy(summary(mod))$statistic[2],
      pval = broom::tidy(summary(mod))$p.value[2],
      rsquared = summary(mod)$r.squared,
      fstat = unname(summary(mod)$fstatistic[1]),
      resstderror = summary(mod)$sigma
      
    ))
    
  } # close j loop
  
  
} 

## ---- results ----

res <- outstat1


#### remove pathways that are explictly found in fungi and/or yeasts
resdf <- res %>% 
  filter(!str_detect(pathway, "yeast")) %>% 
  filter(!str_detect(pathway, "fung")) %>% 
  filter(!str_detect(pathway, "eukaryo")) # only removed 10 pathways 

# adjust p values and get sigs
resdf <- resdf %>% 
  mutate(pfdr = p.adjust(pval, method = "fdr"),
         pbon = p.adjust(pval, method = "bonferroni"))
sigs <- resdf %>% filter(pfdr < 0.05)

