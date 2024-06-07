### parafit and PACO results
# EVS 4/2024

library(ape)
library(vegan)
library(tidyverse)
library(phyloseq)
library(microViz)
library(Biostrings)
library(purrr)


## get trees
load("data/all_otu_subtrees_04-10-2024_nonegbl.RData")

# get Parafit results
load("data/parafit_results_OTUs_nonegbl_04-10-2024.RData")
# get sigs
paradf <- paradf %>% 
  mutate(pfdr = p.adjust(pval, method = "fdr"),
         pbon = p.adjust(pval, method = "bonferroni"))

## get PACo output - these are in individual text files from running on ROAR
dir <- "path/to/codiv_tests_final/PACo_OTU_quasi_p100_nonegbl_04-11-2024/"
myfi <- list.files(dir, pattern = ".txt", full.names = TRUE)
pacodf <- data.frame()
for(i in 1:length(myfi)) {
  fi <- read.table(myfi[i], sep = "\t", header = TRUE)
  pacodf <- rbind(pacodf, fi)
}
pacodf <- pacodf %>% 
  mutate(padj = p.adjust(paco.p, method = "fdr"),
         pbon = p.adjust(paco.p, method = "bonferroni"))

## join
allotu <- paradf %>%
  dplyr::select(otu, para.p = pval, parastat, para.pfdr = pfdr, para.pbon = pbon) %>% 
  full_join(pacodf %>% dplyr::select(otu = subtree, paco.ss, paco.p, paco.pfdr = padj, paco.pbon = pbon)) %>% 
  # calculate paco R2
  mutate(pacor2 = 1-paco.ss)
summary(allotu$pacor2)
## get sigs from both
bsigs <- allotu %>% filter(para.pfdr < 0.05 & paco.pfdr < 0.05)

summary(bsigs$pacor2)
summary(bsigs$parastat) # not as easily interpretable as the r2
## 
# get biggest effect size
bsigs %>% filter(pacor2 > 0.25) # only OTU_144 
quantile(bsigs$pacor2)

# save for supplementary materials
tax <- psf %>% tax_fix() %>% tax_select(allotu$otu, strict_matches = TRUE, n_typos = 0) %>% 
  tt_get() %>% as.data.frame() %>% 
  rownames_to_column(var = "otu") %>% 
  full_join(allotu) %>% 
  dplyr::select(otu, Phylum, Class, Order, Family, Genus, Species,
                parastat, para.pfdr, pacor2, paco.pfdr) %>%
  arrange(desc(pacor2))
write.table(tax, file = "data/allcodiv_results.txt", sep = "\t", row.names = FALSE)


### ----- get fungal information -----

# get phyloseq
load("private/hominid_phyloITS.RData")

## get the rest of the OTUs in bsigs
sigotus <- bsigs$otu
sigps <- psf %>% tax_fix() %>% tax_select(sigotus, strict_matches = TRUE, n_typos = 0) %>% 
  tt_get() %>% as.data.frame() %>% 
  rownames_to_column(var = "otu") %>% full_join(bsigs)
# save
save(sigps, file = "data/final_sigs_otus_parafitpaco.RData")
