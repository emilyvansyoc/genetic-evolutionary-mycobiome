### GTEx analysis
# EVS 1/2024, updated to use the query database 5/17/2024

library(tidyverse)

# get FAVs
ids <- read.table("data/unique_mycAVS_SNPandStructural.txt", header = TRUE)
sigs <- dat %>% filter(ID %in% ids$SNPID) %>% mutate(chr = paste0("chr", X.CHROM))

## ----- convert rsIDs into b38 genomic coordinates -----

## get GTEx lookup table (too big to upload to Github)
# available from GTEx portal 
gtex <- read.table("private/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt", sep = "\t", header = TRUE)

# decrease size 
g1 <- gtex %>% filter(chr %in% sigs$chr)

# format sigsnps for matching
sigs <- sigs %>% 
  mutate(variant_id_b37 = paste(X.CHROM, POS, REF, ALT, "b37", sep = "_")) %>% 
  mutate(pos.flip = paste(X.CHROM, POS, ALT, REF, "b37", sep = "_"))

# match by SNPID
m1 <- g1 %>% 
  inner_join(sigs %>% dplyr::select(ID, pos.flip, myvariant_id_b37 = variant_id_b37), by = c("rs_id_dbSNP151_GRCh38p7" = "ID"))

nom1 <- sigs %>% dplyr::select(ID, variant_id_b37) %>% anti_join(g1, by = c("ID" = "rs_id_dbSNP151_GRCh38p7")) %>% left_join(dat)
# these don't exist in GTEx?
nom2 <- g1 %>% filter(variant_id_b37 %in% nom1$variant_id_b37) # not there

## from the matches, verify that the ref/alt alleles are the same and get b38 IDs
m1 <- m1 %>% 
  mutate(matches = case_when(
    variant_id_b37 == myvariant_id_b37 ~ "perfect",
    variant_id_b37 == pos.flip ~ "flipped",
    variant_id_b37 != myvariant_id_b37 & variant_id_b37 != pos.flip ~ "nomatch"
  )) # all are either a perfect match or have flipped alleles

#### ---- get GTEx data ----

## search for presence of these SNPs in the GTEx data
# downloaded from: https://gtexportal.org/home/downloads/adult-gtex/qtl

## get list of significant eQTL-variant pairs
outdf1 <- data.frame()
myfi1 <- list.files(mypath, pattern = "*signif*", full.names = TRUE)
for(i in 1:length(myfi1)) {
  cat("working on", i, "of", length(myfi1), "\n")
  fi <- read.table(myfi1[i], header = TRUE, sep = "\t")
  fi$tissue <- str_extract(myfi1[i], "[:alpha:]*\\.v8\\.signif")
  outdf1 <- rbind(outdf1, fi)
}


#### ---- get eQTLs for significant snps ----

eqtl <- outdf1 %>% filter(variant_id %in% m1$variant_id)
# join to our data
eqtl <- eqtl %>% full_join(m1) #### THIS MATCHES THE WEB BROWSER RESULTS

length(unique(eqtl$variant_id[!is.na(eqtl$gene_id)])) # 11 unique variants

# save
save(eqtl, file = "data/GTEx_sigSNPs_fromLUtable.RData")
