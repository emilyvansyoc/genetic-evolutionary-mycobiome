### separate 'all variants' into SNPs and structural
# EVS 2/2024 with updated ITS methods
# 3/2024 with snp v structural

library(tidyverse)
library(gridExtra)
library(ggpubr)
#library(qqman)
library(phyloseq)
library(microViz)
library(cowplot)

# set pathway to output files (straight from PLINK)
path <- "gwas-out/updateITS_results_allsnp_allvariant/"


### ---- get effect sizes ----


# read everything in to one dataframe (small enough to do this!)
myfiles1 <- list.files(path, pattern = ".linear$")
alleffects <- data.frame()
for(i in 1:length(myfiles1)) {
  cat("\n now reading file ", i, " of ", length(myfiles1))
  
  fi <- read.table(paste0(path, myfiles1[i]), sep = "\t", header = TRUE, comment.char = "") 
  fi$taxa <- myfiles1[i]
  alleffects <- rbind(alleffects, fi)
}

## get SNPs and significants
snp <- alleffects %>% 
  filter(nchar(REF) == 1 & nchar(ALT) == 1) # 4.6M; matches PLINK output exactly
# adjust for genome-wide significant
snp <- snp %>% 
  group_by(taxa) %>% 
  mutate(gensig = p.adjust(P, method = "fdr"))
snpsig <- snp %>% filter(gensig < 0.05)

## get structural and significant
struc <- alleffects %>% 
  filter(nchar(REF) > 1 | nchar(ALT) > 1)
# adjust for genome-wide significance
struc <- struc %>% 
  group_by(taxa) %>% 
  mutate(gensig = p.adjust(P, method = "fdr"))
strucsig <- struc %>% filter(gensig < 0.05)

# get summary stats
unique(strucsig$taxa)
unique(snpsig$taxa)
length(unique(strucsig$ID)) # 9 
length(unique(snpsig$ID)) # 135

# get effect sizes - SNP
summary(abs(snpsig$BETA)) # highest is 2.8
min(snpsig$gensig)
snpsig %>% filter(BETA > 2.8) %>% dplyr::select(taxa, gensig) # biggest effect size and smallest P value is Pleosporales
snpsig %>% group_by(taxa) %>% count() %>% arrange(desc(n))
snpsig %>% group_by(taxa) %>% 
  summarize(avgef = mean(BETA)) 
# highest is Pleosporales; capnodiales; then Saccharomycetes class and order; then Kazachstania

# get effect sizes - structural
summary(abs(strucsig$BETA)) # highest is 2.95
min(strucsig$gensig) 
strucsig %>% filter(BETA > 2.8) %>% dplyr::select(taxa, gensig) # both Pleosporales
strucsig %>% group_by(taxa) %>% count()
strucsig %>% group_by(taxa) %>% summarize(avgef = mean(BETA)) # highest is Pleosporales, then Capnodiales then Kazachstania

# get unique chromosomes
unique(snpsig$X.CHROM)
unique(strucsig$X.CHROM)

## add together and get summary stats
allsig <- snpsig %>% mutate(Type = "SNP") %>% 
  rbind(strucsig %>% mutate(Type = "Structural")) %>% 
  dplyr::select(taxa, Type, ID, X.CHROM, POS, BETA, A1_FREQ, P, gensig) %>% 
  mutate(taxa = str_remove(taxa, "allsnp_allvariants\\."),
         taxa = str_remove(taxa, "\\.glm\\.linear")) %>% 
  mutate(posID = paste(X.CHROM, POS, sep = ":")) %>% 
  dplyr::select(-c(X.CHROM, POS)) %>% 
  relocate(taxa, Type, ID, posID) %>% 
  arrange(gensig) %>% ungroup()

# get final counts
length(unique(allsig$ID))
length(unique(allsig$ID[allsig$Type == "SNP"])) # 135
length(unique(allsig$ID[allsig$Type == "Structural"])) # 9


names(allsig) <- c("Taxa", "Variant Type", "SNPID", "Chromosome:Position", "Beta", "MAF", "Raw P", "Genome-Wide Q")
write.table(allsig, file = "R/gwas-output/sigtable_SNPandStructural.txt", sep = "\t", row.names = FALSE, quote = FALSE)
save(allsig, file = "R/gwas-output/sigs_SNPandStructural.RData")

# get updated mycavs
myc <- allsig %>% ungroup() %>% dplyr::select(SNPID) %>% distinct()

## ---- format for downstream analysis: SNPNEXUS ----

# write out
write.table(myc, file = "R/gwas-output//unique_mycAVS_SNPandStructural.txt", sep = "\t", row.names = FALSE)

# format for SNPnexus
nex <- myc %>% mutate(Type = "dbsnp") %>% rename(Name = SNPID) %>% relocate(Type)
write.table(nex, file = "R/gwas-output/unique_mycAVs_forSNPNexus_SNPandStructural.txt", sep = "\t", row.names = FALSE, quote = FALSE)





