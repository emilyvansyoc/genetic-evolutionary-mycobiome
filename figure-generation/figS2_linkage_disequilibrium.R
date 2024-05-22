### SNP linkage disequilibrium and generate Figure S2

#EVS 5/1/2024

library(tidyverse)
library(ggpubr)
library(ggtext)

### PLINK code to generate linkage disequilibrium report 
#$PLINK --bfile $BFILE --r2 --ld-snp-list unique_mycAVs.txt --out ld_report_snps --ld_window_r2 0

## ---- ld report ----
# get LD report
ld <- read.table("gwas-out/ld_report_snps_all.ld", sep = "\t", header = TRUE)

# remove duplicates (fixed in plink2 but hasn't been implemented)
ld1 <- ld %>% filter(!SNP_A == SNP_B)

# get only the variants in our mycAVs (not sure why extras are added?)
myids <- read.table("R/gwas-output//unique_mycAVS_SNPandStructural.txt", sep = "\t", header = TRUE) 
ld2 <- ld1 %>% 
  filter(SNP_A %in% myids$SNPID & SNP_B %in% myids$SNPID)

# make horizontal
ldh <- ld2 %>% 
  dplyr::select(SNP_A, SNP_B, R2) %>% 
  pivot_wider(names_from = SNP_B, values_from = R2)

## ---- get ld for each SNP ----
# get all significant FAVs with their taxa
load("R/gwas-output/sigs_SNPandStructural.RData")
allsig <- allsig %>% 
  mutate(chr = sapply(str_split(`Chromosome:Position`, ":"), `[`, 1))

# get pleo correlations
pleo <- allsig %>% filter(Taxa == "Order_Pleosporales") %>% 
  dplyr::select(Taxa, SNPID, VType = `Variant Type`, chr) %>% 
  left_join(ldh, by = c("SNPID" = "SNP_A")) %>% 
  pivot_longer(cols = starts_with("rs"), names_to = "SNP_B", values_to = "R2") %>% 
  drop_na(R2)

# get Saccharomycetales
sacc <- allsig %>% 
  filter(Taxa == "Order_Saccharomycetales") %>% 
  dplyr::select(Taxa, SNPID, VType = `Variant Type`, chr) %>% 
  left_join(ldh, by = c("SNPID" = "SNP_A")) %>% 
  pivot_longer(cols = starts_with("rs"), names_to = "SNP_B", values_to = "R2") %>% 
  drop_na(R2)

# get Capnodiales
cap <- allsig %>% 
  filter(Taxa == "Order_Capnodiales") %>% 
  dplyr::select(Taxa, SNPID, VType = `Variant Type`, chr) %>% 
  left_join(ldh, by = c("SNPID" = "SNP_A")) %>% 
  pivot_longer(cols = starts_with("rs"), names_to = "SNP_B", values_to = "R2") %>% 
  drop_na(R2)

# get Kazachstania
kaz <- allsig %>% 
  filter(Taxa == "Genus_Kazachstania") %>% 
  dplyr::select(Taxa, SNPID, VType = `Variant Type`, chr) %>% 
  left_join(ldh, by = c("SNPID" = "SNP_A")) %>% 
  pivot_longer(cols = starts_with("rs"), names_to = "SNP_B", values_to = "R2") %>% 
  drop_na(R2)

### add together and plot
all <- pleo %>% rbind(sacc) %>% rbind(cap) %>% rbind(kaz) %>% 
  mutate(Taxa = str_replace(Taxa, "_", " ")) %>% 
  mutate(Taxa = if_else(str_detect(Taxa, "Kaz"), "Genus *Kazachstania*", Taxa)) %>% 
  mutate(Taxa = factor(Taxa, ordered = TRUE, levels = c(
    "Order Pleosporales", "Order Saccharomycetales", "Order Capnodiales", "Genus *Kazachstania*"
  ))) %>% 
  mutate(chr = factor(chr, ordered = TRUE, levels = c(1, 4, 5, 10, 16, 18)))

# make plot
gghistogram(all, x = "R2", fill = "chr", facet.by = "Taxa", scales = "free",
            bins = 10, xlab = "R<sup>2</sup>", ylab = "Number of FAVs",
            position = "stack") +
  guides(fill = guide_legend(title = "Chromosome")) +
  theme_pubr(base_size = 14, legend = "right") +
  theme(strip.text = element_markdown(size = 18),
        axis.title.x = element_markdown(),
        strip.background = element_rect(fill = "white"))
ggsave(filename = "figures/ld_report_FAVs.png", dpi = 600)