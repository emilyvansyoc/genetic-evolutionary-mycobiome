#### compare fungal abundance to eQTL allele frequency and gene expression, and generate figure S4
# EVS 4/2024

library(tidyverse)
library(rstatix)
library(phyloseq)
library(microViz)
library(ggpubr)
library(vegan)
library(ggtext)
library(phyloseq)
library(microViz)

## ---- get fungal data ----

### get fungal data
load("data/phylo_ITS_resolvedNA.RData")
psf <- psnewname
# collapse and process via GWAS 
# get function
source("R/fx_myCollapse.R")

# collapse all ranks and remove redundant taxa (e.g.; a Genus with only one OTU in it is redundant with the OTU and should be removed)
allCol <- myCollapse(phylo = psnewname, tax_levels = c("OTU", "Genus", "Family", "Order", "Class", "Phylum"))


# remove any collapse taxa not present in at least 30% of individuals
pa <- decostand(allCol %>% 
                  column_to_rownames(var = "Sample"), method = "pa")  
csum <- colSums(pa)
tokeep <- names(which(csum > (nrow(allCol) * 0.30)))
allColfilt <- allCol %>% 
  select(Sample, all_of(tokeep)) %>% 
  column_to_rownames(var = "Sample")
# re-phyloseq for transformation
rephy <- phyloseq(
  otu_table(allColfilt, taxa_are_rows = FALSE),
  tax_table(data.frame(row.names = names(allColfilt), taxaname = names(allColfilt)) %>% as.matrix())
) %>% 
  tax_transform("clr", zero_replace = "halfmin")
# get data frame
allclr <- rephy %>% otu_get() %>% data.frame()

# get Kazachstania
kz <- allclr %>% dplyr::select(Genus_Kazachstania) %>% rownames_to_column(var = "IID") %>% mutate(IID = as.numeric(IID))


### ---- get gtes ----

# filter to the variants of interest and get SNP IDs
myids <- read.table("data/unique_mycAVS_SNPandStructural.txt", sep = "\t", header = TRUE) 

mygte <- read.table("data/gtex_results.txt", header = TRUE, sep = "\t") %>% 
  rename(Gencode = `Gencode.Id`,
         Gene = `Gene.Symbol`,
         Variant = `Variant.Id`,
         SNPID = `SNP.Id`,
         Pval = `P.Value`) # 67
# filter
mygte <- mygte %>% 
  filter(SNPID %in% myids$SNPID)


## ---- get ref and alt alleles ----

load("data/GTEx_sigSNPs_fromLUtable.RData")
# get sigs
eqtl <- eqtl %>% drop_na()
# one ref and alt alleles match ours
eqtl %>% filter(matches == "perfect")
# others are flipped
eqtl %>% filter(matches == "flipped")
### here, "slope" is equivalent to the GTEx portal "NES" effect size 
# so a positive slope is higher expression in the alt allele
# and a negative slope is higher expression in the ref allele

## add to our GTE results
mygte1 <- mygte %>% 
  full_join(eqtl %>% dplyr::select(Variant = variant_id, ref, alt) %>% distinct())

## ---- get genotypes ----

## these are PLINK output files 
#### RESTRICTED
ped <- read.table("gwas-out/eqtl_genos.ped", header = FALSE)
map <- read.table("gwas-out/eqtl_genos.map", header = FALSE)
names(ped) <- c("IID", "FID", "fam1", "fam2", "fam3", "fam4", "rs10020593_a1", "rs10020593_a2", "rs6842499_a1", "rs6842499_a2", "rs2056279_a1", "rs2056279_a2", "rs2056278_a1", "rs2056278_a2", "rs10003727_a1", "rs10003727_a2", "rs7695409_a1", "rs7695409_a2", "rs11734845_a1", "rs11734845_a2", "rs2162092_a1", "rs2162092_a2", "rs1559699_a1", "rs1559699_a2", "rs28673203_a1", "rs28673203_a2", "rs192598_a1", "rs192598_a2")

## collapse eQTLs for grouped genes
colgte <- mygte %>% dplyr::select(Gene, SNPID, NES) %>% 
  pivot_wider(names_from = Gene, values_from = NES) 
mylabs <- data.frame(SNPID = colgte$SNPID,
                     SNPlab = c("rs192598 (*CDH13*)", "rs10020593 (*OTUD4*, *HHIP*)", "rs6842499 (*OTUD4*, *HHIP*)",
                                "rs2056279 (*OTUD4*, *HHIP*, *HHIP-AS1*)", "rs2056278 (*OTUD4*, *HHIP*, *HHIP-AS1*)",
                                "rs10003727 (*OTUD4*, *HHIP*)", "rs7695409 (*OTUD4*, *HHIP*)", "rs11734845 (*OTUD4*, *HHIP*, *HHIP-AS1*)", "rs2162092 (*OTUD4*, *HHIP*, *HHIP-AS1*)", "rs1559699 (*OTUD4*, *HHIP*, *HHIP-AS1*)", "rs28673203 (*OTUD4*, *HHIP*, *HHIP-AS1*)"),
                     genelab = c("*CDH13* (dec.)", "*OTUD4* and *HHIP* (both inc.)", "*OTUD4* and *HHIP* (both inc.)", "*OTUD4* (inc.), *HHIP* (inc.), and *HHIP-AS1* (dec.)", "*OTUD4* (inc.), *HHIP* (inc.), and *HHIP-AS1* (dec.)", "*OTUD4* and *HHIP* (both inc.)", "*OTUD4* and *HHIP* (both inc.)", "*OTUD4* (inc.), *HHIP* (inc.), and *HHIP-AS1* (dec.)", "*OTUD4* (inc.), *HHIP* (inc.), and *HHIP-AS1* (dec.)", "*OTUD4* (inc.), *HHIP* (inc.), and *HHIP-AS1* (dec.)", "*OTUD4* (inc.), *HHIP* (inc.), and *HHIP-AS1* (dec.)"))

# make full genotypes
gens <- ped %>% 
  dplyr::select(-c(starts_with("fam"), "FID")) %>% 
  mutate(rs10020593_geno = paste0(rs10020593_a1, rs10020593_a2),
         rs6842499_geno = paste0(rs6842499_a1, rs6842499_a2),
         rs2056279_geno = paste0(rs2056279_a1, rs2056279_a2),
         rs2056278_geno = paste0(rs2056278_a1, rs2056278_a2),
         rs10003727_geno = paste0(rs10003727_a1, rs10003727_a2),
         rs7695409_geno = paste0(rs7695409_a1, rs7695409_a2),
         rs11734845_geno = paste0(rs11734845_a1, rs11734845_a2),
         rs2162092_geno = paste0(rs2162092_a1, rs2162092_a2),
         rs1559699_geno = paste0(rs1559699_a1, rs1559699_a2),
         rs28673203_geno = paste0(rs28673203_a1, rs28673203_a2),
         rs192598_geno = paste0(rs192598_a1, rs192598_a2))

# make vertical
genv <- gens %>% 
  dplyr::select(IID, ends_with("geno")) %>% 
  pivot_longer(cols = !IID, names_to = "ID", values_to = "genotype") %>% 
  # add kazachstania
  right_join(kz) %>% 
  # get genotype
  mutate(geno.letters = case_when(
    ID %in% "rs10020593_geno" & genotype %in% "TT" ~ "homo.ref",
    ID %in% "rs10020593_geno" & genotype %in% "TA" ~ "het",
    ID %in% "rs10020593_geno" & genotype %in% "AA" ~ "homo.alt",
    ID %in% "rs6842499_geno" & genotype %in% "AA" ~ "homo.ref",
    ID %in% "rs6842499_geno" & genotype %in% "AT" ~ "het",
    ID %in% "rs6842499_geno" & genotype %in% "TT" ~ "homo.alt",
    ID %in% "rs2056279_geno" & genotype %in% "CC" ~ "homo.ref",
    ID %in% "rs2056279_geno" & genotype == "CT" ~ "het",
    ID == "rs2056279_geno" & genotype == "TT" ~ "homo.alt",
    ID == "rs2056278_geno" & genotype == "CC" ~ "homo.ref",
    ID == "rs2056278_geno" & genotype == "CG" ~ "het",
    ID == "rs2056278_geno" & genotype == "GG" ~ "homo.alt",
    ID == "rs10003727_geno" & genotype == "CC" ~ "homo.ref",
    ID == "rs10003727_geno" & genotype == "CT" ~ "het",
    ID == "rs10003727_geno" & genotype == "TT" ~ "homo.alt",
    ID == "rs7695409_geno" & genotype == "TT" ~ "homo.ref",
    ID == "rs7695409_geno" & genotype == "TC" ~ "het",
    ID == "rs7695409_geno" & genotype == "CC" ~ "homo.alt",
    ID == "rs11734845_geno" & genotype == "GG" ~ "homo.ref",
    ID == "rs11734845_geno" & genotype == "GA" ~ "het",
    ID == "rs11734845_geno" & genotype == "AA" ~ "homo.alt",
    ID == "rs2162092_geno" & genotype == "GG" ~ "homo.ref",
    ID == "rs2162092_geno" & genotype == "GA" ~ "het",
    ID == "rs2162092_geno" & genotype == "AA" ~ "homo.alt",
    ID == "rs1559699_geno" & genotype == "GG" ~ "homo.ref",
    ID == "rs1559699_geno" & genotype == "GC" ~ "het",
    ID == "rs1559699_geno" & genotype == "CC" ~ "homo.alt",
    ID == "rs28673203_geno" & genotype == "AA" ~ "homo.ref",
    ID == "rs28673203_geno" & genotype == "AG" ~ "het",
    ID == "rs28673203_geno" & genotype == "GG" ~ "homo.alt",
    ID == "rs192598_geno" & genotype == "AA" ~ "homo.ref",
    ID == "rs192598_geno" & genotype == "AC" ~ "het",
    ID == "rs192598_geno" & genotype == "CC" ~ "homo.alt"
  )) %>% 
  ## pretty-fy
  mutate(geno.letters = case_when(
    geno.letters == "homo.ref" ~ "REF",
    geno.letters == "het" ~ "Heterozy.",
    geno.letters == "homo.alt"~ "ALT"
  )) %>% 
  mutate(geno.letters = factor(geno.letters, ordered = TRUE, levels = c("REF", "Heterozy.", "ALT"))) %>% 
  # get GTE effect size
  mutate(SNPID = str_remove(ID, "_geno")) %>% 
  left_join(mygte %>% dplyr::select(Gene, SNPID, NES), relationship = "many-to-many") %>% 
  # get directionality of gene expression
  mutate(genexp.dir = if_else(NES > 0, "Inc", "dec")) %>% 
  mutate(label = paste(geno.letters, genotype, sep = "\n")) %>% 
  left_join(mylabs) %>% 
  mutate(genelab = factor(genelab, ordered = TRUE, levels = c(
    "*CDH13* (dec.)", "*OTUD4* and *HHIP* (both inc.)", "*OTUD4* (inc.), *HHIP* (inc.), and *HHIP-AS1* (dec.)"
  )))  %>% 
  mutate(SNPID = factor(SNPID, ordered = TRUE, levels = c("rs192598", "rs10020593", "rs6842499",
                                                          "rs10003727", "rs7695409", "rs2056279", "rs2056278", "rs11734845", "rs2162092", "rs1559699", "rs28673203")))

### plot opposite way with only one representative SNP of the cluster
forplot <- genv %>% filter(SNPID %in% c("rs192598", "rs10020593", "rs2056279")) %>% 
  mutate(SNPkey = case_when(
    genelab %in% "*CDH13* (dec.)" ~ "rs192598",
    genelab %in% "*OTUD4* and *HHIP* (both inc.)" ~ "rs10020593, rs6842499, rs10003727, <br> rs7695409",
    genelab %in% "*OTUD4* (inc.), *HHIP* (inc.), and *HHIP-AS1* (dec.)" ~ "rs2056279, rs2056278, rs1173845, <br> rs2162092, rs1559699, rs28673203"
  )) %>% 
  mutate(SNPkey = factor(SNPkey, ordered = TRUE, levels = c(
    "rs192598", "rs10020593, rs6842499, rs10003727, <br> rs7695409",
    "rs2056279, rs2056278, rs1173845, <br> rs2162092, rs1559699, rs28673203"
  )))
statplot <- forplot %>% group_by(genelab) %>% 
  dunn_test(Genus_Kazachstania ~ geno.letters) %>% 
  add_significance() %>% 
  add_xy_position() %>% 
  mutate(ast = if_else(p.adj < 0.05, "*", "ns")) %>% 
  mutate(y.position = y.position + 0.5)
# plot
ggboxplot(forplot,
          x = "geno.letters", y = "Genus_Kazachstania",
          facet.by = "genelab", fill = "SNPkey",
          size = 1.5,
          xlab = "Genotype", ylab = "*Kazachstania* abundance") +
  ylim(c(-4,9)) +
  stat_pvalue_manual(statplot, label = "ast", hide.ns = TRUE, size = 8, bracket.size = 0.8, tip.length = 0) +
  theme_pubr(base_size = 14) +
  guides(fill = guide_legend(title = "FAV-eQTLs")) +
  theme(
    axis.title.y = element_markdown(),
    strip.text = element_markdown(),
    legend.text = element_markdown(),
    legend.position = "bottom"
  )
ggsave(filename = "figures/eqtl_kaz_v1.png", dpi = 600)
