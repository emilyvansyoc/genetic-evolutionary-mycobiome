### generate Figure 1 
# this is a Frankenscript


library(tidyverse)
library(ggpubr)
library(scales)
library(cowplot)
library(plotly)
library(ggtext) # for element_markdown
library(rcartocolor) # for colorbliend safe palette
library(rstatix)


### ---- generate manhattan plots ----

## get GWAS data for Manhattan plots
load(file = "data/sigtaxa_alleffects.RData")
sigtaxa <- unique(allsig$Taxa)
load(file = "data/sigs_SNPandStructural.RData")


## for each taxa, get Manhattan plot coordinates for chromosome names
plotdf <- data.frame()
axisdf <- data.frame()
for(i in 1:length(sigtaxa)) {
  
  # wrangle
  sub <- df %>% 
    filter(taxa %in% sigtaxa[i]) 
  
  # get coordinates
  coords <- sub %>% 
    # get coordinates
    group_by(X.CHROM) %>% 
    summarize(chr_len = max(POS)) %>% 
    mutate(tot = cumsum(as.numeric(chr_len)) - chr_len) %>% 
    dplyr::select(-chr_len) %>% 
    left_join(sub) %>% 
    arrange(X.CHROM, POS) %>% 
    mutate(BPcum = POS + tot)
  
  # join back
  plotdf <- rbind(plotdf, coords)
  
  # get positions of X axis titles
  ax <- coords %>% group_by(X.CHROM) %>% summarize(center = (max(BPcum) + min(BPcum)) / 2) %>% 
    mutate(taxa = sigtaxa[i])
  axisdf <- rbind(axisdf, ax)
  
}


# reduce size
manplotdf <- plotdf %>% filter(-log10(P) > 1.7) %>% 
  # keep only necessary columns
  dplyr::select(X.CHROM, POS, P, taxa, BPcum) %>% 
  
  ### 4/8/2024; remove Sacc class
  filter(!taxa == "Class_Saccharomycetes") %>% 
  # pretty-fy names
  mutate(plotnames = case_when(
    taxa %in% "Order_Pleosporales" ~ "**A.** Order Pleosporales",
    taxa %in% "Order_Saccharomycetales" ~ "**B.** Order Saccharomycetales",
    #taxa %in% "Class_Saccharomycetes" ~ "**C.** Order Saccharomycetes",
    taxa %in% "Order_Capnodiales" ~ "**C.** Order Capnodiales",
    taxa %in% "Genus_Kazachstania" ~ "**D.** Genus *Kazachstania*"
  )) %>% 
  # order
  mutate(plotnames = factor(plotnames, ordered = TRUE, levels = c("**A.** Order Pleosporales", "**B.** Order Saccharomycetales", "**C.** Order Capnodiales", "**D.** Genus *Kazachstania*")))

# show every chromosome label until 10, then skip every other until 22
axisdf1 <- axisdf %>% mutate(X.CHROM = as.character(X.CHROM)) %>% 
  mutate(newlab = if_else(X.CHROM %in% c(11, 13, 15, 17, 19, 21), "", X.CHROM))

# get SNPID annotations
minids <- allsig %>% 
  mutate(POS = as.numeric(sapply(str_split(`Chromosome:Position`, ":"), `[`, 2)),
         X.CHROM = as.numeric(sapply(str_split(`Chromosome:Position`, ":"), `[`, 1))) %>% 
  filter(!Taxa == "Class_Saccharomycetes") %>% 
  group_by(Taxa, X.CHROM) %>% 
  slice(which.min(`Genome-Wide Q`)) %>% 
  dplyr::select(taxa = Taxa, SNPID, X.CHROM, POS)

# make df
mandf <- manplotdf %>% full_join(minids)

### build plots in a list and ggarrange to make labels uniform
#mandf <- manplotdf %>% filter(!str_detect(plotnames, "Phenogram"))
manplot <- ggplot(mandf, aes(x = BPcum, y = -log10(P))) +
  # add SNPID annotations
  #geom_text(aes(label = SNPID, x = BPcum + 100000, y = -log10(P)),  hjust = "outward", check_overlap = TRUE, size = 4, na.rm = TRUE) +
  facet_wrap(~plotnames, ncol = 2, nrow = 2) +
  # add horizontal line for genome-wide sig SNP
  geom_hline(yintercept = -log10(6.25e-7), color = "black") +
  # add points
  geom_point(aes(color = as.factor(X.CHROM)), alpha = 0.7, size = 3) +
  # add axis breaks
  scale_x_continuous(label = axisdf1$newlab, breaks = axisdf1$center) +
  # make text bigger
  theme_pubr(base_size = 20) +
  # change labels
  labs(x = "", y ="-log<sub>10</sub>(*P*)") +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "white", color = "white"),
        strip.text = element_markdown(hjust = 0, size = 22),
        #strip.text = element_blank(),
        axis.text.x = element_text(angle = 90),
        axis.title.y = element_markdown(),
        plot.background = element_rect(fill = "white", color = NULL))

## ---- generate locus zoom plots ----

# get eQTL data
myids <- read.table("data/unique_mycAVS_SNPandStructural.txt", sep = "\t", header = TRUE) 
mygte <- read.table("data/gtex_results.txt", header = TRUE, sep = "\t") %>% 
  rename(Gencode = `Gencode.Id`,
         Gene = `Gene.Symbol`,
         Variant = `Variant.Id`,
         SNPID = `SNP.Id`,
         Pval = `P.Value`)  %>% 
  filter(SNPID %in% myids$SNPID)
## get list of taxa
tax <- allsig %>% dplyr::select(Taxa, SNPID, `Variant Type`) %>% filter(!Taxa %in% "Class_Saccharomycetes") %>% drop_na() 
# subset the taxa that had SNPs of eQTLs
taxsub <- tax %>% filter(SNPID %in% mygte$SNPID)
# for plot; 
wtax <- mygte %>% 
  # join
  left_join(taxsub, relationship = "many-to-many")


### make an object that stores the parts of the plots that will all be similar

mytheme <- theme_pubr(base_size = 20) +
  theme(legend.position = "none",
        axis.title.y = element_markdown(),
        plot.background = element_rect(fill = "white", color = NULL))

## build base plot that will be the same for all
topplot <- function(df, hexcode) {
  
  myplot <- ggscatter(df %>% filter(is.na(Gene)), x = "mbp", y = "logp",
                      size = 3, alpha = 0.3,
                      color = hexcode) +
    # add QTLs
    geom_point(data = df %>% filter(!is.na(Gene)), aes(x = mbp, y = logp),
               color = "black", size = 4, alpha = 0.6) +
    # add signifiacnce line
    geom_hline(yintercept = -log10(6.25e-7), color = "black") +
    # change y scale
    scale_y_continuous(limits = c(0, 8.8), breaks = c(0, 2, 4, 6, 8)) +
    # change labels
    labs(x = "", y ="-log<sub>10</sub>(*P*)") +
    # add theme
    mytheme
  
  # return plot object to build on
  return(myplot)
}

## make colors match by chromosome
# get the default ggplot2 colors from manhattan plots
cols <- hue_pal()(22)


### ---- plot1: kazachstania chromosome 16 (CDH13) ----

# get kazachstania: CHROMOSOME 16
plotdf1 <- df %>% 
  filter(taxa %in% ("Genus_Kazachstania")) %>% 
  mutate(pfdr = p.adjust(P, method = "fdr")) %>% 
  filter(X.CHROM %in% c(16)) %>% 
  filter(between(POS, 82500000, 84000000)) %>% 
  mutate(logp = -log10(P)) %>% 
  # calculate kb to make easier to interpret
  mutate(kb = POS / 1000,
         mbp = POS / 1000000) %>% 
  # add qtls
  left_join(wtax %>% dplyr::select(SNPID, Gene, Tissue), by = c("ID" = "SNPID")) %>% 
  # make for plot
  mutate(is.qtl = if_else(!is.na(Gene), Gene, "")) 

# get position for text label of CDH13 (from Ensembl)
lab16 <- mean(c(82.660408, 83.830204))

### build top plot
ch16top <- topplot(plotdf1, cols[16]) +
  # get x axis breaks
  scale_x_continuous(breaks = c(82.5, 84)) +
  # add CDH13
  annotate(geom = "segment", x = 82.660408, xend = 83.830204, y = 8,
           arrow = arrow(length = unit(0.3, "cm"), type = "closed"), color = "green", lineend = "round",
           linewidth = 1.5) +
  geom_text(x = lab16, y = 8.4, label = "CDH13", fontface = "italic", check_overlap = TRUE, size = 6) 


## ---- plot2: kazachstania chromosome 4 (ANAPC10) ----


# get eqtls for chr 4 (multiple genes per qtl - just need the list of snps)
c4qtl <- wtax %>% filter(str_detect(Variant, "chr4")) %>% distinct(SNPID)

# get the dataframe
plotdf2 <- df %>% 
  filter(taxa %in% ("Genus_Kazachstania")) %>% 
  mutate(pfdr = p.adjust(P, method = "fdr")) %>% 
  filter(X.CHROM %in% c(4)) %>% 
  filter(between(POS, 145000000, 146500000)) %>% 
  mutate(logp = -log10(P)) %>% 
  # calculate kb to make easier to interpret
  mutate(kb = POS / 1000,
         mbp = POS / 1000000) %>% 
  # add qtls
  #left_join(wtax %>% dplyr::select(SNPID, Gene, Tissue), by = c("ID" = "SNPID")) %>% 
  # make for plot
  mutate(Gene = if_else(ID %in% c4qtl$SNPID, "is.qtl", NA)) 
# get position for text label of ANAPC10 (from Ensembl)
lab4 <- mean(c(145.888264, 146.019693))

### build top plot
ch4top <- topplot(plotdf2, cols[4]) +
  # get x axis breaks
  scale_x_continuous(breaks = c(145.0, 146.5)) +
  # add ANAPC10
  annotate(geom = "segment", x = 145.888264, xend = 146.019693, y = 8,
           arrow = arrow(length = unit(0.3, "cm"), type = "closed",
                         #### ANAPC10 IS ON REVERSE STRAND
                         ends = "first"), color = "green", lineend = "round",
           linewidth = 1.5) +
  geom_text(x = lab4, y = 8.4, label = "ANAPC10", fontface = "italic", check_overlap = TRUE, size = 6) 


## ----- plot3: pleosporales chromosome 1 (PTPRC) ----

# build dataframe
plotdf3 <- df %>% 
  filter(taxa %in% ("Order_Pleosporales")) %>% 
  mutate(pfdr = p.adjust(P, method = "fdr")) %>% 
  filter(X.CHROM %in% c(1)) %>% 
  #filter(between(POS, 82600000, 83800000)) %>% 
  filter(between(POS, 198000000, 199000000)) %>% 
  mutate(logp = -log10(P)) %>% 
  # calculate kb to make easier to interpret
  mutate(kb = POS / 1000,
         mbp = POS / 1000000) %>% 
  # add qtls
  left_join(wtax %>% dplyr::select(SNPID, Gene, Tissue), by = c("ID" = "SNPID")) %>% 
  # make for plot
  mutate(is.qtl = if_else(!is.na(Gene), Gene, ""))

# get position for text label of ANAPC10 (from Ensembl)
lab1 <- mean(c(198.607801, 198.726545))

### build top plot
ch1top <- topplot(plotdf3, cols[1]) +
  # get x axis breaks
  scale_x_continuous(breaks = c(198.0, 198.8)) +
  # add PTPRC
  annotate(geom = "segment", x = 198.607801, xend = 198.726545, y = 8.5,
           arrow = arrow(length = unit(0.3, "cm"), type = "closed"), color = "green", lineend = "round",
           linewidth = 1.5) +
  geom_text(x = lab1, y = 8.8, label = "PTPRC", fontface = "italic", check_overlap = TRUE, size = 6) 

### ---- arrange all ----

zoomplot <- ggarrange(ch1top, ch4top, ch16top, ncol = 3)
### add together 
final.plot <- ggarrange(manplot, zoomplot, 
                        ncol = 1, heights = c(1, 0.6))
ggsave(final.plot, filename = "figures/most_gwas_v1.png", height = 17, width = 20, units = "in", dpi = 300)

## ---- build Kazachstania boxplot ----

# get genotypes
#### RESTRICTED DATA
gen <- read.table("gwas-out/extractsnps.ped", header = FALSE) %>% 
  dplyr::select(V1, V7, V8, V9, V10)
names(gen) <- c("IID", "rs12149890_a1", "rs12149890_a2", "rs12929586_a1", "rs12929586_a2")

gen <- gen %>% 
  mutate(rs12149890_geno = paste0(rs12149890_a1, rs12149890_a2),
         rs12929586_geno = paste0(rs12929586_a1, rs12929586_a2))


### get fungal data
load("data/phylo_ITS_resolvedNA.RData")
psf <- psnewname
# subset gen 
gensub <- gen %>% 
  filter(IID %in% sample_names(psf))

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

# join and plot
forplot <- gensub %>% 
  #full_join(kz) %>% 
  full_join(allclr %>% select(Genus_Kazachstania) %>% rownames_to_column(var = "IID") %>% mutate(IID = as.numeric(IID))) %>% 
  mutate(fullgeno = paste(rs12149890_geno, rs12929586_geno, sep = "/")) %>% 
  mutate(fullgeno = if_else(fullgeno == "TG/AC", "GT/CA", fullgeno)) %>% 
  mutate(fullgeno = factor(fullgeno, ordered = TRUE, levels = c("GG/CC", "GT/CA", "TT/AA")))

## build violin plot
ggviolin(forplot, x = "fullgeno", y = "Genus_Kazachstania", add = "jitter", size = 2, fill = "lightgrey",
         add.params = list(color = "darkblue", alpha = 0.5, size = 5),
         xlab = "rs12149890 / rs12929586", ylab = "*Kazachstania* abundance") +
  stat_pvalue_manual(stats, label = "p.adj.signif", hide.ns = TRUE, size = 10, bracket.size = 1, bracket.nudge.y = 1) +
  theme_pubr(base_size = 20) +
  theme(axis.title.y = element_markdown()) 
ggsave(filename = "figures/kz_phewassnps_violin.png", dpi = 600, height = 6, width = 6, units = "in")