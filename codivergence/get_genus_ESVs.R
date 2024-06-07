#### make fasta files for each genus's exact sequence variants
# EVS 3/2024

library(Biostrings)
library(tidyverse)
library(phyloseq)
library(microViz)
library(ggpubr)

## load phyloseq
load("private/hominid_phyloITS.RData")
psfix <- psf %>% tax_fix()

# get genera
gen <- get_taxa_unique(psfix, "Genus")

# get names of fungal species
tt <- psfix %>% tt_get() %>% as.data.frame() %>% 
  mutate(Species = str_replace_all(Species, " ", "-")) %>% 
  rownames_to_column(var = "OTUID")

## get file to match each read to its assigned OTU
uc <- read.table("data/usearch.uc", sep = "\t", header = FALSE) %>% 
  dplyr::select(V1, V9, V10)
names(uc) <- c("HIT", "READ", "OTU")
dat <- uc %>% 
  filter(HIT == "H") %>% 
  mutate(otulab = str_extract(OTU, "OTU_(\\d){1,10}"))

# get the file of all sequences
allfa <- readDNAStringSet("data/all.fasta")

# get metadata for hominid name
meta <- readxl::read_xlsx("private/hominid_MetadataITS.xlsx") %>% 
  mutate(hominid = case_when(
    str_detect(Group, "Human") ~ "Human",
    str_detect(Group, "Chimp") ~ "Chimp",
    str_detect(Group, "Lowland Gorilla") ~ "LowGorilla",
    str_detect(Group, "Mountain Gorilla") ~ "MountGorilla"
  ))

# create a "master list" of sequence IDs, OTUs, and hominid species
ids <- names(allfa)
iddf <- data.frame(fullname = ids) %>% 
  mutate(FileID = str_remove(fullname, "_S(\\d){1,5}\\.(\\d){1,5};size=(\\d){1,10}")) %>% 
  # join metadata to get hominid
  left_join(meta %>% dplyr::select(FileID, hominid)) %>% 
  # join seqs to get OTU
  left_join(dat, by = c("fullname" = "READ")) %>% 
  # join tax table to get fungi name
  left_join(tt %>% dplyr::select(OTUID, Species), by = c("otulab" = "OTUID")) %>% 
  # drop NA OTUs (some OTUs were removed during filtering etc)
  drop_na(Species) %>% 
  # make new names for output files
  mutate(halfname = str_remove(fullname, ";size=(\\d){1,5}")) %>% 
  mutate(newids = paste(halfname, otulab, Species, hominid, sep = ";")) %>% 
  mutate(halfid = paste(otulab, Species, hominid, sep = ";"))

# set output directory
outdir <- "output-dir/"
if(!dir.exists(outdir)) {dir.create(outdir)}

# ---- make function ----

myHaplo <- function(genus) {
  
  # subset the phyloseq
  sub <- psfix %>% tax_select(genus, "Genus", strict_matches = TRUE, n_typos = 0)
  
  # get OTUs
  otus <- taxa_names(sub)
  
  # subset list of IDs
  newids <- iddf %>% filter(otulab %in% otus)
  
  ### similarly to OTU-level, ensure each hominid is represented at least once (otherwise can't have codiv and are wasting comparisons!) 
  if(length(unique(newids$hominid)) < 4) {
    message("\n Genus ", genus, " does not have all hominids; it has only ", unique(newids$hominid), " ...skipping... \n")
  } else {
    
    # if all four hominids are represented, proceed
    
    # subset sequences
    seqs <- allfa[names(allfa) %in% newids$fullname ]
    # sort to match (ugh)
    seqs <- seqs[match(names(seqs), newids$fullname)]
    
    # rename sequences with fungi, OTU name, and hominid
    names(seqs) <- newids$halfid
    
    # remove spaces from family
    famname <- gsub(" ", "_", genus)
    
    # write out to file
    writeXStringSet(seqs, filepath = paste0(outdir, "/haplotypes_", famname, ".fasta"),
                    append = FALSE)
    
  }
  
}

### ---- run for haplotypes ----

# set output directory
outdir <- "output-dir"
if(!dir.exists(outdir)) {dir.create(outdir)}

# run function
lapply(gen, function(x) myHaplo(x))

### ---- check output statistics ----

# get stats for haplotypes
hapstat <- read.table("stats_genushaplotypes.txt", sep = "\t", header = TRUE) %>% 
  mutate(id = sapply(str_split(file, "haplotypes_"), `[`, 2))

# get stats for derep
derep <- read.table("updated_methods/codiv/stats_derepgenushaplotypes.txt", sep = "\t", header = TRUE) %>% 
  mutate(id = sapply(str_split(file, "haplotypes_"), `[`, 2))

# join
all <- hapstat %>% 
  dplyr::select(num_seqs, sum_len, id) %>% mutate(step = "raw") %>% 
  rbind(derep %>% dplyr::select(num_seqs, sum_len, id) %>% mutate(step = "derep"))

# pivot
allh <- all %>% pivot_wider(names_from = step, values_from = c(num_seqs, sum_len))

# loss from dereplication
summary(allh$num_seqs_derep/allh$num_seqs_raw) # mean 83% retention; min 58%, max 98%
allh %>% mutate(keepsum = num_seqs_derep/num_seqs_raw) %>% 
  filter(keepsum < 0.7) # Dipodascaceae, Pichiaceae, Geotrichum, Trichosporon, and Saccharomycetales Order 

# how big are these files?
ggdensity(all, x = "num_seqs", fill = "step") # most are < ~2500 sequences
all %>% group_by(step) %>% get_summary_stats(c(num_seqs, sum_len), type = "five_number")
# max derep'd is 4946 seqs; will be a BIG alignment  

