## get run IDs for WGS stool of the individuals that we have fungal sequencing data for 
# EVS 3/2024

library(tidyverse)
library(phyloseq)
library(microViz)

## get fungal ID's
ids <- readxl::read_xlsx("SRAMetadata_Nash2017_PRJNA356769.xlsx") %>% 
  filter(!str_detect(`Library Name`, "18S")) %>% 
  select(Run, `Sample Name`) %>% 
  mutate(id = sapply(str_split(`Sample Name`, "_"), `[`, 1)) %>% 
  mutate(visit = sapply(str_split(`Sample Name`, "_"), `[`, 2)) %>% 
  column_to_rownames(var = "Run") %>% 
  mutate(idvisit = paste(id, visit, sep = "_"))

# subset for individual IDs
fids <- ids %>% dplyr::select(id) %>% distinct() %>% rename(RANDSID=id) %>% mutate(RANDSID = as.numeric(RANDSID))

## get WGS metadata
met <- read.table("runinfo_PRJNA48479.csv", sep = ",", header = TRUE) %>% 
  dplyr::select(Run, Sample, SampleName, Subject_ID, Analyte_Type, Model, CenterName, BioSample) %>% 
  column_to_rownames(var = "Run") 
# get just stool
mets <- met %>% 
  filter(Analyte_Type == "G_DNA_Stool")

#### wrangle to match subject ID to RSID
## get sample key (available on dbGaP)
key <- read.table("data/HMP_samplekey.txt", sep = "\t", header = TRUE) 
## match to the BioSample Accession
dat <- mets %>% 
  rownames_to_column(var = "Run") %>% 
  left_join(key, by = c("BioSample" = "BioSample.Accession"))

## match the WGS RANDSID to fungal RANDSID
fmatch <- dat %>% 
  semi_join(fids) # 89 total matches 

length(unique(fmatch$RANDSID)) # 89 unique individuals
length(unique(fmatch$Run)) # 1129
length(unique(fmatch$Sample)) # 182 unique stool samples

# save these Runs   
write.table(fmatch$Run, file = "wgs-analyses/runids_stool_WGS_paired.txt", sep = "\t", quote = FALSE, row.names = FALSE)

## write just the HiSeq samples
hi <- fmatch %>% 
  filter(Model == "Illumina HiSeq 2000") # 959 total runs
length(unique(hi$RANDSID))
# save these
write.table(hi$Run, file = "wgs-analyses/runids_stool_WGS_paired_HiSeq.txt", sep = "\t", quote = FALSE, row.names = FALSE)

## write just the Genome analyzer samples
ga <- fmatch %>% 
  filter(Model == "Illumina Genome Analyzer II")

length(unique(ga$RANDSID)) # 49 unique individuals
length(unique(ga$Run)) # 153 runs
length(unique(ga$Sample)) # 79 unique stool samples

# save
write.table(ga$Run, file = "wgs-analyses/runids_stool_WGS_paired_GAII.txt", sep = "\t", quote = FALSE, row.names = FALSE)

