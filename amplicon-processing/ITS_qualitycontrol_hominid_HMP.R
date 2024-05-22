### quality control for ITS 
# EVS 10/2023

library(dada2)
library(tidyverse)
library(Biostrings)
library(ShortRead)

## pathways
path <- "~/bioinformatics/seq-reads/gomez-hominid-raw-reads/Primate_Data/ITS/"
list.files(path)

# get samples in dada2 format
fnFs <- sort(list.files(path, pattern = "R1_001-ITS.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "R2_001-ITS.fastq.gz", full.names = TRUE))

# list primers
ITS3 <- "GCATCGATGAAGAACGCAGC" # the primer used in the paper
ITS4 <- "TCCTCCGCTTATTGATATGC"
ITS3Meta <- "TCGATGAAGAACGCAGCG" # the primer used in the metadata

# get correct orientation of the primers
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  #require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, 
               Complement = Biostrings::complement(dna), 
               Reverse = Biostrings::reverse(dna),
               RevComp = Biostrings::reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
ITS3.orients <- allOrients(ITS3)
ITS4.orients <- allOrients(ITS4)
ITS3M.orients <- allOrients(ITS3Meta)

## ---- first remove N's -----

fnFs.filtN <- file.path(path, "filtN", basename(fnFs))
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))

# filter
remN <- filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE,
                      verbose = TRUE)

# re-set the path for N-removed files
npath <- "~/bioinformatics/seq-reads/gomez-hominid-raw-reads/Primate_Data/ITS/filtN/"

## ---- get primer hits ----

# get number of primer hits
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

# get primers for all reads
for(i in 1:length(fnFs.filtN)) {
  
  # print sample name
  cat(paste0("\n sample ", fnFs.filtN[i], "\n"))
  
  print(
    # print primers
    rbind(FWD.ForwardReads = sapply(ITS3.orients, primerHits, fn = fnFs.filtN[[i]]), 
          FWD.ReverseReads = sapply(ITS3.orients, primerHits, fn = fnRs.filtN[[i]]),
          REV.ForwardReads = sapply(ITS4.orients, primerHits, fn = fnFs.filtN[[i]]), 
          REV.ReverseReads = sapply(ITS4.orients, primerHits, fn = fnRs.filtN[[i]]),
          FWDM.ForwardReads = sapply(ITS3M.orients, primerHits, fn = fnFs.filtN[[i]]),
          FWDM.ReverseReads = sapply(ITS3M.orients, primerHits, fn = fnRs.filtN[[i]]))
  )
  
}

## this shows that the Rev are in forward orientation in Rev reads
# FwdM primer in Forward orientation in forward reads
# Rev primers are in rev complement in Forward reads
# FWDM primers are in reverse complement in Reverse reads 
# anywhere from a few hundred to several thousand

## the ITS3 primer for the paper is not found -> use the ITS3Meta primer instead
# this mirrors the expected findings from the dada2 tutorial 

## ---- remove primers with cutadapt in R ----

# location of cutadapt 
cutadapt <- "/Users/epb5360/micromamba/envs/bioinfo/bin/cutadapt"
system2(cutadapt, args = "--version") # check this is accessible

# set pathway location
path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs.filtN))
fnRs.cut <- file.path(path.cut, basename(fnRs.filtN))

# set primer directions
FWD.RC <- dada2:::rc(ITS3Meta)
REV.RC <- dada2:::rc(ITS4)

# trim FWD and reverse complement of REV off of forward reads
R1.flags <- paste("-g", ITS3Meta, "-a", REV.RC)
# trim REV and the reverse complement of FWD off of reverse reads
R2.flags <- paste("-G", ITS4, "-A", FWD.RC)

# run cutadapt 
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i], "-j 8"# number of cores
  )) # input files
}

## sanity check: re-run the primer check
for(i in 1:length(fnFs.cut)) {
  
  # print sample name
  cat(paste0("\n sample ", fnFs.cut[i], "\n"))
  
  print(
    # print primers
    rbind(FWD.ForwardReads = sapply(ITS3.orients, primerHits, fn = fnFs.cut[[i]]), 
          FWD.ReverseReads = sapply(ITS3.orients, primerHits, fn = fnRs.cut[[i]]),
          REV.ForwardReads = sapply(ITS4.orients, primerHits, fn = fnFs.cut[[i]]), 
          REV.ReverseReads = sapply(ITS4.orients, primerHits, fn = fnRs.cut[[i]]),
          FWDM.ForwardReads = sapply(ITS3M.orients, primerHits, fn = fnFs.cut[[i]]),
          FWDM.ReverseReads = sapply(ITS3M.orients, primerHits, fn = fnRs.cut[[i]]))
  )
  
}

## looks good!

## ---- do filter and trim -----

# make file paths
filtFs <- file.path(path, "filtered", basename(fnFs))
filtRs <- file.path(path, "filtered", basename(fnRs))

# filter and trim
trims <- filterAndTrim(fnFs.cut, filtFs, fnRs.cut, filtRs,
                       #maxEE = c(2, 2), 
                       #rm.phix = TRUE,
                       truncQ = 30, 
                       #minLen = 50, # enforced AFTER trimming -> default is 20
                       #minQ = 15, # after truncation, remove reads with quality score < 30
                       multithread = 10, verbose = TRUE
)
# how much is lost?
trims <- data.frame(trims) %>% 
  mutate(kept = (reads.out / reads.in)* 100)
head(trims)
sum(trims$reads.out)# not super applicable - includes elephants

## ---- compare number of reads to original paper ----

# get metadata for sample names to remove elephants (& blanks?)

#### get sample metadata
meta <- readxl::read_xlsx("~/bioinformatics/gomez-hominid/EVS_MetadataITS.xlsx")

### add two blank files to the metadata file
bls <- data.frame(SampleID = c("Blank1", "Blank2"),
                  FileID = c("Blank1_009_A09_011", "Blank2_009_A10_011"),
                  SampleName = "Blanks",
                  Group = "Blanks",
                  Captivity_Status = NA,
                  Wonky = "",
                  Dataset = "Blanks")
meta <- meta %>% rbind(bls)

## add
df <- trims %>% 
  rownames_to_column(var = "ID") %>% 
  mutate(FileID = str_remove_all(ID, "_S(\\d){1,4}_R1_001-ITS.fastq.gz")) %>% 
  left_join(meta)
# remove elephants
dat <- df %>% 
  filter(!Group == "Elephant") 
sum(dat$reads.out) # 3,814,090
