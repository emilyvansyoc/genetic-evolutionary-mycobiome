### make neighbor-joining and UPGMA trees for alignments
# make genus-level alignments
# EVS 04/2024

library(ape)
library(tidyverse)
library(doParallel)
library(foreach)
library(parallel)

## ---- build parallel loop for each alignment fasta ----

dir <- "path/to/align_genus_renamed_04-09-2024/"
myfi <- list.files(dir, pattern = "clipkit", full.names = TRUE)

## remove genera that cannot be aligned with ITS: Fusarium, Asperillus, Penicillum, Cortinarius
# Cortinarius are not in the dataset
myfi <- myfi[!str_detect(myfi, "unique_clipkit_derep_haplotypes_Aspergillus.fasta")]
myfi <- myfi[!str_detect(myfi, "unique_clipkit_derep_haplotypes_Fusarium.fasta")]
myfi <- myfi[!str_detect(myfi, "unique_clipkit_derep_haplotypes_Penicillium.fasta")]

#### FOR NOW: remove "big" tree (Pichia) and run on command line
#myfi <- myfi[myfi != c("unique_clipkit_derep_haplotypes_Pichia.fasta")]

## output directories to write tree files
njdir <- "path/to/genus_trees/NJ_JC_04-09-2024/"
if(!dir.exists(njdir)) {dir.create(njdir)}
#dgdir <- "updated_methods/codiv/genus_trees/UPGMA_JC_renamed/"
#if(!dir.exists(dgdir)) {dir.create(dgdir)}


### RUN LOOP IN PARALLEL
n.cores <- parallel::detectCores() - 2
# create a backend
#create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "FORK",
  outfile = "path/to/make-genus-trees-stdout.log" # prints stdout
)

#check cluster definition (optional)
print(my.cluster)
#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)
#check if it is registered (optional)
foreach::getDoParRegistered()

#### RUN LOOP
foreach(i = 1:length(myfi),
        .errorhandling = "pass" # I have my own error catch below so errors in the loop can be passed 
) %dopar% {
  
  myname <- str_remove(myfi[i], "/unique_clipkit_derep_haplotypes_")
  myname <- str_remove(myname, ".fasta")
  
  # print progress
  cat(" \n working on file ", myname, " ... \n")
  
  # read in alignment
  aln <- read.dna(myfi[i], format = "fasta")
  
  # create distance matrix with Jukes Cantor model
  # pairwise deletion removes the site with missing data for all sequences (true gaps)
  # this is necessary because tree building will fail
  dist <- dist.dna(aln, model = "JC69", pairwise.deletion = TRUE)
  
  tryCatch(
    
    {
      # create neighbor joining tree
      # use "njs" to allow missing values 
      tree <- njs(dist)
      
      
      # write to file
      ## pretty-fy file name
      
      write.tree(tree, file = paste0(njdir, "/NJtree_", myname, ".newick"))
      #write.tree(dend, file = paste0(dgdir, "/UPGMAtree_", myname, ".newick"))
      cat("\n done with tree ", myname, "\n")
    },
    error = function(e) {cat("\n Alignment", myfi[i], " cannot form one or both trees \n")}
    
  ) # close tryCatch
  
  
} # close parallel loop

### close cluster
parallel::stopCluster(cl = my.cluster)
