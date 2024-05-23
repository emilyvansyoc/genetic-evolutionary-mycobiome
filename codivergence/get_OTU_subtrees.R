### get OTU subtrees for PACo
# EVS 3/2024

library(ape)
library(tidyverse)

## directory to trees
tdir <- "updated_methods/codiv/genus_trees/NJ_JC_negbl_04-09-2024/"
gtrees <- list.files(tdir, pattern = ".newick", full.names = TRUE)
allgtrees <- list()

# read each tree
alltrees <- list()
for(i in 1:length(gtrees)) {
  
  # get tree
  tr <- read.tree(gtrees[i])
  allgtrees[[i]] <- tr
  
  # get tip labels
  tips <- tr$tip.label
  
  # get OTUs
  otus <- data.frame(
    tiplabs = tr$tip.label
  ) %>% 
    mutate(OTU = str_extract(tiplabs, "OTU_(\\d){1,10}"),
           tip1 = str_remove(tiplabs, "_(\\d){1,10}$"),
           hominid = str_extract(tip1, "[:alpha:]+$"))
  num <- unique(otus$OTU)
  
  subtrees <- list()
  for(j in 1:length(num)) {
    
    # subset
    sub <- otus %>% filter(OTU == num[j])
    # get all 4 hominids
    if(length(unique(sub$hominid)) == 4) {
      
      # subset tree
      subtr <- keep.tip(tr, tip = sub$tiplabs)
      
      # add to list for paco
      subtrees[[j]] <- subtr
      
      names(subtrees)[[j]] <- num[j]
    } # close if statement 
    #else {subtrees[[j]] <- NULL}
    
  } # close j loop
  
  # remove NULL trees
  subtrees <- subtrees[lengths(subtrees) > 0]
  
  ## combine "big" list
  alltrees <- c(alltrees, subtrees)
  
} # close i loop

## save for PACo run
save(alltrees, file = "updated_methods/codiv/all_otu_subtrees_04-10-2024_nonegbl.RData")


