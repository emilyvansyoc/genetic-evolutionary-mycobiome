## replace negative branch lengths in NJ trees
# EVS 4/2024

library(ape)
library(tidyverse)

## directory to NJ trees
njdir <- "~/bioinformatics/r-projects/hominid/updated_methods/codiv/genus_trees/NJ_JC_04-09-2024/"

## output directory to 'fixed' trees
outdir <- "~/bioinformatics/r-projects/hominid/updated_methods/codiv/genus_trees/NJ_JC_negbl_04-09-2024/"
if(!dir.exists(outdir)) {dir.create(outdir)}

## get trees
myfi <- list.files(njdir, full.names = TRUE)

## loop through
for(i in 1:length(myfi)) {
  
  myname <- str_remove(myfi[i], "/Users/epb5360/bioinformatics/r-projects/hominid/updated_methods/codiv/genus_trees/NJ_JC_04-09-2024//NJtree_")
  myname <- str_remove(myname, ".newick")
  
  # read tree
  tr <- read.tree(myfi[i])
  
  # replace negative branch edges
  tr$edge.length[tr$edge.length < 0] <- 0
  
  # resolve these now zero lengths into polytomies
  td <- di2multi(tr)
  
  # print how many internal nodes were resolved
  cat("\n tree", myname, "resolved", tr$Nnode - td$Nnode, "internal nodes \n ")
  
  # write out the resolved tree
  write.tree(td, file = paste0(outdir, "resolved_", myname, ".newick"))
  
  
}
