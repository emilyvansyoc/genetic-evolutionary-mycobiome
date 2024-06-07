## run Parafit
# EVS 4/2024


library(ape)
library(vegan)
library(tidyverse)
#library(phyloseq)
library(microViz)
library(Biostrings)
library(purrr)


## get trees
load("data/all_otu_subtrees_04-10-2024_nonegbl.RData")
# get host cophenetic
load("data/host_cophenetic.RData")

# loop
paralist <- list()
paradf <- data.frame()

for(i in 1:length(alltrees)) {
  
  # print
  cat("\n working on tree ", i, names(alltrees)[i])
  
  # get tree
  mytree <- alltrees[[i]]
  
  # get distances
  mydis <- cophenetic.phylo(mytree)
  
  
  # make the HP matrix
  mat <- data.frame(p = mytree$tip.label) %>% 
    mutate(p = str_remove_all(p, "'"), 
           h = str_extract(p, "[:alpha:]+_[0-9]*$")) %>% 
    ## wonky wrangling to capture hominid name given unique renaming scheme
    mutate(h = if_else(is.na(h), str_extract(p, "[:alpha:]+$"), h)) %>% 
    mutate(h = str_remove(h, "_(\\d){1,5}")) %>% 
    # binary-ize hominid membership
    mutate(Chimp = if_else(h == "Chimp", 1, 0),
           Human = if_else(h == "Human", 1, 0),
           LowGorilla = if_else(h == "LowGorilla", 1, 0),
           MountGorilla = if_else(h == "MountGorilla", 1, 0)) %>% 
    dplyr::select(-h) %>% 
    column_to_rownames(var = "p") %>% 
    as.matrix()
  
  # calculate parafit and skip errors with bad fit
  tryCatch({
    myfit <- parafit(hostd, mydis, mat, correction = "cailliez", #test.links = TRUE, 
                     nperm = 99, silent = TRUE, seed = 395)
    # add to the output list
    paralist[[i]] <- myfit
    names(paralist)[i] <- names(alltrees)[[i]]
    # dataframe-ize
    paradf <- rbind(paradf,
                    data.frame(otu = names(alltrees)[[i]],
                               pval = myfit$p.global,
                               parastat = myfit$ParaFitGlobal))
    
  },
  error = function(e){cat("\n WARNING: Parafit does not fit subtree", i, "\n")})
  
  
}

# save output
save(paradf, file = "data/parafit_results_OTUs_nonegbl_04-10-2024.RData")
