### run PACo in a parallel loop (submitted to an HPC)

library(foreach)
library(doParallel)
library(parallel)
library(tidyverse)
library(ape)
library(paco)

### STDIN: first argument is number of cores, second is log file for stdout

### ---- set up cluster ----

# this sets number of cores
args <- commandArgs(TRUE)
n.cores <- args[[1]]

# start cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "FORK",
  outfile = args[[2]] # prints stdout
)

#check cluster definition (optional)
print(my.cluster)
#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)
#check if it is registered (optional)
foreach::getDoParRegistered()


## ---- get input data -----

# set working directory
dir <- "my-dir"

# get host cophenetic distances
load(paste0(dir, "host_cophenetic.RData"))

# get list of trees
load(paste0(dir, "all_otu_subtrees_04-09-2024.RData"))

# set output directory for PACo results
pacdir <- paste0(dir, "PACo_OTU_quasi_p100_04-09-2024/")
if(!dir.exists(pacdir)) {dir.create(pacdir)}

### ---- RUN LOOP ----
foreach(i = 1:length(alltrees),
        .errorhandling = "pass" # I have my own error catch below so errors in the loop can be passed 
) %dopar% {
  
  
  # get the basename of the tree for matching file
  bname <- names(alltrees)[[i]]
  
  # print progress
  cat(" \n STARTING subtree", i, bname, "...")
  
  # read in biont tree
  btree <- alltrees[[i]]
  
  # calculate cophenetic distance
  bdis <- cophenetic.phylo(btree)
  
  # make presence/absence matrix
  hpmat <- data.frame(
    tiplab = btree$tip.label
  ) %>% 
    mutate(tip1 = str_remove(tiplab, "_(\\d){1,4}$"),
           hominid = str_extract(tip1, "[:alpha:]+$"),
           member = 1) %>% 
    dplyr::select(-tip1) %>% 
    pivot_wider(names_from = hominid, values_from = member, values_fill = 0) %>% 
    column_to_rownames(var = "tiplab") %>% 
    as.matrix() %>% t()
  
  # run through the PACo commands
  tryCatch({
    
    d <- prepare_paco_data(H = hostd, P = bdis, HP = hpmat)
    coord <- add_pcoord(d, correction = "cailliez")
    mypac <- PACo(coord, nperm = 99, symmetric = TRUE, proc.warnings = FALSE, method = "quasiswap")
    # get individual links
    mylink <- paco_links(mypac, proc.warnings = FALSE)
    # get residuals
    res <- residuals_paco(mypac$proc)
    
    # write output to file
    save(mylink, file = paste0(pacdir, "/PACoOutput_OTU_M21_quasi_p99_", bname, ".RData"))
    save(res, file = paste0(pacdir, "/PACoResiduals_OTU_M21_quasi_p99_", bname, ".RData"))
    
    # save global params to text file
    write.table(data.frame(subtree = bname,
                           paco.ss = mypac$gof$ss,
                           paco.p = mypac$gof$p),
                file = paste0(pacdir, "/PACoMetrics_OTU_M21_quasi_p99_", bname, ".txt" ),
                row.names = FALSE, sep = "\t")
    
    ## print finished message to STOUD
    cat("\n FINISHED with Subtree ", i, bname)
    
  },
  error = function(e){cat("\n WARNING: PACo does not fit subtree", i, "\n")})
  
  
} # close parallel loop  

# stop cluster
parallel::stopCluster(cl = my.cluster)