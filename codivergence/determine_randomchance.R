### random chance permutations
# determine the probability of detecting codivergence by random chance
# EVS 4/2024

library(tidyverse)
library(ape)
library(paco)
library(foreach)
library(doParallel)
library(parallel)


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

### ---- set up loop: run parafit first (faster!) ----

# set working directory
dir <- "/storage/group/srb6251/default/emilyvansyoc/hominid_codiv/"

# get host cophenetic distances
load(paste0(dir, "host_cophenetic.RData"))

# get OTU trees
load(paste0(dir, "all_otu_subtrees_04-10-2024_nonegbl.RData"))

# set output directory for results
outdir <- paste0(dir, "random_chance_results/parafit_all")
if(!dir.exists(outdir)) {dir.create(outdir)}

#### RUN PARALLEL LOOP 
foreach(i = 1:length(alltrees),
        .errorhandling = "pass" # I have my own error catch below so errors in the loop can be passed 
) %dopar% {

# print
 cat("\n working on tree ", i, names(alltrees)[[i]])

#### CREATE SHUFFLED TREES
 make tree to shuffle labels
 shuftree <- list()
# subset tree
 tr <- alltrees[[i]]

# permute the shuffling and make big list
 for(j in 1:100) {

# randomly shuffle the tip labels
randtips <- sample(x = tr$tip.label, size = length(tr$tip.label), replace = FALSE)
randtree <- tr
randtree$tip.label <- NULL
randtree$tip.label <- randtips

# save to list
shuftree[[j]] <- randtree
names(shuftree)[[j]] <- paste0(names(alltrees)[[i]], "_iteration", j)

} # close first j loop

#### RUN PARAFIT 
paradf <- data.frame()

for(j in 1:length(shuftree)) {



# get tree
 mytree <- shuftree[[j]]

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
# dataframe-ize
 paradf <- rbind(paradf,
                 data.frame(otu = names(shuftree)[[j]],
                           pval = myfit$p.global,
                          parastat = myfit$ParaFitGlobal))

  },
 error = function(e){cat("\n WARNING: Parafit does not fit subtree", j, "\n")})


  } # close second j loop

# write out the paradf object
 write.table(paradf, file = paste0(outdir, "/parafit_", names(alltrees)[[i]], ".txt"))

# print message
cat("\n PARAFIT FINISHED: ", i, names(alltrees)[[i]])

} # close parallel loop

### ---- set up parallel loop: run PACo -----

# set output directory for results
outdirP <- paste0(dir, "random_chance_results/paco_all")
if(!dir.exists(outdirP)) {dir.create(outdirP)}

#### RUN PARALLEL LOOP 
foreach(i = 1:length(alltrees),
        .errorhandling = "pass" # I have my own error catch below so errors in the loop can be passed 
) %dopar% {
  
  # print
  cat("\n working on tree ", i, names(alltrees)[[i]])
  
  #### CREATE SHUFFLED TREES
  # make tree to shuffle labels
  shuftree <- list()
  # subset tree
  tr <- alltrees[[i]]
  
  # permute the shuffling and make big list
  for(j in 1:100) {
    
    # randomly shuffle the tip labels
    randtips <- sample(x = tr$tip.label, size = length(tr$tip.label), replace = FALSE)
    randtree <- tr
    randtree$tip.label <- NULL
    randtree$tip.label <- randtips
    
    # save to list
    shuftree[[j]] <- randtree
    names(shuftree)[[j]] <- paste0(names(alltrees)[[i]], "_iteration", j)
    
  } # close first j loop
  
  #### RUN PACo
  pacodf <- data.frame()
  
  for(j in 1:length(shuftree)) {
    
    # get tree
    btree <- shuftree[[j]]
    
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
      # save global params to text file
      pacodf <- rbind(pacodf,
                      data.frame(subtree = names(shuftree)[[j]],
                                 paco.ss = mypac$gof$ss,
                                 paco.p = mypac$gof$p) )
      
      
    },
    error = function(e){cat("\n WARNING: PACo does not fit subtree", j, "\n")})
    
    
    
  } # close second j loop
  
  # write out the paco object
  write.table(pacodf, file = paste0(outdirP, "/paco99perm_", names(alltrees)[[i]], ".txt"))
  
  # print message
  cat("\n PACO FINISHED: ", i, names(alltrees)[[i]])
  
} # close parallel loop

# stop cluster
parallel::stopCluster(cl = my.cluster)(base) 