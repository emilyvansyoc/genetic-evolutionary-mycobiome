### function to get p values from topo cong test; host vs dendrogram
# EVS 11/2023

library(tidyverse)
library(ape)

myPval <- function(dist.metric,
                   stoch.filename = "stoch.txt",
                   normalized = TRUE,
                   path) {
  
  ############ get files #################
  
  # get files in PATH
  tabs <- list.files(path, full.names = TRUE)
  
  # get file that matches arguments
  
  myobs <- tabs[str_detect(tabs, paste0("obs_", dist.metric))]
  myst <-tabs[str_detect(tabs, stoch.filename)]
  
  # read files
  obstab <- read.table(myobs, sep = "\t", header = TRUE)
  sttab <- read.table(myst, sep = "\t", header = TRUE)
  
  ########### get scores ####################
  if(normalized == TRUE) {
    
    cat("\n showing NORMALIZED scores and p values \n ")
    
    # subset to get just normalized data
    obssub <- obstab %>% dplyr::select(ends_with("_toUnifAvg"))
    stsub <- sttab %>% dplyr::select(ends_with("_toUnifAvg"))
    
    outdf <- data.frame()
    ### for each metric, get scores
    for(i in 1:length(colnames(stsub))) {
      
      # get observed 
      obs.score <- obssub[,i]
      # get stochastic scores 
      st.sc <- length(which(stsub[,i] <= obs.score))
      
      ## get p values
      p <- st.sc / nrow(stsub)
      
      # print metrics
      cat("\n -------------------------------------------")
      cat("\n total number of stochastic comparisons: ", nrow(stsub))
      cat("\n observed ", colnames(obssub)[i], ": ", obs.score)
      cat("\n p value ",  colnames(obssub)[i], ": ", p)
      cat("\n -------------------------------------------")
      
      ## build results dataframe
      res <- data.frame(
        no.comps = nrow(stsub),
        metric = colnames(obssub)[i],
        obsscore = obs.score,
        pval = p,
        meanstoch = summary(stsub)[,i][4],
        minstoch = summary(stsub)[,i][1],
        maxstoch = summary(stsub)[,i][6]
      )
      
      outdf <- rbind(outdf, res)
    }
  } else {
    
    cat("\n showing UNNORMALIZED scores and p values \n ")
    
    # subset to get just normalized data
    obssub <- obstab %>% dplyr::select(!ends_with("Avg")) %>% 
      dplyr::select(!ends_with("Avg.1")) %>% 
      dplyr::select(-c("No", "RefTree", "Tree"))
    stsub <- sttab %>% dplyr::select(!ends_with("Avg")) %>% 
      dplyr::select(!ends_with("Avg.1")) %>% 
      dplyr::select(-c("No", "RefTree", "Tree"))
    
    outdf <- data.frame()
    ### for each metric, get scores
    for(i in 1:length(colnames(stsub))) {
      
      # get observed 
      obs.score <- obssub[,i]
      # get stochastic scores 
      st.sc <- length(which(stsub[,i] <= obs.score))
      
      ## get p values
      p <- st.sc / nrow(stsub)
      
      # print metrics
      cat("\n -------------------------------------------")
      cat("\n total number of stochastic comparisons: ", nrow(stsub))
      cat("\n observed ", colnames(obssub)[i], ": ", obs.score)
      cat("\n p value ",  colnames(obssub)[i], ": ", p)
      cat("\n -------------------------------------------")
      
      ## build results dataframe
      res <- data.frame(
        no.comps = nrow(stsub),
        metric = colnames(obssub)[i],
        pval = p,
        meanstoch = summary(stsub)[,i][4],
        minstoch = summary(stsub)[,i][1],
        maxstoch = summary(stsub)[,i][6]
      )
      
      outdf <- rbind(outdf, res)
    }
  }
  
  ### return observed scores
  # remove weird text
  outdf <- outdf %>% 
    mutate(across(ends_with("stoch"), ~str_extract(., "(\\d)\\.(\\d){1,10}"))) 
  
  return(outdf)
  
  
}


myPval2 <- function(obs,
                    stoch,
                    normalized = TRUE) {
  
  ############ get files #################
  
  
  # read files
  obstab <- read.table(obs, sep = "\t", header = TRUE)
  sttab <- read.table(stoch, sep = "\t", header = TRUE)
  
  ########### get scores ####################
  
  
  cat("\n showing NORMALIZED scores and p values \n ")
  
  # subset to get just normalized data
  obssub <- obstab %>% dplyr::select(ends_with("_toUnifAvg"))
  stsub <- sttab %>% dplyr::select(ends_with("_toUnifAvg"))
  
  outdf <- data.frame()
  ### for each metric, get scores
  for(i in 1:length(colnames(stsub))) {
    
    # get observed 
    obs.score <- obssub[,i]
    # get stochastic scores 
    st.sc <- length(which(stsub[,i] <= obs.score))
    
    ## get p values
    p <- st.sc / nrow(stsub)
    
    # print metrics
    cat("\n -------------------------------------------")
    cat("\n total number of stochastic comparisons: ", nrow(stsub))
    cat("\n observed ", colnames(obssub)[i], ": ", obs.score)
    cat("\n p value ",  colnames(obssub)[i], ": ", p)
    cat("\n -------------------------------------------")
    
    ## build results dataframe
    res <- data.frame(
      no.comps = nrow(stsub),
      metric = colnames(obssub)[i],
      obsscore = obs.score,
      pval = p,
      meanstoch = summary(stsub)[,i][4],
      minstoch = summary(stsub)[,i][1],
      maxstoch = summary(stsub)[,i][6]
    )
    
    outdf <- rbind(outdf, res)
    
  }
  
  
  ### return observed scores
  # remove weird text
  outdf <- outdf %>% 
    mutate(across(ends_with("stoch"), ~str_extract(., "(\\d)\\.(\\d){1,10}"))) 
  
  return(outdf)
  
  
}
