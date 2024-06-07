## generate figure s1 (QQ plots)
# EVS 1/2024

library(tidyverse)
library(gridExtra)
library(ggpubr)
library(qqman)
library(phyloseq)
library(microViz)
library(cowplot)



# set pathway to output files (from PLINK)
### these files are huge and are not storeable on Github; see final publication for GWAS summary statistics
path <- "private/updateITS_results_indep_snp/"

# read everything in to one dataframe (small enough to do this!)
myfiles <- list.files(path, pattern = ".adjusted$")
alldat <- data.frame()
for(i in 1:length(myfiles)) {
  cat("\n now reading file ", i, " of ", length(myfiles))
  
  fi <- read.table(paste0(path, myfiles[i]), sep = "\t", header = TRUE, comment.char = "") 
  fi$taxa <- myfiles[i]
  alldat <- rbind(alldat, fi)
}

#### plot in loop
mytax <- unique(alldat$taxa)
myplotlist <- c()

for(i in 1:length(mytax)) {
  
  # subset
  sub <- alldat %>% filter(taxa == mytax[i])
  
  # get observed and expected P values
  myp <- data.frame(
    observed = -log10(sort(sub$UNADJ)),
    expected = -log10(ppoints(nrow(sub)))
  )
  # get title
  mytitle <- str_remove_all(str_remove_all(mytax[i],
                                           "\\.glm\\.linear\\.adjusted"), "indep_onlysnp\\.")
  mytitleFancy <- if_else(str_detect(mytitle, "OTU"),
                          # if OTU;
                          paste0(str_split(mytitle, "_")[[1]][1], ": ",
                                 str_split(mytitle, "_")[[1]][3]),
                          paste0(str_split(mytitle, "_")[[1]][1], ": ",
                                 str_split(mytitle, "_")[[1]][2])
  )
  
  # plot with ggplot
  #ggplot(myp, aes(x = expected, y = observed)) +
  # geom_point(alpha = 0.8) +
  # geom_abline(color = "red", linetype = "dashed", alpha = 0.8) +
  # theme_pubr() +
  #labs(x = "-Log10(Expected)", y = "-Log10(Observed)", title = "test")
  # plot with ggpubr
  myplot <- ggscatter(myp, x = "expected", y = "observed", alpha = 0.7,
                      add = "reg.line", add.params = list(color = "red", linetype = "dashed", alpha = 0.8),
                      xlab = "-Log10(Expected)", ylab = "-Log10(Observed)",
                      title = mytitleFancy)
  
  # save to list
  myplotlist[[i]] <- myplot
  
  
}

### arrange and draw
gs <- ggarrange(plotlist = myplotlist, ncol = 6, nrow = 8)

# save
ggsave(plot = gs, dpi = 600, filename = "figures/updateITS_indepSNP_SNP_alltaxa_allQQs.png",
       width = 20, height = 24, units = "in")