## generate figure 3 and Figure S8


library(tidyverse)
library(ape)
library(vegan)
library(ggtree)
library(ggpubr)
library(ggtext)
library(scales)
library(treeio)


## get OTU trees 
load("data/all_otu_subtrees_04-10-2024_nonegbl.RData")

# get sigs
load("data/final_sigs_otus_parafitpaco.RData")

# get colors
source("R/colors.R")


## ---- make function for making circular trees ---

myPlot <- function(mytree) {
  
  # break down hominid representation
  tipdf <- data.frame(tips = mytree$tip.label) %>% 
    mutate(tip1 = str_remove(tips, "_(\\d){1,4}$")) %>% 
    mutate(hominid = str_extract(tip1, "[:alpha:]{1,12}$")) %>% 
    column_to_rownames(var = "tips")
  tipdf %>% group_by(hominid) %>% count()
  
  # add sci name to tipdf
  tipdf <- tipdf %>% 
    mutate(italname = case_when(
      hominid == "Chimp" ~ "*P. troglodyte*",
      hominid == "Human" ~ "*H. sapien*",
      hominid == "LowGorilla" ~ "*G. gorilla*",
      hominid == "MountGorilla" ~ "*G. beringei*"
    )) 
  
  # simplify tip labels
  plottree <- mytree
  plottree$tip.label <- tipdf$italname
  td <- data.frame(rownames = plottree$tip.label,
                   hominid = plottree$tip.label)
  
  # plot circular
  tplot <- plottree %>% ggtree(ladderize = FALSE, layout = "fan", branch.length = "none") %<+% td + 
    geom_tippoint(aes(color = hominid), size = 2) +
    scale_color_manual(values = italic.cols) +
    theme(#legend.text = element_markdown(size = 18),
      #legend.title = element_blank(),
      legend.position = "none",
      plot.margin = margin(t=0, r=0, b=0, l=0, unit = "pt"),
      plot.background = element_rect(fill = "transparent", color = NA),
      panel.background = element_rect(fill = "transparent", color = NA))
  
  # return plot
  return(tplot)
  
}

## ---- plot all 8 significant trees ---

# order by decreasing r2
sigps <- sigps %>% arrange(desc(pacor2))

sigtrees <- alltrees[names(alltrees) %in% sigps$otu]
# do by hand to get ordering right for now
sigtrees <- list(alltrees$OTU_144, alltrees$OTU_299, alltrees$OTU_405, alltrees$OTU_192, alltrees$OTU_242, alltrees$OTU_499, alltrees$OTU_4, alltrees$OTU_224, alltrees$OTU_79, alltrees$OTU_6)
names(sigtrees) <- c("OTU_144", "OTU_299", "OTU_405", "OTU_192", "OTU_242", "OTU_499", "OTU_4", "OTU_224", "OTU_79", "OTU_6")

plist <- lapply(sigtrees, myPlot)
names(plist) <- names(sigtrees)

## 5/3; save individually
for(i in 1:length(plist)){
  ggsave(plist[[i]], filename = paste0("updated_methods/fig_images/codiv_trees/otutree_", names(sigtrees)[i], ".png"), dpi = 600,
         height = 18, width = 30, units = "in")
}

# ---- plot the 2 that are significant but weak (supplementary figure 8) ----

trees <- list(alltrees$OTU_79, alltrees$OTU_6)
names(trees) <- c("OTU_79", "OTU_6")
plist <- lapply(trees, myPlot)

# save together and make labels in Illustrator to color-code to main text figure
ggarrange(plotlist = plist, ncol = 2, nrow = 1)
ggsave(filename = "updated_methods/fig_images/codiv_supplementary.png", dpi = 600)


### ----- get full fungi tree and make schematic ----

ftree <- read.tree("data/li2021_fungaltree/1672taxa_290genes_bb_1.treefile")


# get tip labels
ftips <- data.frame(tiplab = ftree$tip.label) %>% 
  mutate(genus = sapply(str_split(tiplab, "_"), `[`, 1))


# get the genera present in our dataset
load("data/all_genera_list.RData")

# subset
fsub <- ftips %>% 
  filter(genus %in% gen)


#### get Pleosporales nodes
# BLAST suggest these may belong to the Thyridariaceae family
tab <- readxl::read_xlsx(path = "data/li2021_fungaltree/1-s2.0-S0960982221001391-mmc3.xlsx", sheet = "B")
# get Pleosporale
pleo <- tab %>% filter(order == "Pleosporales")
pleosub <- pleo %>% filter(tip_id %in% ftree$tip_id)
pleotree <- keep.tip(ftree, tip = pleo$tip_id)

# make dataframe of node labels by hand
# add labels to nodes
mynodes <- data.frame(
  genus = c("Aureobasidium", "Saturnispora", "Malassezia", "Talaromyces", "Geotrichum", "Xylaria/Nigrospora", "Pleosporales", "Cladosporium"),### add the two that were missing (Pleosporales and Nigrospora) -> Xylaria is now duplicated
  mynode = c("node475", "node247", "node414", "node589", "node327", "node159", "node449", "node467")) %>%  
  left_join(td, by = c("mynode" = "label"))

# get the default colors
cols <- hue_pal()(8)
names(cols) <- mynodes$genus


### build plot (text is added in Illustrator)
tr %>% ggtree() +
  geom_hilight(data = mynodes, aes(node = node, fill = genus),extend = 3.25, alpha = 0.9,  type = "gradient", gradient.direction = "tr") +
  xlim(0, 3.25) +
  rotate() +
  theme(legend.position = "none") +
  scale_fill_manual(values = cols)
# save
ggsave(filename = "updated_methods/fig_images/ftree_finalotus.png", dpi = 600)


