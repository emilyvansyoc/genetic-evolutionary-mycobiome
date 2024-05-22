### plot and get summary stats of mwas sigs
# EVS 4/2024

library(tidyverse)
library(ggpubr)
library(phyloseq)
library(microViz)
library(vegan)
library(ggtext)


# get fungi
load("data/updated/phylo-ITS-filteredrarefied_nosubset.RData")

# get pathways
load("wgs-analyses/wgs-data/stool_pathways.RData")

### get data formatted the same way it was input to the lm

pathpa <- pathf %>% 
  filter(!str_detect(Pathway, "\\|")) %>% 
  mutate(pa = if_else(avgAbund > 0, 1, 0)) %>% 
  group_by(Pathway) %>% 
  summarize(tot = sum(pa)) %>% 
  filter(tot > (83*0.1))
pathfilt <- pathf %>% 
  filter(!str_detect(Pathway, "\\|")) %>% 
  filter(Pathway %in% pathpa$Pathway)
length(unique(pathfilt$Pathway)) # 341

# filter individuals to match shotgun data
psfilt <- psf %>% 
  ps_filter(sample_names(psf) %in% pathfilt$RANDSID)

# make fungi vertical and get genus level abundance data
psclr <- psfilt %>% 
  tax_fix() %>% 
  tax_glom("Genus") %>% 
  # filter for 10% prevalence
  tax_filter(min_prevalence = (83*0.1)) %>% # 26 genera
  tax_transform("clr", zero_replace = "halfmin")
taxa_names(psclr) <- psclr@tax_table[,"Genus"]
fungv <- psclr %>% 
  ps_melt() %>% 
  dplyr::select(Fungi = OTU, Sample, Abundance) %>% 
  mutate(RANDSID = as.numeric(Sample)) # 26 fungal genera

# join
moddf1 <- fungv %>% 
  full_join(pathfilt, relationship = "many-to-many")
# remove pathways with zero total abundance
compa <- moddf1 %>% 
  mutate(pa = if_else(avgAbund > 0, 1, 0)) %>% 
  group_by(Pathway) %>% 
  summarize(tot = sum(pa)) %>% 
  filter(tot > 0) %>% ungroup() # 361
moddf1 <- moddf1 %>% 
  filter(Pathway %in% compa$Pathway)

#### build plots for each fungi: Sarocladium
sar <- moddf1 %>% filter(Fungi == "Sarocladium") %>% 
  # get the sig pathways
  filter(str_detect(Pathway, "PWY-1861")) %>% 
  mutate(logGene = log10(avgAbund + 0.1)) %>% 
  # pretty-fy pathways for axis titles (5/1/2024)
  mutate(pname = "formaldehyde assimilation")


## Leotiomycetes class
leo <- moddf1 %>% filter(Fungi == "Leotiomycetes Class") %>% 
  filter(str_detect(Pathway, "PROPFERM-PWY") | str_detect(Pathway, "PWY-5494"))%>% 
  mutate(logGene = log10(avgAbund + 0.1)) %>% 
  mutate(pname = if_else(str_detect(Pathway, "PROPFERM"), "propionate fermentation from lactate", "proprionate fermentation from L-alanine"))

## Pichia
pic <- moddf1 %>% filter(Fungi == "Pichia") %>% 
  filter(str_detect(Pathway, "PWY-7198"))%>% 
  mutate(logGene = log10(avgAbund + 0.1)) %>% 
  mutate(pname = "pyrimidine deoxyribonucleic acid synthesis")

### plot individually and add together
sarplot <- ggscatter(sar, x = "logGene", y = "Abundance", 
                     add = "reg.line", conf.int = TRUE, add.params = list(color = "darkblue"),
                     size = 6, alpha = 0.7,
                     xlab = "log<sub>10</sub>(formaldehyde assimilation)", ylab = "*Sarocladium* abundance") +
  geom_richtext(x = 0, y = 5.5, label = "R<sup>2</sup>=21.3%, Q=0.025", label.size = 0, size = 6,
                hjust = "left") +
  theme_pubr(base_size = 18) +
  theme(axis.title.y = element_markdown(),
        axis.title.x = element_markdown())
leoplot1 <- ggscatter(leo %>% filter(str_detect(Pathway, "PROPFERM-PWY")), x = "logGene", y = "Abundance", 
                      add = "reg.line", conf.int = TRUE,add.params = list(color = "darkblue"), 
                      size = 6, alpha = 0.7,
                      xlab = "log<sub>10</sub>(proprionate fermentation from L-alanine)", ylab = "Leotiomycetes abundance")+
  geom_richtext(x = -0.5, y = 2, label = "R<sup>2</sup>=20.2%, Q=0.033", label.size = 0, size = 6,
                hjust = "left") +
  theme_pubr(base_size = 18) +
  theme(
    axis.title.x = element_markdown())
leoplot2 <- ggscatter(leo %>% filter(str_detect(Pathway, "PWY-5494")), x = "logGene", y = "Abundance", 
                      add = "reg.line", conf.int = TRUE, add.params = list(color = "darkblue"),
                      size = 6, alpha = 0.7,
                      xlab = "log<sub>10</sub>(proprionate fermentation from lactate)", ylab = "Leotiomycetes abundance")+
  geom_richtext(x = -0.5, y = 2, label = "R<sup>2</sup>=21.2%, Q=0.025", label.size = 0, size = 6,
                hjust = "left") +
  theme_pubr(base_size = 18) +
  theme(
    axis.title.x = element_markdown())
picplot <- ggscatter(pic, x = "logGene", y = "Abundance", 
                     add = "reg.line", conf.int = TRUE, add.params = list(color = "darkblue"),
                     size = 6, alpha = 0.7,
                     xlab = "log<sub>10</sub>(pyrimidine DNA biosynthesis)", ylab = "*Pichia* abundance")+
  geom_richtext(x = 0, y = 3, label = "R<sup>2</sup>=82.0%, Q<0.001", label.size = 0, size = 6,
                hjust = "left") +
  theme_pubr(base_size = 18) +
  theme(axis.title.y = element_markdown(),
        axis.title.x = element_markdown())
# add 
ggarrange(picplot, sarplot, leoplot1, leoplot2, labels = c("A.", "B.", "C.", "D."),
          font.label = list(size = 20, face = "bold"), vjust = 0) +
  theme(plot.margin = margin(t=2, r=0, b=0, l=0, unit = "cm"),
        plot.background = element_rect(fill = "white", color = NULL))
# save
ggsave(filename = "figures/mwas_paneled.png", dpi = 600, bg = "transparent")
