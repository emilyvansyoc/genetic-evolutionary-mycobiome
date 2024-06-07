### generate Figure 2 
#EVS 
# note; there is a bug in ggtree that causes super annoying issues when phyloseq is also loaded into namespace, but OK to ignore 


library(ggpubr)
library(microViz)
library(ape)
library(tidyverse)
library(phyloseq)
library(ggtree)
library(gridExtra)
library(rphylopic)
library(ggimage)
library(cowplot) # for themes 
library(ggtext) # for element markdown
library(gridExtra) # for grid.arrange
library(multcompView) # for letters in boxplot
# to get pairwise distances into a nice dataframe:
#devtools::install_github("kylebittinger/usedist")
library(usedist)
library(ggh4x) # color strip backgrounds

# get colors
source("R/colors.R")


#### ---- Panel A: topological congruency ----

## load hominid dendrogram
load("data/hominidtree_forplot.RData")

# get fungal dendrogram
load("data/fung_bray_dendrogram.RData")
host <- rt

# get phylopic images
h.img <- get_phylopic(uuid = "036b96de-4bca-408e-adb6-2154fcd724ef", preview = TRUE)
g.img <- get_phylopic(uuid = "142e0571-3b5f-443d-a887-b572a224ea22", preview = TRUE)
c.img <- get_phylopic(uuid = "7133ab33-cc79-4d7c-9656-48717359abb4", preview = TRUE)
f.img <- get_phylopic(uuid = "aaecd181-feb8-4203-8c64-f46384257e59", preview = TRUE) # hyphae
f.img1 <- get_phylopic(uuid = "e602729e-044f-4d2b-bab2-64a87a0b48c7", preview = TRUE) # yeast

# build dataframe for plotting in ggtree
labdf <- data.frame(
  row.names = host$tip.label,
  name = host$tip.label,
  uid = c("142e0571-3b5f-443d-a887-b572a224ea22", "142e0571-3b5f-443d-a887-b572a224ea22", ## GORILLA IS TWICE
          "036b96de-4bca-408e-adb6-2154fcd724ef", "7133ab33-cc79-4d7c-9656-48717359abb4"),
  forcolor = host$tip.label
)

# 5/8/2024; build dataframe for fungi 
flabdf <- data.frame(
  row.names = host$tip.label,
  name = host$tip.label,
  hyphae = c("aaecd181-feb8-4203-8c64-f46384257e59"),
  yeast = c("e602729e-044f-4d2b-bab2-64a87a0b48c7"),
  forcolor = host$tip.label
)

# make the hominid plot
plot_hostup <- ggtree(host, size = 6, branch.length = "none"
) %<+% labdf + ## the weird symbol is ggtree "attacher"
  geom_tiplab(aes(image = uid, color = forcolor), geom = "phylopic", 
              offset = 0.3, size = c(0.11, # human 
                                     0.16, 0.16, # gorillas
                                     0.13 # chimp
              )) +
  xlim(NA,6.1) +
  scale_color_manual(values = unders.cols) +
  theme(legend.position = "none",
        plot.margin = margin(r=0))

# make the fungal plot; layer geom_tiplab over each other to create the image we want
plot_fungid <- uf %>% ggtree(size = 5, branch.length = "none")%<+% #+ 
  #geom_tippoint(aes(color = label), size = 15) +
  flabdf + 
  geom_tiplab(aes(image = hyphae, color = forcolor), geom = "phylopic",
              offset = -0.4, size = 0.16, alpha = 0.7)  +
  geom_tiplab(aes(image = yeast, color = forcolor), geom = "phylopic",
              offset = -0.4, size = 0.16) +
  scale_x_reverse(limits = c(4, 0)) +
  # manual colors
  scale_color_manual(values = unders.cols) +
  theme(legend.position = "none",
        plot.margin = margin(l=0)) 

## arrange together
fin.plot <- grid.arrange(plot_hostup, plot_fungid, ncol = 2, widths = c(1, 0.8))

# save plot
ggsave(filename = "figures/topocong_phylopic_v1.png", plot = fin.plot,
       dpi = 600)

### ---- panel B: iris plot ----

# get fungal data (restricted)
load("private/hominid_phyloITS.RData")

# genus level
psfilt <- psf %>% 
  tax_fix() %>% 
  tax_glom("Genus") 

## permanova to check stats
psfilt %>% tax_fix() %>% dist_calc("bray") %>% dist_permanova(variables = "SciName") # significant
# pairwise
library(pairwiseAdonis)
dis <- psfilt %>% tax_fix() %>% dist_calc("bray") %>% microViz::dist_get()
df <- psfilt %>% samdat_tbl()
comp <- pairwise.adonis2(dis ~ SciName, data = df) # all are significant
# check beta dispersion
psfilt %>% tax_fix() %>% dist_calc("bray") %>% dist_bdisp(variables = "SciName", method = "centroid") %>% 
  bdisp_get() # H sapien an G berinei; H sapien and G gorilla; H sapien and P troglodytes
b <- betadisper(dis, group = df$SciName, type = "centroid")
bmod <- anova(b)
plot(b)


# calculate bray
uf <- psfilt %>% tax_fix() %>% dist_calc("bray")

# pretty-fy names
uf <- uf %>% 
  ps_mutate(ItalName = case_when(
    SciName %in% "Gorilla_beringei" ~ "*G. beringei*",
    SciName %in% "Gorilla_gorilla" ~ "*G. gorilla*",
    SciName %in% "Homo_sapien" ~ "*H. sapien*",
    SciName %in% "P_troglodytes_schweinfurthii" ~ "*P. troglodyte*"
  ))

# make iris plot
uf %>% 
  ord_calc() %>% 
  ord_plot_iris(tax_level = "Family", anno_colour = "ItalName", 
                n_taxa = 18,
                anno_colour_style = list(linewidth = 7),
                taxon_renamer = function(tax) {
                  str_replace(tax, "_fam_Incertae_sedis",
                              " (IS)")}) +
  scale_color_manual(values = italic.cols, name = "", guide = NULL) +
  theme(legend.text = element_text(size = 18))

# save
ggsave(filename = "figures/irisplot_BC-genus.png", dpi = 600)

## ---- panel C: within- and between-species ----

# use 'usedist' package to convert to nice dataframe
ufdf <- dist_groups(d = uf, g = psfilt@sam_data$SciName)

# various wrangling to get it into plot format
dat <- ufdf %>% 
  filter(str_detect(Label, "Within")) %>% 
  mutate(newLab = Group1,
         is.within = "Within") %>% 
  dplyr::select(Item1, Item2, is.within, Distance, newLab) %>% 
  rbind(ufdf %>% 
          filter(str_detect(Label, "Between")) %>% 
          mutate(is.within = "Between") %>% 
          dplyr::select(Item1, Item2, Group1, Group2, is.within, Distance) %>% 
          pivot_longer(cols = c(Group1, Group2), names_to = "oldGroup", values_to = "newLab") %>% dplyr::select(-oldGroup)) %>% 
  
  mutate(newLab = factor(newLab, ordered = TRUE, levels = c("H. sapiens", "P. troglodytes", "G. gorilla", "G. beringei"
  )))%>% 
  mutate(is.within = factor(is.within, ordered = TRUE, levels = c("Between", "Within"))) %>% 
  mutate(fulllab = interaction(is.within, newLab, sep = "-", lex.order = TRUE)) %>% 
  ### 5/23; order by phylogeny not median
  mutate(newLab = factor(newLab, ordered = TRUE, levels = c("P. troglodytes", "H. sapiens", "G. gorilla", "G. beringei")))

# stats test
sampdf <- psfilt %>% samdat_tbl()
library(pairwiseAdonis)
pairwise.adonis2(uf ~ SciName, data = sampdf)
mod <- adonis2(uf ~ SciName, data = sampdf)
bd <- betadisper(uf, group = sampdf$SciName)
anova(bd)
TukeyHSD(bd, which = "group")

# confirm with pairwise tests
statv <- dat %>% 
  group_by(newLab) %>% 
  wilcox_test(Distance ~ is.within) %>% 
  add_significance() %>% 
  add_xy_position(x = "is.within") %>% 
  mutate(y.position = 0.1)

# make color vector (hacky)
sci.cols <- italic.cols
names(sci.cols) <- c("Within-G. beringei", "Within-G. gorilla", "Within-H. sapiens", "Within-P. troglodytes")
# add "between" colors (all white)
f.cols <- c(sci.cols, rep("white", 4))
names(f.cols)[5:8] <- c("Between-H. sapiens", "Between-P. troglodytes", "Between-G. gorilla", "Between-G. beringei")

# make plot
p1 <- ggplot(data = dat, aes(x = is.within, y = Distance)) +
  geom_jitter(
    width = 0.2, color = "darkgrey", alpha = 0.3) +
  geom_boxplot(#aes(fill = is.within),
    aes(fill = fulllab),
    size = 1.5, outlier.shape = NA, notch = FALSE, alpha = 0.7) +
  coord_flip() +
  ggh4x::facet_wrap2(~newLab, strip = strip_themed(background_y = elem_list_rect(fill = c("#FDE725FF", "#35B779FF", "#31688EFF", "#440154FF"))), 
                     ncol = 1, strip.position = "left") +
  labs(y = "Bray-Curtis pairwise distances", x = "") +
  stat_pvalue_manual(data = stats, label = "p.signif", hide.ns = TRUE, tip.length = 0, size = 10, coord.flip = TRUE) + # brackets
  scale_fill_manual(values = f.cols) +
  #ylim(c(0, 1.15)) +
  scale_x_discrete(position = "top") +
  theme_pubr(base_size = 22) +
  theme(strip.background = element_rect(fill = c("#FDE725FF", "#35B779FF", "#31688EFF", "#440154FF", alpha = 0.9), color = NA),
        strip.text = element_text(size = 22, face = "italic", color = "white"),
        axis.text.y = element_text(size = 22),
        legend.position = "none",
        panel.spacing = unit(0.05, "cm")) 

## ---- panel D: between each species ----

## pretty-fy names
forplot <- ufdf %>% filter(!str_detect(Label, "Within")) %>% 
  mutate(labplot = case_when(
    Label == "Between G. gorilla and P. troglodytes" ~ "Between *P. troglodytes* and *G. gorilla*",
    Label == "Between H. sapiens and P. troglodytes" ~ "Between *H. sapiens* and *P. troglodytes*",
    Label == "Between G. beringei and P. troglodytes" ~ "Between *P. troglodytes* and *G. beringei*",
    Label == "Between G. gorilla and H. sapiens" ~ "Between *H. sapiens* and *G. gorilla*",
    Label == "Between G. beringei and G. gorilla" ~ "Between *G. gorilla* and *G. beringei*",
    Label == "Between G. beringei and H. sapiens" ~ "Between *H. sapiens* and *G. beringei*"
  ))

#### THIS RELIES ON ANOVA TEST FOR CODING PURPOSES BUT RESULTS ARE IDENTICAL WITH K-W AND DUNN TEST
mod <- aov(sim ~ labplot, data = forplot)
t <- TukeyHSD(mod)

# use multcompview to get letters easily
cld <- multcompLetters4(mod, t)
cldf <- as.data.frame.list(cld$labplot) %>% 
  rownames_to_column(var = "labplot") %>% 
  dplyr::select(labplot, Letters) %>% 
  # get positions
  full_join(forplot %>% group_by(labplot) %>% get_summary_stats(sim, type = "common") %>% dplyr::select(median, labplot)) 

# get text colors
labs <- c("Between <b style='color:#FDE725FF'>*P. troglodytes*</b> and <b style='color:#31688EFF'>*G. gorilla*</b>",
          "Between <b style='color:#FDE725FF'>*P. troglodytes*</b> and <b style='color:#440154FF'>*G. beringei*</b>",
          "Between <b style='color:#31688EFF'>*G. gorilla*</b> and <b style='color:#440154FF'>*G. beringei*</b>",
          "Between <b style='color:#35B779FF'>*H. sapiens*</b> and <b style='color:#FDE725FF'>*P. troglodytes*</b>",
          "Between <b style='color:#35B779FF'>*H. sapiens*</b> and <b style='color:#31688EFF'>*G. gorilla*</b>",
          "Between <b style='color:#35B779FF'>*H. sapiens*</b> and <b style='color:#440154FF'>*G. beringei*</b>")

# add text colors
fp <- forplot %>% 
  mutate(fulllab = case_when(
    labplot == "Between *P. troglodytes* and *G. gorilla*" ~ labs[1],
    labplot == "Between *P. troglodytes* and *G. beringei*" ~ labs[2],
    labplot == "Between *G. gorilla* and *G. beringei*" ~ labs[3],
    labplot == "Between *H. sapiens* and *P. troglodytes*" ~ labs[4],
    labplot == "Between *H. sapiens* and *G. gorilla*" ~ labs[5],
    labplot == "Between *H. sapiens* and *G. beringei*" ~ labs[6]
  ))


### make plot
p2 <- ggplot(data = fp, aes(x = fct_reorder(fulllab, desc(Distance)), y = Distance)) +
  
  geom_point(position = position_jitter(width = 0.2, seed = 123), alpha = 0.5, color = "darkgrey") +
  geom_boxplot(size = 1.5, outlier.shape = NA, notch = FALSE, fill = "lightcyan4", alpha = 0.7) +
  geom_text(data = cldf1, aes(x = fulllab, label = Letters, y = 0.05),
            size = 10) +
  #ylim(-0.1, 1.1) +
  scale_y_continuous(limits = c(-0.1, 1.1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  coord_flip() +
  theme_pubr(base_size = 24) +
  labs(x = "", y = "Bray-Curtis pairwise distances") +
  theme(axis.text.y = element_markdown())

### arrange both and save
ggarrange(p1, p2, ncol = 2, widths = c(0.8, 1))
ggsave("figures/pairwise_both.png", dpi = 600,
       height = 9, width = 25, units = "in")

