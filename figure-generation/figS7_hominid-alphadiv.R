## alpha div and phylo relationships

# EVS 4/2024

library(tidyverse)
library(phyloseq)
library(microViz)
library(vegan)
library(ape)
library(ggpubr)
library(rstatix)
library(ggtext)

# get phyloseq
load("updated_methods/phyloITS.RData")

# rename to scientific name so tip labels match
psfilt <- psf %>% 
  tax_fix() %>% 
  tax_glom("Genus") %>% 
  ps_mutate(SciName = case_when(
    SpeciesCaptive %in% "Wild_Chimp" ~ "*P. troglodytes*",
    SpeciesCaptive %in% "Wild_Lowland Gorilla" ~ "*G. gorilla*",
    SpeciesCaptive %in% "Wild_Mountain Gorilla" ~ "*G. beringei*",
    SpeciesCaptive %in% "Human_Bantu" ~ "*H. sapien*",
    SpeciesCaptive %in% "Human_BaAka" ~ "*H. sapien*"
  ))

# get adiv
adiv <- estimate_richness(psfilt, measures = c("Shannon", "Simpson")) %>% 
  rownames_to_column(var = ".sample_name") %>% 
  full_join(psfilt %>% samdat_tbl())

# compare
mod <- aov(Shannon ~ SciName, data = adiv)
summary(mod)
plot(resid(mod))
qqnorm(resid(mod))
TukeyHSD(mod) # Gorilla-gorilla, homo-gorilla, gorilla b - chimp: no diff between homo-gorilla g, homo-pan, or pan-gorilla gorilla 
mod <- aov(Simpson ~ SciName, data = adiv)
summary(mod)
TukeyHSD(mod) # same relationships


#### ---- make adiv working figure for supplementary ----

# get stats
stats <- adiv %>% 
  tukey_hsd(formula = Shannon ~ SciName) %>% 
  add_xy_position(step.increase = 0.2) %>% 
  mutate(lab = if_else(p.adj < 0.05, "*", "ns"))

# order by increasing means
adiv %>% group_by(SciName) %>% get_summary_stats(Shannon, type = "median_iqr") %>% 
  arrange(desc(median))
adiv <- adiv %>% mutate(SciName = factor(SciName, ordered = TRUE, levels = c("*P. troglodytes*", "*G. gorilla*", "*H. sapien*", "*G. beringei*")))

# plot
p1 <- ggboxplot(adiv, x = "SciName", y = "Shannon",
                add = "jitter", add.params = list(color = "darkgrey", size = 4, alpha = 0.8),
                xlab = "", ylab = "Shannon's H",
                size = 1.5) +
  stat_pvalue_manual(stats, label = "lab", hide.ns = TRUE, tip.length = 0, size = 9) +
  theme_pubr(base_size = 18) +
  theme(axis.text.x = element_markdown())

ggsave(filename = "updated_methods/fig_images/adiv_hominid.png", dpi = 600)

### ---- test for homogeneity of variance ----

library(car)
lt <- leveneTest(Shannon ~ SciName, data = adiv, center = median) # significant
TukeyHSD(lt)

# calculate medians
meds <- aggregate(Shannon ~ SciName, data = adiv, FUN = mean)
colnames(meds)[2] <- "sh_mean"
adiv1 <- adiv %>% full_join(meds)
adiv1$res <- abs(adiv1$Shannon - adiv1$sh_mean)
res.aov <- aov(res ~ SciName, data = adiv1) # same p, F, and DF as leveneTest
TukeyHSD(res.aov) # diff between human-chimp, human-gorilla, marginal mountain-chimp

# plot
s1 <- adiv1 %>% 
  tukey_hsd(formula = res ~ SciName) %>% 
  add_xy_position(step.increase = 0.2) %>% 
  mutate(lab = if_else(p.adj < 0.05, "*", "ns"))

p2 <- ggboxplot(adiv1, x = "SciName", y = "res",
                add = "jitter", add.params = list(color = "darkblue", size = 4, alpha = 0.8),
                xlab = "", ylab = "Residuals",
                size = 1.5) +
  stat_pvalue_manual(s1, label = "lab", hide.ns = TRUE, tip.length = 0, size = 9) +
  theme_pubr(base_size = 18) +
  theme(axis.text.x = element_markdown())

## combine
ggarrange(p1, p2, ncol = 1, labels = c("A.", "B."), font.label = list(size = 20, face = "bold"))
ggsave(filename = "updated_methods/fig_images/adiv_paneled.png", dpi = 600, height = 10, width =10, units = "in")
