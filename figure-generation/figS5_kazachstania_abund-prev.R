
# get kazachstania abundance
library(microViz)
library(phyloseq)
load("data/updated/phylo_ITS_resolvedNA.RData")
psf <- psnewname
psg <- psf %>% tax_fix() %>% tax_glom("Genus")
taxa_names(psg) <- psg@tax_table[,"Genus"]

# get only the top 25 taxa
topn <- tax_top(psg, n = 25, rank = "Genus")

# get relative abundance of all genera
ra <- psg %>% 
  tax_select(topn, "Genus", strict_matches = TRUE) %>% 
  tax_transform("compositional") %>% 
  ps_melt() %>% dplyr::select(Sample, Abundance, OTU) %>% 
  mutate(OTU = factor(OTU, ordered = TRUE, levels = rev(topn))) %>% 
  mutate(iskaz = if_else(OTU == "Kazachstania", "kaz", "not"))



# get prevalence of all genera
kzp <- psg %>% 
  tax_select(topn, "Genus", strict_matches = TRUE) %>% 
  tax_transform("pa") %>% ps_melt() %>%
  dplyr::select(Sample, Pres = Abundance, OTU) %>% 
  # calculate total prevalence
  group_by(OTU) %>% summarize(tot = sum(Pres)) %>% mutate(prev = (tot/125)*100) %>% 
  ungroup() %>% 
  mutate(OTU = factor(OTU, ordered = TRUE, levels = rev(topn))) %>% 
  mutate(iskaz = if_else(OTU == "Kazachstania", "kaz", "not"))

## plot
# boxplot for RA
## super hacky way to change Kazachsatnia color to red
labs <- paste("<span style = 'color: ",
              ifelse(topn == "Kazachstania", "red", "black"),
              ";'>",
              topn,
              "</span>", sep = "")
ra <- ra %>% 
  mutate(mylab = paste("<span style = 'color: ",
                       ifelse(iskaz == "kaz", "red", "black"),
                       ";'>",
                       OTU,
                       "</span>", sep = "")) %>% 
  mutate(mylab = factor(mylab, ordered = TRUE, levels = rev(labs)))
pa <- ggplot(data = ra, aes(x = Abundance, y = mylab, color = iskaz)) +
  geom_boxplot() + geom_jitter(height = 0.1) +
  theme_pubr(base_size = 16) +
  xlab("Relative Abundance") + ylab("") +
  theme(axis.text.y = element_markdown(face = "italic"),
        legend.position = "none") +
  scale_color_manual(values = c("red", "black"))

# barplot for prevalence
pb <- ggplot(data = kzp, aes(x = prev, y = OTU, fill = iskaz)) +
  geom_bar(stat = "identity") +
  #scale_x_reverse(limits = c(100, 0)) +
  xlim(c(0, 100)) +
  theme_pubr(base_size = 16) +
  xlab("Prevalence (%)") + ylab("") +
  scale_fill_manual(values = c("red", "lightgrey")) +
  #guides(y = "none", y.sec = "axis") + # move axis line to the right side
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

## arrange
ggarrange(pa,pb, ncol = 2, widths = c(1, 0.5))
ggsave(filename = "figures/kaz_ra_prev.png", dpi = 600)
