# run topological congruency tests
# EVS 2/2024

library(tidyverse)
library(phyloseq)
library(microViz)
library(vegan)
library(ape)

# get phyloseq
load("private/hominid_phyloITS.RData")


#### ----- run TC: OTU level -----

# set path to TreeCMP
TC <- "dir/to/TreeCMp"

# set path to output files 
outpath <- "my-output"
if(!dir.exists(outpath)) {dir.create(outpath)}


# run TreeCMP on the random trees versus the host tree 
system(paste0("java -jar ", TC, "/treeCmp.jar -r ", 
              "data/hominid-tree.newick -i data/random_hominid.newick -d rc mc tt mt rf -N -o ", 
              outpath, "/host_v_random.txt"))

# rename to scientific name so tip labels match
psfilt <- psf %>% 
  tax_fix() %>% 
  #tax_glom("Genus") %>% 
  ps_mutate(SciName = case_when(
    SpeciesCaptive %in% "Wild_Chimp" ~ "P_troglodytes_schweinfurthii",
    SpeciesCaptive %in% "Wild_Lowland Gorilla" ~ "Gorilla_gorilla",
    SpeciesCaptive %in% "Wild_Mountain Gorilla" ~ "Gorilla_beringei",
    SpeciesCaptive %in% "Human_Bantu" ~ "Homo_sapien",
    SpeciesCaptive %in% "Human_BaAka" ~ "Homo_sapien"
  ))

## collapse to scientific name
#taxa_names(psfilt) <- psfilt@tax_table[,"Genus"]
df <- psfilt %>% ps_melt() 
# get mean and median
dfc <- df %>% 
  group_by(OTU, SciName) %>% 
  summarize(avg = round(mean(Abundance), 0),
            med = round(median(Abundance), 0))

# re-phyloseq
# median
psmed <- phyloseq(
  sample_data(data.frame(
    row.names = c("Gorilla_beringei", "Gorilla_gorilla", "Homo_sapien", "P_troglodytes_schweinfurthii"),
    SciName = c("Gorilla_beringei", "Gorilla_gorilla", "Homo_sapien", "P_troglodytes_schweinfurthii")
  )),
  otu_table(dfc %>% dplyr::select(-med) %>% 
              pivot_wider(names_from = OTU, values_from = avg) %>% 
              column_to_rownames(var = "SciName"), taxa_are_rows = FALSE),
  tax_table(df %>% 
              dplyr::select(OTU, Kingdom, Phylum, Class, Order, Family, Genus) %>% 
              distinct() %>% 
              column_to_rownames(var = "OTU") %>% 
              as.matrix())
)

## jaccard
write.tree(unroot(psmed %>% dist_calc("jaccard", binary = TRUE) %>% dist_get() %>% hclust(method = "average") %>% as.phylo()), file = paste0(outpath, "/fung_jaccard.newick"))

## bray
write.tree(unroot(psmed %>% dist_calc("bray") %>% dist_get() %>% hclust(method = "average") %>% as.phylo()),
           file = paste0(outpath, "/fung_bray.newick"))

### run TreeCMP in a loop
dists <- c("jaccard", "bray")

for(i in 1:length(dists)) {
  
  cat("\n now running distance metric ", dists[i])
  
  # run
  system(paste0("java -jar ", TC, "/treeCmp.jar -r ", 
                "data/hominid-tree.newick -i ",  outpath, 
                "/fung_", dists[i], ".newick -d rc mc  -N -o ",
                outpath, "/obs_", dists[i], ".txt"))
}
#### get results
outpath <- "my-dir"
source("R/fx_myPvalue.R")
# jaccard
jac <- myPval2(obs = paste0(outpath, "/obs_jaccard.txt"), stoch = paste0(outpath, "/host_v_random.txt"),
               normalized = TRUE) # no sig
# bray-curtis
bray <- myPval2(obs = paste0(outpath, "/obs_bray.txt"), stoch = paste0(outpath, "/host_v_random.txt"),
                normalized = TRUE) # no sig


## ----- run TC: genus level ----

# set path to TreeCMP
TC <- "dir/to/TreeCMP"

# set path to output files 
outpath <- "my-dir"
if(!dir.exists(outpath)) {dir.create(outpath)}


# run TreeCMP on the random trees versus the host tree 
system(paste0("java -jar ", TC, "/treeCmp.jar -r ", 
              "data/hominid-tree.newick -i data/random_hominid.newick -d rc mc tt mt rf -N -o ", 
              outpath, "/host_v_random.txt"))

# rename to scientific name so tip labels match
psfiltg <- psf %>% 
  tax_fix() %>% 
  tax_glom("Genus") %>% 
  ps_mutate(SciName = case_when(
    SpeciesCaptive %in% "Wild_Chimp" ~ "P_troglodytes_schweinfurthii",
    SpeciesCaptive %in% "Wild_Lowland Gorilla" ~ "Gorilla_gorilla",
    SpeciesCaptive %in% "Wild_Mountain Gorilla" ~ "Gorilla_beringei",
    SpeciesCaptive %in% "Human_Bantu" ~ "Homo_sapien",
    SpeciesCaptive %in% "Human_BaAka" ~ "Homo_sapien"
  ))
## collapse to scientific name
taxa_names(psfiltg) <- psfiltg@tax_table[,"Genus"]
df <- psfiltg %>% ps_melt() 
# get mean and median
dfc <- df %>% 
  group_by(OTU, SciName) %>% 
  summarize(avg = round(mean(Abundance), 0),
            med = round(median(Abundance), 0))

# re-phyloseq
# median
psmed <- phyloseq(
  sample_data(data.frame(
    row.names = c("Gorilla_beringei", "Gorilla_gorilla", "Homo_sapien", "P_troglodytes_schweinfurthii"),
    SciName = c("Gorilla_beringei", "Gorilla_gorilla", "Homo_sapien", "P_troglodytes_schweinfurthii")
  )),
  otu_table(dfc %>% dplyr::select(-med) %>% 
              pivot_wider(names_from = OTU, values_from = avg) %>% 
              column_to_rownames(var = "SciName"), taxa_are_rows = FALSE),
  tax_table(df %>% 
              dplyr::select(OTU, Kingdom, Phylum, Class, Order, Family, Genus) %>% 
              distinct() %>% 
              column_to_rownames(var = "OTU") %>% 
              as.matrix())
)

## jaccard
write.tree(unroot(psmed %>% dist_calc("jaccard", binary = TRUE) %>% dist_get() %>% hclust(method = "average") %>% as.phylo()), file = paste0(outpath, "/fung_jaccard.newick"))

## bray
write.tree(unroot(psmed %>% dist_calc("bray") %>% dist_get() %>% hclust(method = "average") %>% as.phylo()),
           file = paste0(outpath, "/fung_bray.newick"))


### run TreeCMP in a loop
dists <- c("jaccard", "bray")

for(i in 1:length(dists)) {
  
  cat("\n now running distance metric ", dists[i])
  
  # run
  system(paste0("java -jar ", TC, "/treeCmp.jar -r ", 
                "data/hominid-tree.newick -i ",  outpath, 
                "/fung_", dists[i], ".newick -d rc mc  -N -o ",
                outpath, "/obs_", dists[i], ".txt"))
}
#### get results
source("R/fx_myPval.R")
outpath <- "my-dir"
# jaccard
jac <- myPval2(obs = paste0(outpath, "/obs_jaccard.txt"), stoch = paste0(outpath, "/host_v_random.txt"),
               normalized = TRUE) 
# bray-curtis
bray <- myPval2(obs = paste0(outpath, "/obs_bray.txt"), stoch = paste0(outpath, "/host_v_random.txt"),
                normalized = TRUE) 


