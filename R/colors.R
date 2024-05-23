## set colors for the species
# EVS 12/2023

# use viridis palette
library(viridis)

#show_col(viridis(4))

# get hex codes
cols <- c("#440154FF", "#31688EFF", "#35B779FF", "#FDE725FF")

# assign to names (shortened)
names(cols) <- c("G. beringei", "G. gorilla", "H. sapien", "P. troglodyte")


# assign to fullnames
fullname <- c("Gorilla beringei", "Gorilla gorilla", "Homo sapien", "Pan troglodyte")
fullname.col <- cols
names(fullname.col) <- fullname

# version with underscore for phylo trees
unders <- c("Gorilla_beringei", "Gorilla_gorilla", "Homo_sapien", "P_troglodytes_schweinfurthii")
unders.cols <- cols
names(unders.cols) <- unders

# version with pretty-fied names
italic <- c("*G. beringei*", "*G. gorilla*", "*H. sapien*", "*P. troglodyte*")
italic.cols <- cols
names(italic.cols) <- italic
