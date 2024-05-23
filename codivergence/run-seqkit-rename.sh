#!/bin/bash

### add unique identifier to alignment or fasta files to run phylogeny

# set input directory 
IN=$1
echo "####### looking for "clipkit" prefix, edit bash file if different #########"

# set output directory
OUT=$2
mkdir -p $OUT

## need seqkit from micromamba environment
source ~/.bashrc
eval "$(micromamba shell hook --shell bash)"
micromamba activate iqtree

# run loop 
for file in "$IN"/clipkit_*; do

    seqkit rename -n $file -o $OUT/unique_$(basename "$file")
done