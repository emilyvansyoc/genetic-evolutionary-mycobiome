#!/bin/bash

## loop through fasta files in a directory and create an alignment for each one

source ~/miniconda3/etc/profile.d/conda.sh 
conda activate muscle
#set -uex

## set input directory
IN=$1
## set output directory for alignments
OUT=$2
mkdir -p $OUT


# loop for alignments in MAFFT
for file in "$IN"/*.fasta; do
    echo "WORKING ON: $(basename "$file")"
    output_file=$OUT/align_$(basename "$file")
    clipout=$OUT/clipkit_$(basename "$file")
    # check if file exists
    if [ -e "$output_file" ] & [ -e "$clipout" ]; then
        echo "SKIPPING: $output_file already exists"
        continue # skip to the next iteration
    fi

    # run alignment in mafft
    mafft --auto $file > $output_file

    # run clipkit
    clipkit $output_file -o $clipout
done

echo "#### done with alignments #####"
