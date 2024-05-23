#!/bin/bash

## dereplicate exact sequence variants by sequence name and identity

source ~/.bashrc
eval "$(micromamba shell hook --shell bash)"
micromamba activate vsearch
set -uex

## directory to input files
IN=$1

## directory to output files
OUT=$2
mkdir -p $OUT

## run each sample
for file in "$IN"/*.fasta; do  
    echo "WORKING ON: $(basename "$file")"
    output_file=$OUT/derep_$(basename "$file")
    # run vsearch
    vsearch --derep_id $file --output $output_file

done


