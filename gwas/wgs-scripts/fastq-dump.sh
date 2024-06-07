#!/bin/bash

### submitted as a job

## run fastq-dump on raw reads

source ~/.bashrc
eval "$(micromamba shell hook --shell bash)"
micromamba activate sratools
set -uex

### set working directory
DIR=/my-pathway/

## name of IDs file
IDS=$DIR/runids_hmp48479.txt

## output directory for raw reads
RAW=~/scratch/raw/
mkdir -p $RAW

## ----- get reads from SRA with fastq-dump ----

# paired-end reads
cat $DIR/norun | parallel -j $SLURM_CPUS_PER_TASK fastq-dump {} --split-files --outdir $RAW

# single-end reads
cat $SIDS | parallel -j $SLURM_CPUS_PER_TASK fastq-dump {} --outdir $RAW