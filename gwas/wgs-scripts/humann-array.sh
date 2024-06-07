#!/bin/bash

## this is run as an array on an HPC

### run humann3
source ~/.bashrc
conda activate humenv
set -uex

# set up array files
# input directory (CANNOT take paired end reads)
DATADIR=~/scratch/tmp_dup_hiseq/
arrayoffiles=($(ls $DATADIR | grep sort.trimmed.1.fastq | cut -f 1 -d . | sort | uniq))

# output directory; make in scratch and move over since temp files are often not removed
OUT=~/scratch/humann-out_hiseq-nodup/
mkdir -p $OUT

#### print name of file to stout ####
echo "########################## WORKING ON SAMPLE"
echo $DATADIR/${arrayoffiles[${SLURM_ARRAY_TASK_ID}]}
echo "##########################"

## concatenate all fastq files together
cat $DATADIR/${arrayoffiles[${SLURM_ARRAY_TASK_ID}]}*.fastq > $DATADIR/${arrayoffiles[${SLURM_ARRAY_TASK_ID}]}_all.fastq


## run humann3
humann --input $DATADIR/${arrayoffiles[${SLURM_ARRAY_TASK_ID}]}_all.fastq --output $OUT --threads $SLURM_CPUS_PER_TASK --input-format fastq --nucleotide-database /refs/human-choco/chocophlan/ --protein-database /refs/uniref/ --metaphlan-options "--bowtie2db /refs/choco/"

## check for the existence of output files
if [ -e "$OUT/${arrayoffiles[${SLURM_ARRAY_TASK_ID}]}_all_pathabundance.tsv" ]; then
    echo "Humann ran successfully and the output file ${arrayoffiles[${SLURM_ARRAY_TASK_ID}]} pathabundance exists"
    else
        # print message
        echo "Humann did not run successfully: the pathabundance file does not exist"
fi