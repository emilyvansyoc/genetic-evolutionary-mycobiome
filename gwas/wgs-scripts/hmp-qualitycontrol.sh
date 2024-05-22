#!/bin/bash

### submitted as a job array to HPC
# Hiseq and Genome analyzer files were run separately and then merged - no other parameters changed between the two
### run CD-HIT duplicate removal (to mimic the Picard HMP tool) before proceeding with HMP trimming process

module load parallel

## set up array files (hiseq)
DATADIR=~/scratch/raw_stool/
arrayoffiles=($(ls $DATADIR | grep _1 | cut -f 1 -d _ | sort | uniq))

## make output directory to duplicated files
DUP=~/scratch/dups_hiseq/
mkdir -p $DUP

## make output directory to temporary files for quality control
TMP=~/scratch/tmp_dup_hiseq/
mkdir -p $TMP

# location of trimBWAstyle executable
TRIM=/path/to/bin/

# directory to cd-hit-dups executable (part of cd-hit auxtools -> not on conda)
CDHIT=/path/to/bin/cd-hit-v4.8.1-2019-0228/cd-hit-auxtools/

### print name of file to stout ####
echo "########################## WORKING ON SAMPLE"
echo ${arrayoffiles[${SLURM_ARRAY_TASK_ID}]}
echo "##########################"

#### ----- run duplicate removal ----

## skip this step if the output files already exist
if [ -e "$DUP/${arrayoffiles[${SLURM_ARRAY_TASK_ID}]}_1.nodup.fastq" ]; then
    echo "SKIPPING: $DUP/${arrayoffiles[${SLURM_ARRAY_TASK_ID}]}_1.nodup.fastq already exists"
    else
        # run CDHIT
        $CDHIT/cd-hit-dup -i $DATADIR/${arrayoffiles[${SLURM_ARRAY_TASK_ID}]}_1.fastq -i2 $DATADIR/${arrayoffiles[${SLURM_ARRAY_TASK_ID}]}_2.fastq -o $DUP/${arrayoffiles[${SLURM_ARRAY_TASK_ID}]}_1.nodup.fastq -o2 $DUP/${arrayoffiles[${SLURM_ARRAY_TASK_ID}]}_2.nodup.fastq -e 0.03
fi
echo "######### done with duplicate removal #############"

## ----- run HMP quality control ----

source ~/.bashrc
eval "$(micromamba shell hook --shell bash)"
micromamba activate samtools

## convert fastq to BAM
if [ -e "$TMP/${arrayoffiles[${SLURM_ARRAY_TASK_ID}]}.bam" ]; then 
    echo "SKIPPING: $TMP/${arrayoffiles[${SLURM_ARRAY_TASK_ID}]}.bam already exists"
    else
        # run picard
        picard FastqToSam -F1 $DUP/${arrayoffiles[${SLURM_ARRAY_TASK_ID}]}_1.nodup.fastq -F2 $DUP/${arrayoffiles[${SLURM_ARRAY_TASK_ID}]}_2.nodup.fastq -O $TMP/${arrayoffiles[${SLURM_ARRAY_TASK_ID}]}.nodup.bam -SAMPLE_NAME ${arrayoffiles[${SLURM_ARRAY_TASK_ID}]}.nodup
fi
echo "######## done with FastqToSam ##########"

## sort BAM 
if [ -e "$TMP/${arrayoffiles[${SLURM_ARRAY_TASK_ID}]}.nodup.sort.bam" ]; then 
    echo "SKIPPING: -e $TMP/${arrayoffiles[${SLURM_ARRAY_TASK_ID}]}.nodup.sort.bam already exists"
    else
        # run picard
        picard SortSam -I $TMP/${arrayoffiles[${SLURM_ARRAY_TASK_ID}]}.nodup.bam -O $TMP/${arrayoffiles[${SLURM_ARRAY_TASK_ID}]}.nodup.sort.bam -SORT_ORDER queryname
fi
echo "########## done with SortSam ############"

## trim
cd $TMP
perl $TRIM/trimBWAstyle.usingBam.pl -f $TMP/${arrayoffiles[${SLURM_ARRAY_TASK_ID}]}.nodup.sort.bam -q 3 -o 65 

echo "##### done with trimBWA #####"