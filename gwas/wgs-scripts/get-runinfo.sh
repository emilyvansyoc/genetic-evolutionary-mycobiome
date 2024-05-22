## run locally

## script to get runinfo and test accession number

# activate conda
source ~/miniconda3/etc/profile.d/conda.sh
conda activate myenv

# catch errors
set -uex

## ---- CHANGE THESE VARIABLES ----

# accession number
ACC=PRJNA48479

# input directory
DIR=/storage/group/srb6251/default/emilyvansyoc/hmp-wgs/

# general output name
OUTNAME=hmp48479

# set prefix (European is ERR, American is SRR)
PREF=SRR

## ---- workflow: get runinfo and test accession ----

# get runinfo
esearch -db sra -query $ACC | efetch -format runinfo > $DIR/runinfo_$OUTNAME.csv

# get runids
cat $DIR/runinfo_$OUTNAME.csv | cut -f 1 -d ',' | grep $PREF > $DIR/runids_$OUTNAME.txt

# get one sample with 10,000 reads to test connection
cat $DIR/runids_$OUTNAME.txt | head -1 | parallel fastq-dump -X 10000 --split-files --outdir $DIR {}

# print stats for the reads
seqkit stat $DIR/*.fastq

# print finished message
echo "done!"

