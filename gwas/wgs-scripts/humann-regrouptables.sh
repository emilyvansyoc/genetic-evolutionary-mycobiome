#!/bin/bash

# submitted as a batch job

### get humann regrouped tables
source ~/.bashrc
conda activate humenv
set -uex

## set input directory
DIR=$1

#### regroup metacyc
humann_regroup_table -i $DIR/all_genefamilies.tsv --groups uniref90_rxn -o $DIR/all_metacyc.tsv

echo "##### done with Metacyc ########"

### regroup uniref EC names
humann_regroup_table -i $DIR/all_genefamilies.tsv --groups uniref90_level4ec -o $DIR/all_l4ec.tsv

echo "###### done with UnirefL4EC ########"

## add names to feature IDs (can increase file size)
humann_rename_table -i $DIR/all_genefamilies.tsv -n uniref90 -o $DIR/all_genefamilies_renamed.tsv

echo "##### done with Renaming #####"