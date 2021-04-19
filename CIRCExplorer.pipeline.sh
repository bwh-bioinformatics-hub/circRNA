#!/bin/bash

## This script takes pair-ended RNAseq data and the prefix of output files as input
## to run STAR alignments to SARS-Cov-2 virus genome and identify circular RNAs with CIRCExplorer2.
## Prequiste: do bash run-star.sh first
## NOTE: Need to use absolute path for reference directory
##
#### USAGE: bash CIRCExplorer.pipeline.sh <output_dir> <sample_name> <ref_dir> 
## Current project reference: "human_vzv_ref"

output_dir=$1
sample_name=$2
ref_dir=$3

cd "$output_dir/$sample_name/" || exit

# Convert Chimeric.out.junction to back_spliced_junction.bed
[ ! -f back_spliced_junction.bed ] && \
CIRCexplorer2 parse -t STAR Chimeric.out.junction 

# Identify circular RNAs with CIRCExplorer2
[ -s "$ref_dir" ] && [ ! -f circularRNA_known.txt ] && \
CIRCexplorer2 annotate \
    -r "$ref_dir/genome.refFlat" \
    -g "$ref_dir/genome.fa" \
    -b back_spliced_junction.bed \
    -o circularRNA_known.txt \
    --low-confidence \
    --no-fix 

