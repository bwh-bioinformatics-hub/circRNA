#!/bin/bash

## This script takes air-ended RNAseq data, the genome index directory path, and the prefix of output files as input
## to run STAR alignments to SARS-Cov-2 virus genome.
##
#### USAGE: bash run-star.sh <output_dir> <sample_name> <ref_dir> <read1> <read2> 
## Current project reference: "human_vzv_ref"

output_dir=$1
sample_name=$2
ref_dir=$3
read1=$4
read2=$5

# Load required modules
module load glibc/2.14
module load star/2.7.3
module load fastqc

[ ! -s "$output_dir/$sample_name" ] && mkdir -p "$output_dir/$sample_name"

# Generate genome if genome parameters not existed
[ ! -s "$ref_dir/genomeParameters.txt" ] && \
STAR \
    --runMode genomeGenerate \
    --genomeFastaFiles "${ref_dir}/genome.fa" \
    --genomeDir "${ref_dir}" \
    --runThreadN 4 \
    --sjdbGTFfile "${ref_dir}/genome.gtf"

# Run genome alignment
fastqc --outdir "$output_dir/$sample_name" --extract -t 2 "$read1" "$read2" 

[ ! -s "$output_dir/$sample_name/Chimeric.out.junction" ] && \
STAR \
    --genomeDir "$ref_dir" \
    --runMode alignReads \
    --twopassMode Basic \
    --outFileNamePrefix "$output_dir/$sample_name/" \
    --readFilesIn "$read1" "$read2" \
    --outReadsUnmapped Fastx \
    --quantMode GeneCounts \
    --runThreadN 4 \
    --chimSegmentMin 10 \
    --chimJunctionOverhangMin 10 \
    --chimScoreMin 1 \
    --alignIntronMax 1000000 \
    --outFilterMismatchNoverReadLmax 0.02 \
    --alignTranscriptsPerReadNmax 100000 \
    --outSAMtype BAM SortedByCoordinate \
    --chimOutType Junctions SeparateSAMold \
    --outFilterMultimapNmax 2 \

