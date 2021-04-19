#!/bin/bash

## This script takes pair-ended RNAseq data and the prefix of output files as input
## to identify circular RNA candidate for SARS-Cov-2 virus using find_circ.
##
#### USAGE: bash find_circ.pipeline.sh <output_dir> <sample_name> <read1> <read2>

# read1="data/CR2020.human/Calu3_totalRNA-S2-12h-A_R1.fastq"
# read2="data/CR2020.human/Calu3_totalRNA-S2-12h-A_R2.fastq"
# output_directory="output/CR2020.human/"
# sample_name="Calu3_totalRNA-S2-12h-A"

output_dir=$1
sample_name=$2
ref_dir=$3
read1=$4
read2=$5

num_cpu=4
# Load required modules
module load python/2.7.3 
module load pysam
module load bowtie2
module load samtools/1.4.1

# Create output directory if not existed
[ -s "$output_dir/$sample_name" ] || mkdir -p "$output_dir/$sample_name"

# Build genome index
bowtie2-build "$ref_dir/genome.fa" "$ref_dir"

# Mapping; Generate sorted alignments
[ "$read2" == "" ] && \
bowtie2 \
    -p $num_cpu \
    --very-sensitive \
    --score-min=C,-15,0 --mm \ # f(x) = -15 as the min score function; x = read length
    -x "$ref_dir/genome" \
    -q \
    -1 $read1 2> "$output_dir/$sample_name/bowtie2.log" | \
    samtools view -@ $num_cpu -hbuS - | \
    samtools sort -@ $num_cpu - -o "$output_dir/$sample_name/bowtie2.bam"

[ "$read2" != "" ] && \
bowtie2 \
    -p $num_cpu \
    --very-sensitive \
    --score-min=C,-15,0 --mm \
    -x "$ref_dir/genome" \
    -q \
    -1 $read1 -2 "$read2" 2> "$output_dir/$sample_name/bowtie2.log" | \
    samtools view -@ $num_cpu -hbuS - | \
    samtools sort -@ $num_cpu - -o "$output_dir/$sample_name/bowtie2.bam"

# Get the unmapped reads 
samtools view -h -f 4 -@ $num_cpu "$output_dir/$sample_name/bowtie2.bam" | \ # -f4 filter unmapped reads
	samtools view -Sb -@ $num_cpu - > "$output_dir/$sample_name/bowtie2_unmapped.bam"

# Pipe through unmapped2anchors.py to extract anchors
unmapped2anchors.py "$output_dir/$sample_name/bowtie2_unmapped.bam" | \
	gzip > "$output_dir/$sample_name/anchors.fastq.gz"

# Second comparison - screen both 
bowtie2  -p $num_cpu --score-min=C,-15,0 --reorder --mm \
        -q -U "$output_dir/$sample_name/anchors.fastq.gz" \
        -x "$ref_dir/genome" | \
    find_circ.py \
    --genome="$ref_dir/genome.fa" \
    --prefix="${sample_name}_" \
    --stats="$output_dir/$sample_name/stats.txt" \
    --reads="$output_dir/$sample_name/spliced_reads.fa" \
     > "$output_dir/$sample_name/splice_sites.bed"
                    
grep CIRCULAR "$output_dir/$sample_name/splice_sites.bed" | \
    # awk '$5>=2' | \ # more then 2 reads -> filter later when merging results together 
    grep UNAMBIGUOUS_BP | \
    grep ANCHOR_UNIQUE \
    > "$output_dir/$sample_name/circRNA_output.txt"
