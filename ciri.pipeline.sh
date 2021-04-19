#!/bin/bash

## This script takes the sam file output from BWA-MEM and the prefix of output files as input
## to identify circular RNAs with CIRI2.
## It then computes the alternative splicing and visualize identified circular RNA candidates with CIRI_vis.  
##
#### USAGE: bash ciri.pipeline.sh <input-fastq-1> <input-fastq-2> <output-prefix> <ref-fasta> <annotation-gtf>

output_dir=$1
sample_name=$2
ref_dir=$3
read1=$4
read2=$5

num_cpu=4

[ ! -s "$output_dir/$sample_name" ] && mkdir -p "$output_dir/$sample_name"

# Alignment with BWA MEM - what is min alignment score(=19)?
bwa mem -T 19 -t $num_cpu \
	"$ref_dir/genome.fa" "$read1"  "$read2" \
	1> "$output_dir/$sample_name/bwa.sam" 2> "$output_dir/$sample_name/bwa.log"

# Identify circular RNAs with CIRI2
CIRI2.pl \
	-I "$output_dir/$sample_name/bwa.sam" \
	-O "$output_dir/$sample_name/CIRI2.txt" \
   	-F "$ref_dir/genome.fa" \
   	-A "$ref_dir/genome.gtf" \
   	-T $num_cpu

# Computes the alternative splicing
CIRI_AS.pl \
	-S "$output_dir/$sample_name/bwa.sam" \
   	-C "$output_dir/$sample_name/CIRI2.txt" \
   	-F "$ref_dir/genome.fa" \
   	-A "$ref_dir/genome.gtf" \
   	-O "$output_dir/$sample_name/CIRI_AS" \
   	-D yes

# Visualize identified circular RNA candidates with CIRI_vis
java -jar src/CIRI_vis.jar \
	-i "$output_dir/${sample_name}/CIRI_AS_jav.list" \
   	-l "$output_dir/${sample_name}/CIRI_AS_library_length.list" \
   	-r "$ref_dir/genome.fa" \
   	-min 1 \
   	-o "$output_dir/$sample_name/CIRI_vis"


