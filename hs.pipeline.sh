#!/bin/bash

## This script takes the sam file output from BWA-MEM and the prefix of output files as input
## and identify circular RNAs with in-house R-script.
##
#### USAGE: bash hs.pipeline.sh <input-sam-file> <output-table-prefix>

## extract clipped junction reads in virus - CIGAR line: xx[SH]xxM | xxMxx[SH]
awk '($3 == "MN908947.3" && $7 == "="){if ($6 ~ /^([0-9]*M[0-9]*[SH])$/ ){print "E\t" $0} if( $6 ~ /^([0-9]*[SH][0-9]*M)$/) {print "S\t" $0}}' "${1}" > "${2}".virus.chimeric.sam 

## Run R Script to process reads
Rscript sortReads.r "${2}".virus.chimeric.sam "${2}".output.txt
rm "${2}".virus.chimeric.sam 

