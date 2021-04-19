#!/bin/bash

## This script takes directory to files with sample names as prefix output from STAR and the prefix of output files as input
## to identify circular RNAs with circRNA_finder.
##
#### USAGE: bash circRNA_finder.pipeline.sh <star-output-directory> <output-prefix>

/data/bioinformatics/external_data/circrna/src/circRNA_finder/postProcessStarAlignment.pl \
    --starDir "$1" \
    --outDir "$2" \
    --minLen 19

