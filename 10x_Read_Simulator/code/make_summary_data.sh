#!/bin/bash

## USAGE: make_summary_data.sh depth_coverage_file.gz chromosome_list.txt <number of rows>
## make_summary_data.sh data/depth13Genome.depth.gz hg38_chrom_list.txt 1000 > output_file.txt

## DESCRIPTION: This script will ouput a subset of rows from a depth file 
## for use as a summary data set

input_file="$1" # depth coverage file .gz (needs .tbi present)
chrom_list="$2" # list of chrom's to use
nrows="$3" # number of rows to save

cat $chrom_list | while read i; do
if [[ ! -z "$i" ]]; then
chrom="${i}"
tabix $input_file $chrom | head -$nrows


fi
done
