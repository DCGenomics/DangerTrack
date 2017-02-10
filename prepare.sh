#!/bin/bash


# bash commands to prepare, clean, and bin the data


####################


# download and format references


# genome FASTA
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz
tar -xvzf chromFa.tar.gz
cat chr*.fa > hg19.fa
rm chr*.fa

# chrom sizes
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes

# mappability or uniqueness of reference genome from ENCODE
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign50mer.bigWig
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign100mer.bigWig

# add chr to GIAB_H002 BED
cat GIAB_H002.bed | awk -F $'\t' 'BEGIN {OFS=FS} {print "chr"$1,$2,$2,$3}' > GIAB_H002.chr.bed

# add chr to dbvar BED
cat dbvar_estd219.uniq.bed | awk -F $'\t' 'BEGIN {OFS=FS} {print "chr"$1,$2,$2,$3}' > dbvar_estd219.uniq.chr.bed

# download GRC issues (GRCh37.p13_issues.gff3)
wget ftp://ftp.ncbi.nlm.nih.gov/pub/grc/human/GRC/Issue_Mapping/GRCh37.p13_issues.gff3

# convert GRC issues to BED
cat GRCh37.p13_issues.gff3 | grep -v "^#" | grep -v "chr=Un" \
| awk -F $'\t' 'BEGIN {OFS=FS} {print $9,$4-1,$5}' \
| perl -pe 's/Name.*chr=(.*?);.*?\t/chr\1\t/' | LC_ALL=C sort -k1,1 -k2,2n \
> GRC_issues.bed

# create a merged GRC issues BED
bedtools merge -i GRC_issues.bed > GRC_issues.merged.bed

# download ENCODE DAC Blacklisted Regions
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeDacMapabilityConsensusExcludable.bed.gz

# convert ENCODE DAC Blacklisted Regions to BED
gunzip -c wgEncodeDacMapabilityConsensusExcludable.bed.gz \
| cut -f 1-4 | LC_ALL=C sort -k1,1 -k2,2n \
> ENCODE_DAC_blacklisted.bed

# create a merged ENCODE DAC Blacklisted Regions BED
bedtools merge -i ENCODE_DAC_blacklisted.bed > ENCODE_DAC_blacklisted.merged.bed

# combined GRC and ENCODE issues
cat GRC_issues.bed ENCODE_DAC_blacklisted.bed | cut -f 1-3 | LC_ALL=C sort -k1,1 -k2,2n | bedtools merge > issues.bed


####################


# size for grouping events into windows/bins


bin_size="5000"
bin_bed="hg19.bin.${bin_size}.bed"

# divide genome to $bin_size windows and add interval-style regions name
bedtools makewindows -g hg19.chrom.sizes -w "$bin_size" | grep -v "_" \
| LC_ALL=C sort -k1,1 -k2,2n | awk -F $'\t' 'BEGIN {OFS=FS} {print $0,$1":"$2"-"$3}' \
> "$bin_bed"


####################


# summarize events by bin
# summary files for hg19 5kb bins are about 20 MB each


# bedtools coverage outputs all original BED file columns (4 columns here) and adds:
# 1) The number of features in B that overlapped the A interval.
# 2) The number of bases in A that had non-zero coverage.
# 3) The length of the entry in A.
# 4) The fraction of bases in A that had non-zero coverage.

# GRC known issues (fraction of window with any issues)
echo -e "#BIN\tGRC_issues" > summary.GRC_issues.${bin_size}.txt
bedtools coverage -a "$bin_bed" -b GRC_issues.bed | cut -f 4,8 >> summary.GRC_issues.${bin_size}.txt

# ENCODE DAC Blacklisted Regions (fraction of window with any issues)
echo -e "#BIN\tENCODE_DAC_blacklisted" > summary.ENCODE_DAC_blacklisted.${bin_size}.txt
bedtools coverage -a "$bin_bed" -b ENCODE_DAC_blacklisted.bed | cut -f 4,8 >> summary.ENCODE_DAC_blacklisted.${bin_size}.txt

# GIAB break points
echo -e "#BIN\tevents_GIAB" > summary.GIAB_H002.${bin_size}.txt
bedtools coverage -a "$bin_bed" -b GIAB_H002.chr.bed | cut -f 4,5 >> summary.GIAB_H002.${bin_size}.txt

# 1KG break points
echo -e "#BIN\tevents_1KG" > summary.dbvar_estd219.${bin_size}.txt
bedtools coverage -a "$bin_bed" -b dbvar_estd219.uniq.chr.bed | cut -f 4,5 >> summary.dbvar_estd219.${bin_size}.txt

# average mappability per bin (using bigWigAverageOverBed from UCSC)
bigWigAverageOverBed wgEncodeCrgMapabilityAlign50mer.bigWig $bin_bed wgEncodeCrgMapabilityAlign50mer.${bin_size}.txt
bigWigAverageOverBed wgEncodeCrgMapabilityAlign100mer.bigWig $bin_bed wgEncodeCrgMapabilityAlign100mer.${bin_size}.txt

# average 50bp mappability per bin in comparable format
echo -e "#BIN\tMapability50mer" > summary.map50mer.${bin_size}.txt
cat wgEncodeCrgMapabilityAlign50mer.${bin_size}.txt | cut -f 1,5 \
| tr ':' '\t' | tr '-' '\t' | LC_ALL=C sort -k1,1 -k2,2n \
| awk -F $'\t' 'BEGIN {OFS=FS} {print $1":"$2"-"$3,$4}' >> summary.map50mer.${bin_size}.txt

# average 100bp mappability per bin in comparable format
echo -e "#BIN\tMapability100mer" > summary.map100mer.${bin_size}.txt
cat wgEncodeCrgMapabilityAlign100mer.${bin_size}.txt | cut -f 1,5 \
| tr ':' '\t' | tr '-' '\t' | LC_ALL=C sort -k1,1 -k2,2n \
| awk -F $'\t' 'BEGIN {OFS=FS} {print $1":"$2"-"$3,$4}' >> summary.map100mer.${bin_size}.txt


####################



# end
