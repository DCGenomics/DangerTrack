#!/bin/bash


# bash commands for preparing, cleaning, and binning the data


####################


# download and format references


# genome FASTA
rsync -a -P rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz .
tar -xvzf chromFa.tar.gz
cat chr*.fa > hg19.fa
rm chr*.fa

# genome 2bit
rsync -a -P rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit .

# chrom sizes
rsync -a -P rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes .

# repeatmasker
rsync -a -P rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/rmsk.txt.gz .
cat rmsk.txt.gz | gunzip | cut -f 6,7,8,12 > rmsk.bed

# Mappability or Uniqueness of Reference Genome from ENCODE
rsync -a -P rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign50mer.bigWig .

# add chr to GIAB_H002 BED
cat GIAB_H002.bed | awk -F $'\t' 'BEGIN {OFS=FS} {print "chr"$1,$2,$2,$3}' > GIAB_H002.chr.bed

# add chr to dbvar BED
cat dbvar_estd219.uniq.bed | awk -F $'\t' 'BEGIN {OFS=FS} {print "chr"$1,$2,$2,$3}' > dbvar_estd219.uniq.chr.bed


####################


# size for grouping events into windows/bins


bin_size="5000"
bin_bed="hg19.bin.${bin_size}.bed"

# divide genome to $bin_size windows and add interval-style regions name
bedtools makewindows -g hg19.chrom.sizes -w "$bin_size" | grep -v "_" | LC_ALL=C sort -k1,1 -k2,2n \
| awk -F $'\t' 'BEGIN {OFS=FS} {print $0,$1":"$2"-"$3}' > "$bin_bed"


####################


# summarize events by bin
# summary files for hg19 5kb bins are about 20 MB each


# bedtools coverage output includes original BED file columns (4 columns here)
# 1) The number of features in B that overlapped the A interval.
# 2) The number of bases in A that had non-zero coverage.
# 3) The length of the entry in A.
# 4) The fraction of bases in A that had non-zero coverage.

# fraction of bases in repeatmasker per bin
echo -e "#BIN\trmsk_fraction" > summary.rmsk.${bin_size}.txt
bedtools coverage -a "$bin_bed" -b rmsk.bed | cut -f 4,8 >> summary.rmsk.${bin_size}.txt

# GIAB_H002 break points
echo -e "#BIN\tGIAB_H002_events" > summary.GIAB_H002.${bin_size}.txt
bedtools coverage -a "$bin_bed" -b GIAB_H002.chr.bed | cut -f 4,5 >> summary.GIAB_H002.${bin_size}.txt

# dbvar_estd219 break points
echo -e "#BIN\tdbvar_estd219_events" > summary.dbvar_estd219.${bin_size}.txt
bedtools coverage -a "$bin_bed" -b dbvar_estd219.uniq.chr.bed | cut -f 4,5 >> summary.dbvar_estd219.${bin_size}.txt

# GC content (assumes 4-column BED file)
echo -e "#BIN\tGC_fraction" > summary.gc.${bin_size}.txt
bedtools nuc -fi hg19.fa -bed "$bin_bed" | cut -f 1,2,3,6 | grep -v "#" \
| awk -F $'\t' 'BEGIN {OFS=FS} {print $1":"$2"-"$3,$4}' >> summary.gc.${bin_size}.txt

# N content (assumes 4-column BED file)
echo -e "#BIN\tN_fraction" > summary.n.${bin_size}.txt
bedtools nuc -fi hg19.fa -bed "$bin_bed" | cut -f 1,2,3,11,13 | grep -v "#" \
| awk -F $'\t' 'BEGIN {OFS=FS} {print $1":"$2"-"$3,$4/$5}' >> summary.n.${bin_size}.txt

# average mapability per bin (using bigWigAverageOverBed from UCSC)
bigWigAverageOverBed wgEncodeCrgMapabilityAlign50mer.bigWig $bin_bed wgEncodeCrgMapabilityAlign50mer.${bin_size}.txt

# average mapability per bin in comparable format
echo -e "#BIN\tMapability50mer" > summary.map50mer.${bin_size}.txt
cat wgEncodeCrgMapabilityAlign50mer.${bin_size}.txt | cut -f 1,5 \
| tr ':' '\t' | tr '-' '\t' | LC_ALL=C sort -k1,1 -k2,2n \
| awk -F $'\t' 'BEGIN {OFS=FS} {print $1":"$2"-"$3,$4}' >> summary.map50mer.${bin_size}.txt


####################



# end
