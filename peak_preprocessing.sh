#!/bin/bash

# Craig Bielski
# combine narrowPeak files, sort, merge overlapping peaks, and convert to SAF

# combine and sort peaks with hacky implementation of sort -V for OSX
awk 'FNR>1 || NR==1' ~/peak_data/*/*/*narrowPeak | sed -e 's/chrX/chr23/' -e 's/chrY/chr24/' | sort -k1.4 -n -k2,2n | sed -e 's/chr23/chrX/' -e 's/chr24/chrY/' > all_peaks.sort

# bedtools merge, default 10bp threshold
bedtools merge -i all_peaks.sort > all_peaks.merge.sort.bed

# convert to SAF format for featureCounts
awk -v OFS="\t" 'BEGIN {print "GeneID","Chr","Start","End","Strand"} { print "Peak_"NR,$1,$2,$3,"."}' all_peaks.merge.sort.bed > all_peaks.merge.sort.saf

