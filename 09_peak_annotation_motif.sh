#!/bin/bash

set -euo pipefail

############################################
# ATAC-seq Peak Annotation and Motif Analysis
# Arabidopsis thaliana TAIR10
############################################

############################################
# HOMER installation
############################################

echo "=== Installing HOMER ==="

wget http://homer.ucsd.edu/homer/configureHomer.pl

perl configureHomer.pl -install

perl configureHomer.pl -list

perl configureHomer.pl -install tair10


############################################
# Create BED file from CSAW significant peaks
############################################

echo "=== Creating BED file ==="

awk 'BEGIN{OFS="\t"}
NR>1{
    gsub("NC_003070.9","1",$1)
    gsub("NC_003071.7","2",$1)
    gsub("NC_003074.8","3",$1)
    gsub("NC_003075.7","4",$1)
    gsub("NC_003076.8","5",$1)
    print $1, ($2-1), $3, $4
}' \
treat_vs_control_csaw_DA-windows_significant.txt \
| sort -k1,1 -k2,2n \
> PEAKS.bed


############################################
# Peak annotation using HOMER
############################################

echo "=== Annotating peaks using HOMER ==="

annotatePeaks.pl \
    ALL_PEAKS.bed \
    tair10 \
    > ANNOTATED_ALL_PEAKS.txt


############################################
# Create BED file for upregulated peaks
############################################

echo "=== Creating BED file for upregulated peaks ==="

awk 'BEGIN{OFS="\t"}
NR>1 && $14=="up"{
    gsub("NC_003070.9","1",$1)
    gsub("NC_003071.7","2",$1)
    gsub("NC_003074.8","3",$1)
    gsub("NC_003075.7","4",$1)
    gsub("NC_003076.8","5",$1)
    print $1, ($2-1), $3
}' \
treat_vs_control_csaw_DA-windows_significant.txt \
> UP_PEAKS.bed


############################################
# Motif analysis of upregulated peaks
############################################

echo "=== Performing motif analysis for upregulated peaks ==="

findMotifsGenome.pl \
    UP_PEAKS.bed \
    tair10 \
    UP_MOTIFS/ \
    -size given


############################################
# Create BED file for downregulated peaks
############################################

echo "=== Creating BED file for downregulated peaks ==="

awk 'BEGIN{OFS="\t"}
NR>1 && $14=="down"{
    gsub("NC_003070.9","1",$1)
    gsub("NC_003071.7","2",$1)
    gsub("NC_003074.8","3",$1)
    gsub("NC_003075.7","4",$1)
    gsub("NC_003076.8","5",$1)
    print $1, ($2-1), $3
}' \
treat_vs_control_csaw_DA-windows_significant.txt \
> DOWN_PEAKS.bed


############################################
# Motif analysis of downregulated peaks
############################################

echo "=== Performing motif analysis for downregulated peaks ==="

findMotifsGenome.pl \
    DOWN_PEAKS.bed \
    tair10 \
    DOWN_MOTIFS/ \
    -size given


echo "=========================================="
echo "Peak annotation and motif analysis completed."
echo "=========================================="
