#!/bin/bash

set -euo pipefail

############################################
# ATAC-seq Blacklist Generation
############################################

BLACKLIST_DIR="blacklist"

mkdir -p "${BLACKLIST_DIR}"
mkdir -p logs


############################################
# Input annotation
############################################

GFF="${BLACKLIST_DIR}/GCF_000001735.4_TAIR10.1_genomic.gff"


############################################
# Create annotation blacklist
# Repeat regions and transposable elements
############################################

echo "=== Creating annotation blacklist ==="

grep -i \
    "repeat_region\|transposable_element" \
    "${GFF}" \
    | awk -F'\t' \
    '{print $1"\t"$4-1"\t"$5}' \
    > "${BLACKLIST_DIR}/annotation_blacklist.bed"


############################################
# Create combined BAM
############################################

echo "=== Merging filtered BAM files ==="

samtools merge \
    "${BLACKLIST_DIR}/all.bam" \
    alignment/filter/control1.fil.sorted.bam \
    alignment/filter/control2.fil.sorted.bam \
    alignment/filter/heat1.fil.sorted.bam \
    alignment/filter/heat2.fil.sorted.bam


############################################
# Genome coverage
############################################

echo "=== Calculating genome coverage ==="

bedtools genomecov \
    -ibam "${BLACKLIST_DIR}/all.bam" \
    -bg \
    > "${BLACKLIST_DIR}/coverage.bedgraph"


############################################
# Remove zero coverage
############################################

awk '$4 > 0' \
    "${BLACKLIST_DIR}/coverage.bedgraph" \
    > "${BLACKLIST_DIR}/coverage.filtered.bedgraph"


############################################
# Calculate mean and SD
############################################

echo "=== Calculating coverage threshold ==="

awk '
{
    sum += $4
    sum2 += $4*$4
}
END {
    mean = sum/NR
    sd = sqrt(sum2/NR - mean*mean)
    print "Mean coverage =", mean
    print "SD =", sd
    print "Threshold =", mean + 3*sd
}' "${BLACKLIST_DIR}/coverage.filtered.bedgraph" \
> "${BLACKLIST_DIR}/coverage_statistics.txt"


############################################
# Extract threshold automatically
############################################

THRESHOLD=$(awk '/Threshold/ {print $NF}' \
    "${BLACKLIST_DIR}/coverage_statistics.txt")

echo "Coverage blacklist threshold: ${THRESHOLD}"


############################################
# High coverage blacklist
############################################

awk -v threshold="${THRESHOLD}" \
    '$4 > threshold' \
    "${BLACKLIST_DIR}/coverage.filtered.bedgraph" \
    > "${BLACKLIST_DIR}/coverage_blacklist.bed"


############################################
# Merge high coverage regions
############################################

bedtools merge \
    -i "${BLACKLIST_DIR}/coverage_blacklist.bed" \
    > "${BLACKLIST_DIR}/coverage_blacklist_merged.bed"


############################################
# Final blacklist
############################################

cat \
    "${BLACKLIST_DIR}/annotation_blacklist.bed" \
    "${BLACKLIST_DIR}/coverage_blacklist_merged.bed" \
    | sort -k1,1 -k2,2n \
    | bedtools merge \
    > "${BLACKLIST_DIR}/final_blacklist.bed"


############################################
# Completion message
############################################

echo "=========================================="
echo "Blacklist generation completed."
echo "Final blacklist:"
echo "${BLACKLIST_DIR}/final_blacklist.bed"
echo "=========================================="
