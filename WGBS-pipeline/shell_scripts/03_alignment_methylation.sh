#!/bin/bash

###############################################################
# WGBS Computational Protocol
# Script: 03_alignment_methylation.sh
#
# Description:
#   Performs reference genome preparation, bisulfite read
#   alignment, BAM deduplication, methylation extraction,
#   and generation of methylation reports.
#
# Input:
#   reference_genome/
#       Arabidopsis.fa
#
#   trimmed/
#       *_trimmed.fq.gz
#
# Output:
#   alignment/
#   results/
#
# Software:
#   Bismark v0.25.1
#   Bowtie2 v2.5.1
#   SAMtools v1.22
#
# Usage:
#   bash shell_scripts/03_alignment_methylation.sh
#
# Author:
#   Surapuram Aswini et al.
###############################################################

set -e
set -o pipefail

echo "=================================================="
echo "Step 3 : Alignment and Methylation Calling"
echo "Started : $(date)"
echo "=================================================="

###############################
# Check input directories
###############################

if [ ! -d "reference_genome" ]; then
    echo "ERROR: reference_genome directory not found."
    exit 1
fi

if [ ! -d "trimmed" ]; then
    echo "ERROR: trimmed directory not found."
    exit 1
fi

###############################
# Create output directories
###############################

mkdir -p alignment
mkdir -p results

###############################
# Genome Preparation
###############################

echo ""
echo "Preparing Bismark genome..."

bismark_genome_preparation \
    --bowtie2 \
    reference_genome/

###############################
# Read Alignment
###############################

echo ""
echo "Running Bismark alignment..."

for file in trimmed/*_trimmed.fq.gz
do

    echo "Processing $(basename "$file")"

    bismark \
        --genome reference_genome/ \
        --output_dir alignment \
        --bowtie2 \
        "$file"

done

###############################
# Deduplication
###############################

echo ""
echo "Removing PCR duplicates..."

for bam in alignment/*_bismark_bt2.bam
do

    deduplicate_bismark \
        "$bam"

done

###############################
# Methylation Extraction
###############################

echo ""
echo "Extracting methylation information..."

for bam in alignment/*deduplicated.bam
do

    bismark_methylation_extractor \
        --single-end \
        --comprehensive \
        --bedGraph \
        --cytosine_report \
        --gzip \
        --genome_folder reference_genome \
        --output results \
        "$bam"

done

###############################
# Verify Outputs
###############################

echo ""
echo "Generated files:"

ls results

###############################
# Finish
###############################

echo ""
echo "=================================================="
echo "Alignment and methylation calling completed."
echo "Finished : $(date)"
echo "=================================================="