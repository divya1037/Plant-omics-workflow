#!/bin/bash

###############################################################
# WGBS Computational Protocol
# Script: 02_preprocessing.sh
#
# Description:
#   Performs quality assessment of raw sequencing reads,
#   adapter trimming, and post-trimming quality assessment.
#
# Input:
#   raw_data/*.fastq.gz
#
# Output:
#   qc/
#       ├── Raw_FastQC/
#       ├── Raw_MultiQC/
#       ├── Trimmed_FastQC/
#       └── Trimmed_MultiQC/
#
#   trimmed/
#       ├── *_trimmed.fq.gz
#       └── *_trimming_report.txt
#
# Software:
#   FastQC v0.12.1
#   MultiQC v1.20
#   Trim Galore! v0.6.10
#   Cutadapt v4.4
#
# Usage:
#   bash shell_scripts/02_preprocessing.sh
#
# Author:
#   Surapuram Aswini et al.
###############################################################

set -e
set -o pipefail

echo "=================================================="
echo "Step 2 : Quality Control and Preprocessing"
echo "Started : $(date)"
echo "=================================================="

###############################
# Check input directory
###############################

if [ ! -d "raw_data" ]; then
    echo "ERROR: raw_data directory not found."
    exit 1
fi

###############################
# Create output directories
###############################

mkdir -p qc/Raw_FastQC
mkdir -p qc/Raw_MultiQC
mkdir -p qc/Trimmed_FastQC
mkdir -p qc/Trimmed_MultiQC
mkdir -p trimmed

###############################
# Check input FASTQ files
###############################

FASTQ_FILES=$(find raw_data -name "*.fastq.gz")

if [ -z "$FASTQ_FILES" ]; then
    echo "ERROR: No FASTQ files found."
    exit 1
fi

echo ""
echo "Detected FASTQ files:"
echo "$FASTQ_FILES"
echo ""

###############################
# Step 2.1 Raw FastQC
###############################

echo "Running FastQC on raw reads..."

fastqc \
    raw_data/*.fastq.gz \
    --outdir qc/Raw_FastQC

###############################
# Step 2.2 Raw MultiQC
###############################

echo "Generating Raw MultiQC report..."

multiqc \
    qc/Raw_FastQC \
    --outdir qc/Raw_MultiQC

###############################
# Step 2.3 Adapter Trimming
###############################

echo "Running Trim Galore..."

for file in raw_data/*.fastq.gz
do

    trim_galore \
        --quality 20 \
        --length 20 \
        --fastqc \
        --output_dir trimmed \
        "$file"

done

###############################
# Step 2.4 FastQC on Trimmed Reads
###############################

echo "Running FastQC on trimmed reads..."

fastqc \
    trimmed/*.fq.gz \
    --outdir qc/Trimmed_FastQC

###############################
# Step 2.5 MultiQC
###############################

echo "Generating Trimmed MultiQC report..."

multiqc \
    qc/Trimmed_FastQC \
    --outdir qc/Trimmed_MultiQC

###############################
# Verify outputs
###############################

echo ""
echo "Trimmed files:"

ls -lh trimmed/*.fq.gz

###############################
# Finish
###############################

echo ""
echo "=================================================="
echo "Preprocessing completed successfully."
echo "Finished : $(date)"
echo "=================================================="