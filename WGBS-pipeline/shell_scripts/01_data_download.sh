#!/bin/bash

###############################################################
# WGBS Computational Protocol
# Script: 01_data_download.sh
#
# Description:
#   Downloads whole-genome bisulfite sequencing (WGBS) data
#   from the NCBI Sequence Read Archive (SRA), converts the
#   downloaded files into compressed FASTQ format, and prepares
#   them for downstream analysis.
#
# Input:
#   SRA accession numbers
#
# Output:
#   raw_data/
#      ├── SRR5494752.fastq.gz
#      └── SRR5494755.fastq.gz
#
# Software:
#   SRA Toolkit (v3.0 or later)
#
# Tested Version:
#   SRA Toolkit v3.0.7
#
# Author:
#   Surapuram Aswini et al.
###############################################################

set -e
set -o pipefail

echo "=============================================="
echo "Step 1 : Downloading WGBS Data"
echo "Started : $(date)"
echo "=============================================="

###############################
# Create directories
###############################

mkdir -p raw_data

###############################
# Move to working directory
###############################

cd raw_data

###############################
# Download SRA datasets
###############################

echo "Downloading SRR5494752 ..."
prefetch SRR5494752

echo "Downloading SRR5494755 ..."
prefetch SRR5494755

###############################
# Convert SRA to FASTQ
###############################

echo "Converting SRR5494752 ..."
fasterq-dump SRR5494752

echo "Converting SRR5494755 ..."
fasterq-dump SRR5494755

###############################
# Compress FASTQ files
###############################

gzip SRR5494752.fastq
gzip SRR5494755.fastq

###############################
# Verify output
###############################

echo ""
echo "Downloaded FASTQ files:"

ls -lh *.fastq.gz

###############################
# Finish
###############################

echo ""
echo "=============================================="
echo "Data download completed successfully."
echo "Finished : $(date)"
echo "=============================================="