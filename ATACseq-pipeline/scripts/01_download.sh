#!/bin/bash

set -euo pipefail

############################################
# ATAC-seq Download Script
# Arabidopsis thaliana TAIR10.1
############################################

mkdir -p raw_reads
mkdir -p reference
mkdir -p blacklist

############################################
# Download ATAC-seq FASTQ files
############################################

echo "=== Downloading control1 ==="

wget -P raw_reads/ \
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR826/DRR826752/DRR826752_1.fastq.gz

wget -P raw_reads/ \
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR826/DRR826752/DRR826752_2.fastq.gz

mv raw_reads/DRR826752_1.fastq.gz raw_reads/control1_R1.fastq.gz
mv raw_reads/DRR826752_2.fastq.gz raw_reads/control1_R2.fastq.gz


echo "=== Downloading control2 ==="

wget -P raw_reads/ \
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR826/DRR826753/DRR826753_1.fastq.gz

wget -P raw_reads/ \
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR826/DRR826753/DRR826753_2.fastq.gz

mv raw_reads/DRR826753_1.fastq.gz raw_reads/control2_R1.fastq.gz
mv raw_reads/DRR826753_2.fastq.gz raw_reads/control2_R2.fastq.gz


echo "=== Downloading heat1 ==="

wget -P raw_reads/ \
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR826/DRR826754/DRR826754_1.fastq.gz

wget -P raw_reads/ \
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR826/DRR826754/DRR826754_2.fastq.gz

mv raw_reads/DRR826754_1.fastq.gz raw_reads/heat1_R1.fastq.gz
mv raw_reads/DRR826754_2.fastq.gz raw_reads/heat1_R2.fastq.gz


echo "=== Downloading heat2 ==="

wget -P raw_reads/ \
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR826/DRR826755/DRR826755_1.fastq.gz

wget -P raw_reads/ \
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR826/DRR826755/DRR826755_2.fastq.gz

mv raw_reads/DRR826755_1.fastq.gz raw_reads/heat2_R1.fastq.gz
mv raw_reads/DRR826755_2.fastq.gz raw_reads/heat2_R2.fastq.gz


############################################
# Download Arabidopsis TAIR10.1 reference
############################################

echo "=== Downloading reference genome ==="

wget -P reference/ \
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.fna.gz

gunzip -f reference/GCF_000001735.4_TAIR10.1_genomic.fna.gz


############################################
# Download genome annotation
############################################

echo "=== Downloading genome annotation ==="

wget -P blacklist/ \
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.gff.gz

gunzip -f blacklist/GCF_000001735.4_TAIR10.1_genomic.gff.gz


echo "=========================================="
echo "Download completed successfully."
echo "=========================================="
