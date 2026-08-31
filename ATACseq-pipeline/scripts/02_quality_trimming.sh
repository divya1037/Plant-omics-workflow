#!/bin/bash

set -euo pipefail

############################################
# ATAC-seq Quality Control and Trimming
############################################

mkdir -p raw_fastqc_reports
mkdir -p filtered_reads
mkdir -p filtered_reads/filtered_fastqc_reports
mkdir -p logs


############################################
# Raw read FastQC
############################################

echo "=== Running FastQC on raw reads ==="

fastqc raw_reads/*.fastq.gz \
    -o raw_fastqc_reports/ \
    --threads 8


############################################
# Raw MultiQC
############################################

echo "=== Running MultiQC on raw reads ==="

mkdir -p raw_fastqc_reports/raw_multiqc_reports

multiqc raw_fastqc_reports/ \
    -o raw_fastqc_reports/raw_multiqc_reports/


############################################
# Adapter trimming and quality filtering
############################################

echo "=== Trimming control1 ==="

trim_galore \
    -q 20 \
    --fastqc \
    --paired \
    raw_reads/control1_R1.fastq.gz \
    raw_reads/control1_R2.fastq.gz \
    -o filtered_reads/


echo "=== Trimming control2 ==="

trim_galore \
    -q 20 \
    --fastqc \
    --paired \
    raw_reads/control2_R1.fastq.gz \
    raw_reads/control2_R2.fastq.gz \
    -o filtered_reads/


echo "=== Trimming heat1 ==="

trim_galore \
    -q 20 \
    --fastqc \
    --paired \
    raw_reads/heat1_R1.fastq.gz \
    raw_reads/heat1_R2.fastq.gz \
    -o filtered_reads/


echo "=== Trimming heat2 ==="

trim_galore \
    -q 20 \
    --fastqc \
    --paired \
    raw_reads/heat2_R1.fastq.gz \
    raw_reads/heat2_R2.fastq.gz \
    -o filtered_reads/


############################################
# MultiQC after trimming
############################################

echo "=== Running MultiQC on trimmed reads ==="

mkdir -p filtered_reads/filtered_multiqc_reports

multiqc filtered_reads/ \
    -o filtered_reads/filtered_multiqc_reports/


echo "=========================================="
echo "Quality control and trimming completed."
echo "=========================================="
