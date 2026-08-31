#!/bin/bash

set -euo pipefail

############################################
# ATAC-seq Alignment
# Reference: Arabidopsis TAIR10.1
############################################

REFERENCE="reference/GCF_000001735.4_TAIR10.1_genomic.fna"
INDEX="reference/bowtie2_index/index"

mkdir -p reference/bowtie2_index
mkdir -p alignment
mkdir -p logs


############################################
# Build Bowtie2 index
############################################

echo "=== Building Bowtie2 index ==="

if [ ! -f "${INDEX}.1.bt2" ] && [ ! -f "${INDEX}.1.bt2l" ]; then

    bowtie2-build \
        "${REFERENCE}" \
        "${INDEX}"

else

    echo "Bowtie2 index already exists. Skipping."

fi


############################################
# Alignment
############################################

echo "=== Aligning control1 ==="

bowtie2 \
    --very-sensitive \
    -X 1000 \
    -x "${INDEX}" \
    -1 filtered_reads/control1_R1_val_1.fq.gz \
    -2 filtered_reads/control1_R2_val_2.fq.gz \
    2> logs/control1_bowtie2.log \
    | samtools view -bS - \
    -o alignment/control1.bam


echo "=== Aligning control2 ==="

bowtie2 \
    --very-sensitive \
    -X 1000 \
    -x "${INDEX}" \
    -1 filtered_reads/control2_R1_val_1.fq.gz \
    -2 filtered_reads/control2_R2_val_2.fq.gz \
    2> logs/control2_bowtie2.log \
    | samtools view -bS - \
    -o alignment/control2.bam


echo "=== Aligning heat1 ==="

bowtie2 \
    --very-sensitive \
    -X 1000 \
    -x "${INDEX}" \
    -1 filtered_reads/heat1_R1_val_1.fq.gz \
    -2 filtered_reads/heat1_R2_val_2.fq.gz \
    2> logs/heat1_bowtie2.log \
    | samtools view -bS - \
    -o alignment/heat1.bam


echo "=== Aligning heat2 ==="

bowtie2 \
    --very-sensitive \
    -X 1000 \
    -x "${INDEX}" \
    -1 filtered_reads/heat2_R1_val_1.fq.gz \
    -2 filtered_reads/heat2_R2_val_2.fq.gz \
    2> logs/heat2_bowtie2.log \
    | samtools view -bS - \
    -o alignment/heat2.bam


############################################
# Sort and index
############################################

echo "=== Sorting and indexing BAM files ==="

for sample in control1 control2 heat1 heat2
do

    samtools sort \
        alignment/${sample}.bam \
        -o alignment/${sample}.sorted.bam

    samtools index \
        alignment/${sample}.sorted.bam

done


############################################
# Mapping statistics
############################################

echo "=== Generating idxstats ==="

for sample in control1 control2 heat1 heat2
do

    samtools idxstats \
        alignment/${sample}.sorted.bam \
        > alignment/${sample}_idxstats.txt

done


############################################
# Remove mitochondrial and chloroplast reads
############################################

echo "=== Removing mitochondrial and chloroplast reads ==="

for sample in control1 control2 heat1 heat2
do

    samtools view -h \
        alignment/${sample}.sorted.bam \
        | grep -v 'NC_037304.1' \
        | grep -v 'NC_000932.1' \
        | samtools sort \
        -O bam \
        -o alignment/${sample}.noMtChlo.bam

    samtools index \
        alignment/${sample}.noMtChlo.bam

    samtools idxstats \
        alignment/${sample}.noMtChlo.bam \
        > alignment/${sample}_noMtChlo_idxstats.txt

done


echo "=========================================="
echo "Alignment completed successfully."
echo "=========================================="
