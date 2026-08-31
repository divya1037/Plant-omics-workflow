#!/bin/bash

set -euo pipefail

############################################
# ATAC-seq Post-alignment Processing
############################################

mkdir -p alignment/filter
mkdir -p alignment/filter/atac_format
mkdir -p logs


############################################
# Add Read Groups
############################################

echo "=== Adding read groups ==="

for sample in control1 control2 heat1 heat2
do

    java -jar picard.jar AddOrReplaceReadGroups \
        I=alignment/${sample}.noMtChlo.bam \
        O=alignment/${sample}.noMtChlo.RG.bam \
        RGID=${sample} \
        RGLB=lib1 \
        RGPL=ILLUMINA \
        RGPU=unit1 \
        RGSM=${sample}

done


############################################
# Remove duplicates
############################################

echo "=== Removing duplicates ==="

for sample in control1 control2 heat1 heat2
do

    java -jar picard.jar MarkDuplicates \
        I=alignment/${sample}.noMtChlo.RG.bam \
        O=alignment/${sample}.noMtChlo.nodup.bam \
        M=alignment/${sample}.dup_metrics.txt \
        REMOVE_DUPLICATES=true

done


############################################
# Filter:
# MAPQ >= 20
# properly paired reads
############################################

echo "=== Filtering BAM files ==="

for sample in control1 control2 heat1 heat2
do

    samtools view \
        -b \
        -f 2 \
        alignment/${sample}.noMtChlo.nodup.bam \
        -o alignment/filter/${sample}.fil.bam

    samtools sort \
        alignment/filter/${sample}.fil.bam \
        -o alignment/filter/${sample}.fil.sorted.bam

    samtools index \
        alignment/filter/${sample}.fil.sorted.bam

done


############################################
# ATAC-seq formatting
############################################

echo "=== Sorting BAM files by read name ==="

for sample in control1 control2 heat1 heat2
do

    samtools sort \
        -n \
        -o alignment/filter/atac_format/${sample}.namesorted.bam \
        alignment/filter/${sample}.fil.sorted.bam

done


############################################
# Fix read mates
############################################

echo "=== Fixing read mates ==="

for sample in control1 control2 heat1 heat2
do

    samtools fixmate \
        alignment/filter/atac_format/${sample}.namesorted.bam \
        alignment/filter/atac_format/${sample}.fixed.bam

done


############################################
# Convert BAM to BEDPE
############################################

echo "=== Converting BAM to BEDPE ==="

for sample in control1 control2 heat1 heat2
do

    samtools view \
        -bf 0x2 \
        alignment/filter/atac_format/${sample}.fixed.bam \
        | bedtools bamtobed \
            -i stdin \
            -bedpe \
        > alignment/filter/atac_format/${sample}.fixed.bedpe

done


############################################
# Tn5 shifting
############################################

echo "=== Tn5 shifting ==="

for sample in control1 control2 heat1 heat2
do

    bash tools/bedpeTn5shift.sh \
        alignment/filter/atac_format/${sample}.fixed.bedpe \
        > alignment/filter/atac_format/${sample}.tn5.bedpe

done


############################################
# Convert BEDPE to MACS2-compatible format
############################################

echo "=== Converting BEDPE format for MACS2 ==="

for sample in control1 control2 heat1 heat2
do

    bash tools/bedpeMinimalConvert.sh \
        alignment/filter/atac_format/${sample}.tn5.bedpe \
        > alignment/filter/atac_format/${sample}.minimal.bedpe

done


echo "=========================================="
echo "Post-alignment processing completed."
echo "=========================================="
