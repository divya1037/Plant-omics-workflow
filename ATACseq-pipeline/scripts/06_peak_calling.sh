#!/bin/bash

set -euo pipefail

############################################
# ATAC-seq Peak Calling
############################################

ATAC_DIR="alignment/filter/atac_format"
PEAK_DIR="${ATAC_DIR}/peakcalling"
BLACKLIST_DIR="blacklist"

mkdir -p "${PEAK_DIR}"
mkdir -p "${BLACKLIST_DIR}"
mkdir -p logs


############################################
# Individual peak calling
############################################

echo "=== Calling peaks for individual replicates ==="

for sample in control1 control2 heat1 heat2
do

    macs2 callpeak \
        -t "${ATAC_DIR}/${sample}.minimal.bedpe" \
        -f BEDPE \
        -n "${sample}" \
        --outdir "${PEAK_DIR}" \
        -g 1.191e8 \
        --broad \
        --broad-cutoff 0.05 \
        --keep-dup all

done


############################################
# Pooled control peaks
############################################

echo "=== Calling pooled control peaks ==="

macs2 callpeak \
    -t \
    "${ATAC_DIR}/control1.minimal.bedpe" \
    "${ATAC_DIR}/control2.minimal.bedpe" \
    -f BEDPE \
    -n control_pool \
    --outdir "${PEAK_DIR}" \
    -g 1.191e8 \
    --broad \
    --broad-cutoff 0.05 \
    --keep-dup all


############################################
# Pooled heat peaks
############################################

echo "=== Calling pooled heat peaks ==="

macs2 callpeak \
    -t \
    "${ATAC_DIR}/heat1.minimal.bedpe" \
    "${ATAC_DIR}/heat2.minimal.bedpe" \
    -f BEDPE \
    -n heat_pool \
    --outdir "${PEAK_DIR}" \
    -g 1.191e8 \
    --broad \
    --broad-cutoff 0.05 \
    --keep-dup all




############################################
# Filter individual peaks
############################################

echo "=== Filtering peaks using blacklist ==="

for sample in control1 control2 heat1 heat2
do

    bedtools intersect \
        -v \
        -a "${PEAK_DIR}/${sample}_peaks.broadPeak" \
        -b "${BLACKLIST_DIR}/final_blacklist.bed" \
        > "${PEAK_DIR}/${sample}_peaks.filt.broadPeak"

done


############################################
# Filter pooled peaks
############################################

bedtools intersect \
    -v \
    -a "${PEAK_DIR}/heat_pool_peaks.broadPeak" \
    -b "${BLACKLIST_DIR}/final_blacklist.bed" \
    > "${PEAK_DIR}/heat_pool_peaks.filt.broadPeak"


bedtools intersect \
    -v \
    -a "${PEAK_DIR}/control_pool_peaks.broadPeak" \
    -b "${BLACKLIST_DIR}/final_blacklist.bed" \
    > "${PEAK_DIR}/control_pool_peaks.filt.broadPeak"


############################################
# Reproducible peaks
############################################

echo "=== Finding reproducible heat peaks ==="

bash tools/naiveOverlapBroad.sh \
    "${PEAK_DIR}/heat1_peaks.filt.broadPeak" \
    "${PEAK_DIR}/heat2_peaks.filt.broadPeak" \
    "${PEAK_DIR}/heat_pool_peaks.filt.broadPeak" \
    > "${PEAK_DIR}/heat_overlap_peaks.filt.broadPeak"


echo "=== Finding reproducible control peaks ==="

bash tools/naiveOverlapBroad.sh \
    "${PEAK_DIR}/control1_peaks.filt.broadPeak" \
    "${PEAK_DIR}/control2_peaks.filt.broadPeak" \
    "${PEAK_DIR}/control_pool_peaks.filt.broadPeak" \
    > "${PEAK_DIR}/control_overlap_peaks.filt.broadPeak"


echo "=========================================="
echo "Peak calling completed successfully."
echo "=========================================="
