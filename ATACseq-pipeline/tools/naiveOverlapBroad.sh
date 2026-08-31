#!/bin/bash

############################################
# Reproducible broad-peak identification
#
# Adapted from:
# Jake Reske, Michigan State University, 2019
#
# Purpose:
# Identify pooled broad peaks supported by both
# individual biological replicates.
#
# A peak is considered overlapping when the overlap
# represents at least 50% of either member of the
# overlapping peak pair.
#
# Inputs:
#   1. Replicate 1 broadPeak
#   2. Replicate 2 broadPeak
#   3. Pooled broadPeak
#
# Usage:
# bash naiveOverlapBroad.sh \
#     replicate1.broadPeak \
#     replicate2.broadPeak \
#     pooled.broadPeak \
#     > reproducible.broadPeak
############################################

REP1="$1"
REP2="$2"
POOL="$3"

if [[ $# -ne 3 ]]; then
    echo "Usage: bash naiveOverlapBroad.sh REP1.broadPeak REP2.broadPeak POOL.broadPeak" >&2
    exit 1
fi


############################################
# First require overlap between pooled peaks
# and replicate 1.
############################################

bedtools intersect -wo \
    -a "${POOL}" \
    -b "${REP1}" |
awk '
BEGIN {
    FS = "\t"
    OFS = "\t"
}
{
    peak1_size = $3 - $2
    peak2_size = $12 - $11

    if (($19 / peak1_size >= 0.5) ||
        ($19 / peak2_size >= 0.5)) {
        print $0
    }
}' |
cut -f 1-9 |
sort |
uniq |
bedtools intersect -wo \
    -a stdin \
    -b "${REP2}" |
awk '
BEGIN {
    FS = "\t"
    OFS = "\t"
}
{
    peak1_size = $3 - $2
    peak2_size = $12 - $11

    if (($19 / peak1_size >= 0.5) ||
        ($19 / peak2_size >= 0.5)) {
        print $0
    }
}' |
cut -f 1-9 |
sort |
uniq
