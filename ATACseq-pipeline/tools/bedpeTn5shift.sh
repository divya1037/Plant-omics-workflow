#!/bin/bash

############################################
# BEDPE Tn5 insertion-site correction
#
# Adapted from:
# Jake Reske, Michigan State University, 2019
#
# Purpose:
# Adjust ATAC-seq BEDPE fragment coordinates to
# account for the 9-bp offset introduced by Tn5.
#
# The strand-specific adjustment is:
#   + strand: left coordinate +4, right coordinate -5
#   - strand: left coordinate -5, right coordinate +4
#
# Usage:
# bash bedpeTn5shift.sh input.bedpe > output.tn5.bedpe
############################################

INPUT_BEDPE="$1"

if [[ -z "${INPUT_BEDPE}" ]]; then
    echo "Usage: bash bedpeTn5shift.sh input.bedpe" >&2
    exit 1
fi

awk -F $'\t' '
BEGIN {
    OFS = FS
}
{
    if ($9 == "+") {
        $2 = $2 + 4
        $6 = $6 - 5
    }
    else if ($9 == "-") {
        $3 = $3 - 5
        $5 = $5 + 4
    }

    print $0
}
' "${INPUT_BEDPE}"
