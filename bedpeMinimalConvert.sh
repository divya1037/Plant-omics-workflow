#!/bin/bash

############################################
# BEDPE conversion for MACS2
#
# Adapted from:
# Jake Reske, Michigan State University, 2019
#
# Purpose:
# Convert the standard BEDPE representation produced
# by bedtools into the minimal BEDPE representation
# required for MACS2 paired-end peak calling.
#
# Usage:
# bash bedpeMinimalConvert.sh input.bedpe > output.minimal.bedpe
############################################

INPUT_BEDPE="$1"

if [[ -z "${INPUT_BEDPE}" ]]; then
    echo "Usage: bash bedpeMinimalConvert.sh input.bedpe" >&2
    exit 1
fi


############################################
# Sort paired fragment coordinates
############################################

coordSort()
{
    while read -r line
    do
        ary=(${line})

        for ((i=0; i<${#ary[@]}; ++i))
        do
            for ((j=i+1; j<${#ary[@]}; ++j))
            do
                if (( ary[i] > ary[j] ))
                then
                    key=${ary[i]}
                    ary[i]=${ary[j]}
                    ary[j]=${key}
                fi
            done
        done

        echo "${ary[@]}"

    done < "$1"
}


############################################
# Extract fragment coordinates
############################################

COORD_FILE="coords.$(basename "${INPUT_BEDPE}")"

awk 'BEGIN{OFS="\t"}
{
    printf "%s\t%s\t%s\t%s\n",$2,$3,$5,$6
}' "${INPUT_BEDPE}" > "${COORD_FILE}"


############################################
# Generate MACS2-compatible BEDPE
############################################

coordSort "${COORD_FILE}" |
paste - "${INPUT_BEDPE}" |
awk 'BEGIN{OFS="\t"}
{
    printf "%s\t%s\t%s\t%s\n",$5,$1,$4,$11
}' -


############################################
# Remove temporary coordinate file
############################################

rm -f "${COORD_FILE}"
