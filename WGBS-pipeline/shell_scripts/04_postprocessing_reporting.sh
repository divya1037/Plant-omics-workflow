#!/bin/bash

###############################################################
# WGBS Computational Protocol
# Script: 04_postprocessing_reporting.sh
#
# Description:
#   Organizes methylation analysis results, generates Bismark
#   summary reports, and prepares files for downstream
#   statistical analysis and visualization.
#
# Input:
#   results/
#
# Output:
#   reports/
#   example_results/
#
# Software:
#   Bismark v0.25.1
#
# Usage:
#   bash shell_scripts/04_postprocessing_reporting.sh
#
# Author:
#   Surapuram Aswini et al.
###############################################################

set -e
set -o pipefail

echo "=================================================="
echo "Step 4 : Post-processing and Reporting"
echo "Started : $(date)"
echo "=================================================="

###############################
# Check results directory
###############################

if [ ! -d "results" ]; then
    echo "ERROR: results directory not found."
    exit 1
fi

###############################
# Create output directories
###############################

mkdir -p reports
mkdir -p example_results
mkdir -p example_results/QC
mkdir -p example_results/Figures
mkdir -p example_results/Tables

###############################
# Generate Bismark HTML Report
###############################

echo ""
echo "Generating Bismark report..."

bismark2report

###############################
# Generate Bismark Summary
###############################

echo ""
echo "Generating Bismark summary..."

bismark2summary

###############################
# Organize Reports
###############################

echo ""
echo "Moving reports..."

mv *.html reports/ 2>/dev/null || true
mv *.txt reports/ 2>/dev/null || true

###############################
# Copy Important Results
###############################

echo ""
echo "Preparing example results..."

cp results/*.cov.gz example_results/ 2>/dev/null || true

cp results/*.bedGraph.gz example_results/ 2>/dev/null || true

cp results/*splitting_report.txt example_results/QC/ 2>/dev/null || true

cp results/*M-bias.txt example_results/QC/ 2>/dev/null || true

###############################
# Display Generated Files
###############################

echo ""
echo "Generated Reports"

ls reports

echo ""

echo "QC Files"

ls example_results/QC

###############################
# Finish
###############################

echo ""
echo "=================================================="
echo "Pipeline Completed Successfully"
echo "Finished : $(date)"
echo "=================================================="

echo ""
echo "Next Step:"
echo "Run the R scripts for downstream analysis:"
echo ""
echo "Rscript R_scripts/01_data_processing.R"
echo "Rscript R_scripts/02_visualization.R"
echo "Rscript R_scripts/03_differential_methylation.R"
echo "Rscript R_scripts/04_quality_assessment.R"
echo ""