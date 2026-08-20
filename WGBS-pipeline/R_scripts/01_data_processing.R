###############################################################
# WGBS Computational Protocol
# Script : 01_data_processing.R
#
# Description:
#   Loads Bismark coverage files, performs quality filtering,
#   calculates methylation levels, merges samples and saves
#   processed methylation data.
#
# Input:
#   *.bismark.cov.gz
#
# Output:
#   processed_methylation_data.RData
#   processed_methylation_data.tsv
#
# Author:
#   Surapuram Aswini et al.
###############################################################

###############################################################
# Load Required Packages
###############################################################

required_packages <- c("data.table")

for(pkg in required_packages){

  if(!require(pkg, character.only = TRUE)){

    install.packages(pkg)

    library(pkg, character.only = TRUE)

  }

}

###############################################################
# Parameters
###############################################################

coverage_threshold <- 5

###############################################################
# Create Output Directories
###############################################################

dirs <- c("results","figures","logs")

for(d in dirs){

  if(!dir.exists(d))
    dir.create(d)

}

###############################################################
# Start Log
###############################################################

cat("\n=========================================\n")
cat("WGBS Data Processing\n")
cat("Started :", as.character(Sys.time()), "\n")
cat("=========================================\n")

###############################################################
# Locate Coverage Files
###############################################################

cov_files <- list.files(

  pattern="\\.bismark\\.cov\\.gz$",

  full.names=TRUE

)

if(length(cov_files) < 2){

  stop(
    "At least two Bismark coverage files (*.bismark.cov.gz) are required."
  )

}

cat("\nCoverage files detected:\n")

print(basename(cov_files))

###############################################################
# Load Coverage Files
###############################################################

control <- fread(cov_files[1])

treatment <- fread(cov_files[2])

###############################################################
# Rename Columns
###############################################################

column_names <- c(

"Chromosome",

"Start",

"End",

"Methylation_Percentage",

"Methylated_Reads",

"Unmethylated_Reads"

)

setnames(control,column_names)

setnames(treatment,column_names)

###############################################################
# Dataset Summary
###############################################################

cat("\nControl dimensions : ",dim(control)[1]," x ",dim(control)[2],"\n")

cat("Treatment dimensions : ",dim(treatment)[1]," x ",dim(treatment)[2],"\n")

###############################################################
# Coverage Filtering
###############################################################

control <- control[

(Methylated_Reads + Unmethylated_Reads) >= coverage_threshold

]

treatment <- treatment[

(Methylated_Reads + Unmethylated_Reads) >= coverage_threshold

]

###############################################################
# Calculate Methylation Level
###############################################################

control[,Methylation_Level:=

Methylated_Reads/

(Methylated_Reads+Unmethylated_Reads)

]

treatment[,Methylation_Level:=

Methylated_Reads/

(Methylated_Reads+Unmethylated_Reads)

]

###############################################################
# Remove Missing Values
###############################################################

control <- na.omit(control)

treatment <- na.omit(treatment)

###############################################################
# Merge Samples
###############################################################

merged_data <- merge(

control[,.(Chromosome,

Start,

Control=Methylation_Level)],

treatment[,.(Chromosome,

Start,

Treatment=Methylation_Level)],

by=c("Chromosome","Start")

)

###############################################################
# Differential Methylation
###############################################################

merged_data[,Difference:=Treatment-Control]

###############################################################
# Save Outputs
###############################################################

save(

control,

treatment,

merged_data,

file="results/processed_methylation_data.RData"

)

fwrite(

merged_data,

"results/processed_methylation_data.tsv",

sep="\t"

)

###############################################################
# Save Session Information
###############################################################

writeLines(

capture.output(sessionInfo()),

"logs/sessionInfo.txt"

)

###############################################################
# Completion Message
###############################################################

cat("\n-----------------------------------------\n")

cat("Processing Completed Successfully\n")

cat("-----------------------------------------\n")

cat("Coverage threshold :",coverage_threshold,"\n")

cat("Control cytosines :",nrow(control),"\n")

cat("Treatment cytosines :",nrow(treatment),"\n")

cat("Merged cytosines :",nrow(merged_data),"\n")

cat("\nOutput Files\n")

cat("results/processed_methylation_data.RData\n")

cat("results/processed_methylation_data.tsv\n")

cat("logs/sessionInfo.txt\n")

cat("\nFinished :",as.character(Sys.time()),"\n")

cat("-----------------------------------------\n")