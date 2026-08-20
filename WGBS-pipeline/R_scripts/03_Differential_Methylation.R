###############################################################
# WGBS Computational Protocol
# Script : 03_differential_methylation.R
#
# Description:
#   Identifies putative differentially methylated cytosines
#   (DMCs) based on methylation difference threshold.
#
# NOTE:
#   Demonstration workflow using one representative sample
#   per condition. Statistical DMR calling requires
#   biological replicates (e.g., DSS, methylKit).
#
# Author:
#   Surapuram Aswini et al.
###############################################################

##############################
# Load Packages
##############################

required_packages <- c("data.table")

for(pkg in required_packages){

  if(!require(pkg,character.only=TRUE)){

    install.packages(pkg)

    library(pkg,character.only=TRUE)

  }

}

##############################
# Parameters
##############################

difference_threshold <- 0.30

##############################
# Create Output Directories
##############################

dirs <- c("results","figures","logs")

for(d in dirs){

  if(!dir.exists(d))
    dir.create(d)

}

##############################
# Load Processed Data
##############################

if(!file.exists("results/processed_methylation_data.RData")){

  stop("Run 01_data_processing.R first.")

}

load("results/processed_methylation_data.RData")

##############################
# Identify Putative DMCs
##############################

putative_DMCs <- merged_data[
  abs(Difference) >= difference_threshold
]

##############################
# Classify DMCs
##############################

putative_DMCs$Status <- ifelse(

  putative_DMCs$Difference > 0,

  "Hypermethylated",

  "Hypomethylated"

)

##############################
# Separate Hyper/Hypo Sites
##############################

hyper_DMCs <- putative_DMCs[
  Status=="Hypermethylated"
]

hypo_DMCs <- putative_DMCs[
  Status=="Hypomethylated"
]

##############################
# Save Results
##############################

fwrite(

putative_DMCs,

"results/Putative_DMCs.tsv",

sep="\t"

)

fwrite(

hyper_DMCs,

"results/Hypermethylated_DMCs.tsv",

sep="\t"

)

fwrite(

hypo_DMCs,

"results/Hypomethylated_DMCs.tsv",

sep="\t"

)

##############################
# Summary Statistics
##############################

summary_df <- data.frame(

Metric=c(

"Total Cytosines",

"Putative DMCs",

"Hypermethylated",

"Hypomethylated",

"Threshold"

),

Value=c(

nrow(merged_data),

nrow(putative_DMCs),

nrow(hyper_DMCs),

nrow(hypo_DMCs),

difference_threshold

)

)

write.csv(

summary_df,

"results/DMC_summary.csv",

row.names=FALSE

)

##############################
# Volcano-style Plot
##############################

png(

"figures/Figure5_Differential_Methylation.png",

width=1200,

height=900,

res=300

)

plot(

merged_data$Control,

merged_data$Difference,

pch=16,

cex=0.4,

col=rgb(0,0,0,0.15),

xlab="Control Methylation Level",

ylab="Methylation Difference",

main="Putative Differentially Methylated Cytosines"

)

abline(

h=c(-difference_threshold,difference_threshold),

col="red",

lty=2,

lwd=2

)

dev.off()

##############################
# Pie Chart
##############################

png(

"figures/Figure6_DMC_Composition.png",

width=900,

height=900,

res=300

)

pie(

c(

nrow(hyper_DMCs),

nrow(hypo_DMCs)

),

labels=c(

paste0("Hyper\n",nrow(hyper_DMCs)),

paste0("Hypo\n",nrow(hypo_DMCs))

),

col=c("tomato","steelblue"),

main="Putative DMC Composition"

)

dev.off()

##############################
# Save Session Info
##############################

writeLines(

capture.output(sessionInfo()),

"logs/sessionInfo_DMC.txt"

)

##############################
# Finished
##############################

cat("\n=====================================\n")
cat("Differential Methylation Completed\n")
cat("=====================================\n")
cat("Threshold :",difference_threshold,"\n")
cat("Total DMCs :",nrow(putative_DMCs),"\n")
cat("Results saved in results/\n")
cat("Figures saved in figures/\n")
cat("=====================================\n")