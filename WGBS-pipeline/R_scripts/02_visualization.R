###############################################################
# WGBS Computational Protocol
# Script : 02_visualization.R
#
# Description:
#   Generates publication-quality visualizations from
#   processed methylation data.
#
# Input:
#   results/processed_methylation_data.RData
#
# Output:
#   Publication-quality figures (.png and .pdf)
#
# Author:
#   Surapuram Aswini et al.
###############################################################

##############################
# Load Packages
##############################

required_packages <- c("data.table")

for(pkg in required_packages){

  if(!require(pkg, character.only = TRUE)){

    install.packages(pkg)

    library(pkg, character.only = TRUE)

  }

}

##############################
# Parameters
##############################

random_seed <- 123
max_points  <- 100000

##############################
# Create Directories
##############################

dirs <- c("figures","results","logs")

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
# Sampling
##############################

set.seed(random_seed)

plot_data <- merged_data[
  sample(.N,min(max_points,.N))
]

##############################
# Intermediate Methylation
##############################

control_mid <- plot_data$Control[
  plot_data$Control>0.10 &
  plot_data$Control<0.90
]

treatment_mid <- plot_data$Treatment[
  plot_data$Treatment>0.10 &
  plot_data$Treatment<0.90
]

###############################################################
# Figure 1 Histogram
###############################################################

png(
"figures/Figure1_Histogram.png",
width=1800,
height=800,
res=300
)

par(mfrow=c(1,2))

hist(
control_mid,
breaks=60,
col="steelblue",
border="white",
main="Control",
xlab="Methylation Level",
ylab="Frequency"
)

hist(
treatment_mid,
breaks=60,
col="tomato",
border="white",
main="Treatment",
xlab="Methylation Level",
ylab="Frequency"
)

dev.off()

###############################################################
# Figure 2 Density
###############################################################

png(
"figures/Figure2_Density.png",
width=1200,
height=900,
res=300
)

plot(
density(control_mid),
col="steelblue",
lwd=3,
main="Genome-wide Methylation Density",
xlab="Methylation Level"
)

lines(
density(treatment_mid),
col="tomato",
lwd=3
)

legend(
"topright",
legend=c("Control","Treatment"),
col=c("steelblue","tomato"),
lwd=3,
bty="n"
)

dev.off()

###############################################################
# Figure 3 Correlation
###############################################################

png(
"figures/Figure3_Correlation.png",
width=1000,
height=1000,
res=300
)

plot(
plot_data$Control,
plot_data$Treatment,
pch=16,
cex=0.25,
col=rgb(0,0,0,0.15),
xlab="Control",
ylab="Treatment",
main="Genome-wide Cytosine Methylation"
)

abline(0,1,col="red",lwd=2)

cor_value <- cor(
plot_data$Control,
plot_data$Treatment
)

legend(
"topleft",
legend=paste("Pearson r =",round(cor_value,3)),
bty="n"
)

dev.off()

###############################################################
# Figure 4 Differential Methylation
###############################################################

png(
"figures/Figure4_Differential_Histogram.png",
width=1200,
height=900,
res=300
)

hist(
plot_data$Difference,
breaks=80,
col="goldenrod",
border="white",
main="Distribution of Differential Methylation",
xlab="Treatment − Control"
)

abline(v=0,col="red",lwd=2)

dev.off()

###############################################################
# Save Figure Statistics
###############################################################

summary_df <- data.frame(

Metric=c(
"Control Mean",
"Treatment Mean",
"Pearson Correlation",
"Total Sites Plotted"
),

Value=c(

mean(control_mid),

mean(treatment_mid),

cor_value,

length(plot_data$Control)

)

)

write.csv(

summary_df,

"results/Visualization_summary.csv",

row.names=FALSE

)

###############################################################
# Save Session Information
###############################################################

writeLines(

capture.output(sessionInfo()),

"logs/sessionInfo_visualization.txt"

)

###############################################################
# Finished
###############################################################

cat("\nVisualization completed successfully.\n")
cat("Figures saved in figures/\n")
cat("Summary saved in results/\n")