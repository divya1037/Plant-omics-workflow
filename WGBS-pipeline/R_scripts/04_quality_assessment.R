###############################################################
# WGBS Computational Protocol
# Script : 04_quality_assessment.R
#
# Description:
#   Performs quality assessment of Bismark methylation
#   extraction results.
#
#   - Automatically reads M-bias files
#   - Automatically parses Bismark splitting reports
#   - Generates QC plots
#   - Produces methylation summary tables
#
# Input
#   *.M-bias.txt
#   *_splitting_report.txt
#
# Output
#   QC summary table
#   M-bias plot
#   Genome-wide methylation composition
#
# Author:
#   Surapuram Aswini et al.
###############################################################

##############################
# Load Packages
##############################

required_packages <- c("data.table")

for(pkg in required_packages){

  if(!require(pkg, character.only=TRUE)){

    install.packages(pkg)

    library(pkg, character.only=TRUE)

  }

}

##############################
# Create Directories
##############################

dirs <- c("results","figures","logs")

for(d in dirs){

  if(!dir.exists(d))
    dir.create(d)

}

###############################################################
# Locate Files
###############################################################

mbias_files <- list.files(

pattern="M-bias.txt$",

full.names=TRUE

)

report_files <- list.files(

pattern="splitting_report.txt$",

full.names=TRUE

)

if(length(mbias_files)==0)
stop("No M-bias files found.")

if(length(report_files)==0)
stop("No Bismark splitting reports found.")

###############################################################
# Function to Plot M-bias
###############################################################

plot_mbias <- function(file,outfile){

lines <- readLines(file)

idx <- grep("^position",lines)

block1 <- read.table(
text=lines[(idx[1]+1):(idx[2]-2)],
header=FALSE
)

block2 <- read.table(
text=lines[(idx[2]+1):(idx[3]-2)],
header=FALSE
)

block3 <- read.table(
text=lines[(idx[3]+1):length(lines)],
header=FALSE
)

colnames(block1)<-
colnames(block2)<-
colnames(block3)<-
c(
"Position",
"Methylated",
"Unmethylated",
"Percent",
"Coverage"
)

png(
outfile,
width=1200,
height=900,
res=300
)

plot(
block1$Position,
block1$Percent,
type="l",
lwd=3,
col="steelblue",
ylim=c(0,max(c(
block1$Percent,
block2$Percent,
block3$Percent
))),
xlab="Read Position",
ylab="Methylation (%)",
main=basename(file)
)

lines(
block2$Position,
block2$Percent,
col="darkgreen",
lwd=3
)

lines(
block3$Position,
block3$Percent,
col="tomato",
lwd=3
)

legend(
"topright",
legend=c("CG","CHG","CHH"),
col=c("steelblue","darkgreen","tomato"),
lwd=3,
bty="n"
)

dev.off()

}

###############################################################
# Generate M-bias Plots
###############################################################

for(i in seq_along(mbias_files)){

outfile <- paste0(

"figures/Mbias_",

i,

".png"

)

plot_mbias(

mbias_files[i],

outfile

)

}

###############################################################
# Parse Bismark Reports
###############################################################

extract_percentage <- function(lines,pattern){

value <- grep(pattern,lines,value=TRUE)

value <- sub(".*:\\s*","",value)

value <- gsub("%","",value)

as.numeric(value)

}

summary_list <- list()

for(f in report_files){

lines <- readLines(f)

cg <- extract_percentage(
lines,
"C methylated in CpG context"
)

chg <- extract_percentage(
lines,
"C methylated in CHG context"
)

chh <- extract_percentage(
lines,
"C methylated in CHH context"
)

summary_list[[basename(f)]] <-

data.frame(

Sample=gsub(
"_splitting_report.txt",
"",
basename(f)
),

CG=cg,

CHG=chg,

CHH=chh

)

}

summary_table <-

rbindlist(summary_list)

###############################################################
# Save QC Summary
###############################################################

write.csv(

summary_table,

"results/Methylation_Summary.csv",

row.names=FALSE

)

###############################################################
# Plot Genome-wide Methylation
###############################################################

meth <- t(as.matrix(summary_table[,2:4]))

png(

"figures/Genomewide_Methylation.png",

width=1200,

height=900,

res=300

)

barplot(

meth,

beside=TRUE,

col=c(

"steelblue",

"darkgreen",

"tomato"

),

ylim=c(

0,

max(meth)+5

),

ylab="Methylation (%)",

main="Genome-wide Cytosine Methylation"

)

legend(

"topright",

legend=rownames(meth),

fill=c(

"steelblue",

"darkgreen",

"tomato"

),

bty="n"

)

dev.off()

###############################################################
# Save Session Information
###############################################################

writeLines(

capture.output(sessionInfo()),

"logs/sessionInfo_QC.txt"

)

###############################################################
# Finished
###############################################################

cat("\n=====================================\n")

cat("Quality Assessment Completed\n")

cat("=====================================\n")

cat("M-bias files processed :",length(mbias_files),"\n")

cat("Reports processed :",length(report_files),"\n")

cat("Results saved in results/\n")

cat("Figures saved in figures/\n")

cat("=====================================\n")