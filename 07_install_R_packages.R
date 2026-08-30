############################################
# ATAC-seq R Package Installation
############################################

############################################
# CRAN packages
############################################

install.packages(c(
  "ggplot2",
  "pheatmap"
))


############################################
# Install BiocManager
############################################

install.packages("BiocManager")


############################################
# Bioconductor packages
############################################

BiocManager::install(c(
  "csaw",
  "GenomicRanges",
  "Rsamtools",
  "edgeR"
))