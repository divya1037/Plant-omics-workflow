```r
############################################
# ATAC-seq Genomic Distribution Plot
# Arabidopsis thaliana TAIR10.1
#
# Input:
#   ANNOTATED_ALL_PEAKS.txt
#
# Output:
#   Genomic distribution pie chart
############################################

############################################
# Load annotation file
############################################

annot <- read.delim(
  "annotation/ANNOTATED_ALL_PEAKS.txt",
  sep = "\t",
  header = TRUE,
  quote = ""
)


############################################
# Extract annotation column
############################################

anno <- annot$Annotation


############################################
# Classify genomic regions
############################################

category <- rep("Other", length(anno))

category[grep(
  "^promoter-TSS",
  anno,
  ignore.case = TRUE
)] <- "Promoter-TSS"

category[grep(
  "^TTS",
  anno,
  ignore.case = TRUE
)] <- "TTS"

category[grep(
  "exon",
  anno,
  ignore.case = TRUE
)] <- "Exon"

category[grep(
  "intron",
  anno,
  ignore.case = TRUE
)] <- "Intron"

category[grep(
  "intergenic",
  anno,
  ignore.case = TRUE
)] <- "Intergenic"


############################################
# Count genomic categories
############################################

anno.count <- table(category)

print(anno.count)


############################################
# Calculate percentages
############################################

pct <- round(
  100 * anno.count / sum(anno.count),
  1
)

annotation.summary <- data.frame(
  Category = names(anno.count),
  Count = as.numeric(anno.count),
  Percent = pct
)

print(annotation.summary)


############################################
# Create labels
############################################

labels <- paste0(
  names(anno.count),
  "\n",
  pct,
  "%"
)


############################################
# Create output directory
############################################

dir.create(
  "results/figures",
  recursive = TRUE,
  showWarnings = FALSE
)


############################################
# Generate genomic distribution pie chart
############################################

png(
  "results/figures/genomic_distribution_ATAC_peaks.png",
  width = 2400,
  height = 2400,
  res = 300
)

pie(
  anno.count,
  labels = labels,
  col = c(
    "#E41A1C",
    "#377EB8",
    "#4DAF4A",
    "#984EA3",
    "#FF7F00"
  ),
  main = "Genomic Distribution of ATAC-seq Peaks"
)

dev.off()


############################################
# Save annotation summary
############################################

write.table(
  annotation.summary,
  "results/genomic_annotation_summary.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

############################################
# Completed
############################################

cat(
  "==========================================\n",
  "Genomic annotation plot completed.\n",
  "==========================================\n"
)
```
