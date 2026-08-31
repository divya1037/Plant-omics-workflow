
############################################
# CSAW Differential Accessibility Analysis
# Arabidopsis thaliana TAIR10.1
############################################

############################################
# Load libraries
############################################

library(GenomicRanges)
library(csaw)
library(Rsamtools)
library(edgeR)
library(ggplot2)
library(pheatmap)

############################################
# Directory paths
############################################

peak.dir <- "alignment/filter/atac_format/peakcalling"
bam.dir <- "alignment/filter"
blacklist.file <- "blacklist/final_blacklist.bed"

############################################
# Read MACS2 filtered broadPeak files
############################################

treat1.peaks <- read.table(
  file.path(peak.dir, "heat1_peaks.filt.broadPeak"),
  sep = "\t"
)[,1:3]

treat2.peaks <- read.table(
  file.path(peak.dir, "heat2_peaks.filt.broadPeak"),
  sep = "\t"
)[,1:3]

control1.peaks <- read.table(
  file.path(peak.dir, "control1_peaks.filt.broadPeak"),
  sep = "\t"
)[,1:3]

control2.peaks <- read.table(
  file.path(peak.dir, "control2_peaks.filt.broadPeak"),
  sep = "\t"
)[,1:3]

colnames(treat1.peaks) <- c("chrom", "start", "end")
colnames(treat2.peaks) <- c("chrom", "start", "end")
colnames(control1.peaks) <- c("chrom", "start", "end")
colnames(control2.peaks) <- c("chrom", "start", "end")

############################################
# Read reproducible broadPeak files
############################################

treat.overlap.peaks <- read.table(
  file.path(peak.dir, "heat_overlap_peaks.filt.broadPeak"),
  sep = "\t"
)[,1:3]

control.overlap.peaks <- read.table(
  file.path(peak.dir, "control_overlap_peaks.filt.broadPeak"),
  sep = "\t"
)[,1:3]

colnames(treat.overlap.peaks) <- c("chrom", "start", "end")
colnames(control.overlap.peaks) <- c("chrom", "start", "end")

############################################
# Convert to GRanges objects
############################################

treat1.peaks <- GRanges(treat1.peaks)
treat2.peaks <- GRanges(treat2.peaks)

control1.peaks <- GRanges(control1.peaks)
control2.peaks <- GRanges(control2.peaks)

treat.overlap.peaks <- GRanges(treat.overlap.peaks)
control.overlap.peaks <- GRanges(control.overlap.peaks)

############################################
# Define consensus peakset
#
# Union of reproducible heat and control peaks
############################################

all.peaks <- union(
  treat.overlap.peaks,
  control.overlap.peaks
)

############################################
# Specify paired-end BAM files
############################################

pe.bams <- c(
  file.path(bam.dir, "control1.fil.sorted.bam"),
  file.path(bam.dir, "control2.fil.sorted.bam"),
  file.path(bam.dir, "heat1.fil.sorted.bam"),
  file.path(bam.dir, "heat2.fil.sorted.bam")
)

############################################
# Read Arabidopsis TAIR10.1 blacklist
############################################

blacklist <- read.table(
  blacklist.file,
  sep = "\t"
)

colnames(blacklist) <- c(
  "chrom",
  "start",
  "end"
)

blacklist <- GRanges(blacklist)

############################################
# Define read parameters
############################################

standard.chr <- seqlevels(
  BamFile(
    file.path(
      bam.dir,
      "control1.fil.sorted.bam"
    )
  )
)

param <- readParam(
  max.frag = 1000,
  pe = "both",
  discard = blacklist,
  restrict = standard.chr
)

############################################
# Count reads in consensus peaks
############################################

peak.counts <- regionCounts(
  pe.bams,
  all.peaks,
  param = param
)

############################################
# Filter low-abundance peaks
############################################

peak.abundances <- aveLogCPM(
  asDGEList(peak.counts)
)

peak.counts.filt <- peak.counts[
  peak.abundances > -3,
]

############################################
# Paired-end fragment size distribution
############################################

control1.pe.sizes <- getPESizes(
  file.path(bam.dir, "control1.fil.sorted.bam")
)

control2.pe.sizes <- getPESizes(
  file.path(bam.dir, "control2.fil.sorted.bam")
)

treat1.pe.sizes <- getPESizes(
  file.path(bam.dir, "heat1.fil.sorted.bam")
)

treat2.pe.sizes <- getPESizes(
  file.path(bam.dir, "heat2.fil.sorted.bam")
)

gc()

############################################
# Fragment size distribution plots
############################################

hist(
  treat1.pe.sizes$sizes,
  col = "red",
  main = "Heat1 Fragment Size Distribution",
  xlab = "Fragment Size (bp)"
)

hist(
  treat2.pe.sizes$sizes,
  col = "orange",
  main = "Heat2 Fragment Size Distribution",
  xlab = "Fragment Size (bp)"
)

hist(
  control1.pe.sizes$sizes,
  col = "blue",
  main = "Control1 Fragment Size Distribution",
  xlab = "Fragment Size (bp)"
)

hist(
  control2.pe.sizes$sizes,
  col = "darkgreen",
  main = "Control2 Fragment Size Distribution",
  xlab = "Fragment Size (bp)"
)

############################################
# Count BAM reads in 150 bp windows
############################################

counts <- windowCounts(
  pe.bams,
  width = 150,
  param = param
)

############################################
# Filter windows by local enrichment
#
# Local background estimator: 2 kb
############################################

neighbor <- suppressWarnings(
  resize(
    rowRanges(counts),
    width = 2000,
    fix = "center"
  )
)

wider <- regionCounts(
  pe.bams,
  regions = neighbor,
  param = param
)

filter.stat <- filterWindowsLocal(
  counts,
  wider
)

counts.local.filt <- counts[
  filter.stat$filter > log2(3),
]

############################################
# Count background bins for TMM normalization
############################################

binned <- windowCounts(
  pe.bams,
  bin = TRUE,
  width = 10000,
  param = param
)

############################################
# NORMALIZATION
############################################

peak.counts.tmm <- peak.counts.filt

peak.counts.tmm <- normFactors(
  binned,
  se.out = peak.counts.tmm
)

############################################
# DIFFERENTIAL ACCESSIBILITY ANALYSIS
############################################

working.windows <- peak.counts.tmm

############################################
# Setup edgeR object
############################################

y <- asDGEList(
  working.windows
)

colnames(y$counts) <- c(
  "control1",
  "control2",
  "treat1",
  "treat2"
)

rownames(y$samples) <- c(
  "control1",
  "control2",
  "treat1",
  "treat2"
)

y$samples$group <- c(
  "control",
  "control",
  "treat",
  "treat"
)

############################################
# Design matrix
############################################

design <- model.matrix(
  ~0 + group,
  data = y$samples
)

colnames(design) <- c(
  "control",
  "treat"
)

############################################
# Dispersion estimation
############################################

y <- estimateDisp(
  y,
  design
)

############################################
# Fit quasi-likelihood model
############################################

fit <- glmQLFit(
  y,
  design,
  robust = TRUE
)

############################################
# Test for differential accessibility
############################################

results <- glmQLFTest(
  fit,
  contrast = makeContrasts(
    treat - control,
    levels = design
  )
)

############################################
# Add differential statistics
############################################

rowData(working.windows) <- cbind(
  rowData(working.windows),
  results$table
)

############################################
# Merge nearby windows
#
# Tolerance: 500 bp
# Maximum merged width: 5000 bp
############################################

merged.peaks <- mergeWindows(
  rowRanges(working.windows),
  tol = 500L,
  max.width = 5000L
)

############################################
# Select most significant window
# for each merged region
############################################

tab.best <- getBestTest(
  merged.peaks$id,
  results$table
)

############################################
# Create final merged peak table
############################################

final.merged.peaks <- GRanges(
  cbind(
    as.data.frame(merged.peaks$region),
    results$table[
      tab.best$rep.test,
      -4
    ],
    tab.best[,-c(7:8)]
  )
)

############################################
# Sort by FDR
############################################

final.merged.peaks <- final.merged.peaks[
  order(
    final.merged.peaks@elementMetadata$FDR
  ),
]

final.merged.peaks

############################################
# Filter significant DA regions
############################################

FDR.thresh <- 0.05

final.merged.peaks.sig <- final.merged.peaks[
  final.merged.peaks@elementMetadata$FDR < FDR.thresh,
]

final.merged.peaks.sig

############################################
# Write CSAW results
############################################

write.table(
  final.merged.peaks,
  "treat_vs_control_csaw_DA-windows_all.txt",
  sep = "\t",
  quote = F,
  col.names = T,
  row.names = F
)

write.table(
  final.merged.peaks.sig,
  "treat_vs_control_csaw_DA-windows_significant.txt",
  sep = "\t",
  quote = F,
  col.names = T,
  row.names = F
)

############################################
# MA plot
############################################

final.merged.peaks$sig <- "n.s."

final.merged.peaks$sig[
  final.merged.peaks$FDR < FDR.thresh
] <- "significant"

ggplot(
  data = data.frame(final.merged.peaks),
  aes(
    x = logCPM,
    y = logFC,
    col = factor(
      sig,
      levels = c(
        "n.s.",
        "significant"
      )
    )
  )
) +
  geom_point() +
  scale_color_manual(
    values = c(
      "black",
      "red"
    )
  ) +
  geom_smooth(
    inherit.aes = F,
    aes(
      x = logCPM,
      y = logFC
    ),
    method = "loess"
  ) +
  geom_hline(
    yintercept = 0
  ) +
  labs(
    col = NULL
  )

############################################
# SAMPLE CORRELATION HEATMAP
############################################

logcpm <- cpm(
  asDGEList(peak.counts),
  log = TRUE,
  prior.count = 1
)

############################################
# Calculate Pearson correlation
############################################

sample.cor <- cor(
  logcpm,
  method = "pearson"
)

############################################
# Rename samples
############################################

colnames(sample.cor) <- c(
  "Control_1",
  "Control_2",
  "Heat_1",
  "Heat_2"
)

rownames(sample.cor) <- colnames(sample.cor)

############################################
# Plot correlation heatmap
############################################

pheatmap(
  sample.cor,
  display_numbers = TRUE,
  main = "Pearson correlation heatmap"
)

############################################
# End of CSAW workflow
############################################