# ATAC-seq Analysis: A Reproducible Workflow for Chromatin Accessibility Analysis

## Overview

This repository provides the scripts and supporting files used for the analysis of **ATAC-seq data** from *Arabidopsis thaliana* using the **TAIR10.1 reference genome**.

The workflow includes raw sequencing data acquisition, quality assessment, adapter trimming, read alignment, post-alignment processing, blacklist generation, peak calling, reproducibility assessment, differential accessibility analysis, peak annotation, motif enrichment analysis, and visualization.

---

## Workflow

```text
ENA/NCBI Sequencing Data
          │
          ▼
Quality Control (FastQC)
          │
          ▼
Quality Summary (MultiQC)
          │
          ▼
Adapter & Quality Trimming (Trim Galore)
          │
          ▼
Reference Genome (TAIR10.1)
          │
          ▼
Read Alignment (Bowtie2)
          │
          ▼
BAM Processing & Duplicate Removal
          │
          ▼
ATAC-seq Fragment Processing
          │
          ▼
Tn5 Shifting
          │
          ▼
Genomic Blacklist Generation
          │
          ▼
Broad Peak Calling (MACS2)
          │
          ▼
Blacklist Filtering & Reproducible Peaks
          │
          ▼
Differential Accessibility Analysis
        (csaw + edgeR)
          │
          ▼
Peak Annotation & Motif Analysis
        (HOMER)
          │
          ▼
Visualization
```

---

## Repository Structure

```text
ATAC-seq/
│
├── README.md
├── LICENSE
├── environment.yml
├── 07_install_R_packages.R
│
├── scripts/
│   ├── 01_download.sh
│   ├── 02_quality_trimming.sh
│   ├── 03_alignment.sh
│   ├── 04_post_alignment.sh
│   ├── 05_blacklist.sh
│   ├── 06_peak_calling.sh
│   ├── 08_CSAW_analysis.R
│   ├── 09_peak_annotation_motif.sh
│   ├── 10_annotation_plot.R
│   │
│   └── tools/
│       ├── bedpeMinimalConvert.sh
│       ├── bedpeTn5shift.sh
│       └── naiveOverlapBroad.sh
│
├── metadata/
│   └── sample.tsv
│
├── annotation/
│   └── ANNOTATED_ALL_PEAKS.txt
│
├── blacklist/
│   └── final_blacklist.bed
│
└── results/
```

---

## Case Study Dataset

The repository uses publicly available paired-end ATAC-seq datasets from *Arabidopsis thaliana*.

| Sample   | Condition | Replicate | SRA/ENA Accession |
| -------- | --------- | --------- | ----------------- |
| control1 | Control   | 1         | DRR826752         |
| control2 | Control   | 2         | DRR826753         |
| heat1    | Heat      | 1         | DRR826754         |
| heat2    | Heat      | 2         | DRR826755         |

The paired-end FASTQ files are downloaded by `01_download.sh` and renamed as:

```text
DRR826752_1.fastq.gz → control1_R1.fastq.gz
DRR826752_2.fastq.gz → control1_R2.fastq.gz

DRR826753_1.fastq.gz → control2_R1.fastq.gz
DRR826753_2.fastq.gz → control2_R2.fastq.gz

DRR826754_1.fastq.gz → heat1_R1.fastq.gz
DRR826754_2.fastq.gz → heat1_R2.fastq.gz

DRR826755_1.fastq.gz → heat2_R1.fastq.gz
DRR826755_2.fastq.gz → heat2_R2.fastq.gz
```

---

## Reference Genome

*Arabidopsis thaliana* TAIR10.1

NCBI Assembly: `GCF_000001735.4`

The reference genome and genome annotation are downloaded by `01_download.sh` and are not included in the repository.

---

## Software Requirements

| Software      | Purpose                                       |
| ------------- | --------------------------------------------- |
| FastQC        | Sequencing quality assessment                 |
| MultiQC       | Quality-control report aggregation            |
| Trim Galore   | Adapter and quality trimming                  |
| Bowtie2       | Read alignment                                |
| SAMtools      | BAM processing and statistics                 |
| BEDTools      | Genomic interval processing                   |
| Picard        | Read-group assignment and duplicate removal   |
| MACS2         | Broad peak calling                            |
| HOMER         | Peak annotation and motif enrichment analysis |
| R             | Downstream analysis and visualization         |
| csaw          | Differential accessibility analysis           |
| GenomicRanges | Genomic interval manipulation                 |
| Rsamtools     | BAM file processing                           |
| edgeR         | Statistical analysis                          |
| ggplot2       | Data visualization                            |
| pheatmap      | Sample correlation heatmap                    |

The command-line software environment is provided in `environment.yml`.

The required R packages are provided in `install_R_packages.R`.

---

## Installation

Clone the repository:

```bash
git clone https://github.com/<username>/ATAC-seq.git
cd ATAC-seq
```

Create the Conda environment:

```bash
conda env create -f environment.yml
conda activate atacseq_analysis
```

Install the required R packages:

```bash
Rscript install_R_packages.R
```

HOMER and the TAIR10 genome configuration are installed separately as described in `08_peak_annotation_motif.sh`.

---

## Pipeline Execution

### Step 1. Download sequencing data and reference genome

```bash
bash scripts/01_download.sh
```

### Step 2. Perform quality control and preprocessing

```bash
bash scripts/02_quality_trimming.sh
```

### Step 3. Perform read alignment

```bash
bash scripts/03_alignment.sh
```

### Step 4. Perform post-alignment processing

```bash
bash scripts/04_post_alignment.sh
```

### Step 5. Generate genomic blacklist

```bash
bash scripts/05_blacklist.sh
```

### Step 6. Perform peak calling and identify reproducible peaks

```bash
bash scripts/06_peak_calling.sh
```

### Step 7. Perform differential accessibility analysis

```bash
Rscript scripts/07_csaw_workflow.R
```

### Step 8. Perform peak annotation and motif analysis

```bash
bash scripts/08_peak_annotation_motif.sh
```

### Step 9. Generate genomic distribution plot

```bash
Rscript scripts/09_annotation_plot.R
```

---

## Expected Outputs

The pipeline generates:

* FastQC reports
* MultiQC reports
* Trimmed FASTQ files
* Bowtie2 alignment BAM files
* Sorted and indexed BAM files
* Alignment statistics
* Duplicate-removal metrics
* Filtered BAM files
* BEDPE fragment files
* Tn5-shifted BEDPE files
* MACS2-compatible BEDPE files
* Genomic blacklist
* MACS2 broad peak files
* Blacklist-filtered peak files
* Reproducible peak files
* Differential accessibility tables
* Fragment-size distribution plots
* MA plot
* Sample correlation heatmap
* Differentially accessible peak BED files
* Peak annotation results
* Upregulated and downregulated peak BED files
* HOMER motif enrichment results
* Genomic distribution pie chart

Representative results and figures can be placed in the **results/** directory.

---

## Reproducibility

The workflow was developed and tested in a **Linux/WSL2 environment**.

The numbered scripts should be executed sequentially as described in the pipeline execution section. Large sequencing files, reference genome files, alignment files, and intermediate files are not included in the repository.

The software environment is specified in `environment.yml`, and the required R packages are specified in `install_R_packages.R`.

---

## Contact

For questions regarding the workflow or the associated book chapter, please contact the corresponding author.
