# Computational Epigenomics Protocols: A Reproducible Workflow for Whole-Genome Bisulfite Sequencing (WGBS) Analysis

## Overview

This repository accompanies the book chapter **"Computational Approaches in Epigenomics"** and provides a reproducible computational workflow for the analysis of **Whole-Genome Bisulfite Sequencing (WGBS)** data.

The workflow demonstrates a complete end-to-end WGBS analysis pipeline, including raw data acquisition, quality assessment, adapter trimming, bisulfite alignment, methylation extraction, quality evaluation, differential methylation analysis, and visualization using publicly available *Arabidopsis thaliana* datasets.

---

## Workflow

```
NCBI SRA Data
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
Reference Genome Preparation (Bismark)
      │
      ▼
Bisulfite Alignment
      │
      ▼
PCR Duplicate Removal
      │
      ▼
Methylation Extraction
      │
      ▼
Coverage Files & Cytosine Reports
      │
      ▼
Downstream Statistical Analysis (R)
      │
      ▼
Visualization & Quality Assessment
```

---

# Repository Structure

```
WGBS_Pipeline/

├── README.md
├── LICENSE
├── environment.yml
├── requirements.md
├── sample_metadata.csv
│
├── shell_scripts/
│   ├── 01_data_download.sh
│   ├── 02_preprocessing.sh
│   ├── 03_alignment_methylation.sh
│   └── 04_postprocessing_reporting.sh
│
├── R_scripts/
│   ├── 01_data_processing.R
│   ├── 02_visualization.R
│   ├── 03_differential_methylation.R
│   └── 04_quality_assessment.R
│
├── docs/
│
├── example_results/
│
├── figures/
│
└── reports/
```

---

# Case Study Dataset

The repository uses publicly available **Arabidopsis thaliana** WGBS datasets.

| Sample | SRA Accession | Condition | GEO Accession |
|---------|---------------|-----------|---------------|
| Control | SRR5494752 | Unstressed (4-week-old Col-0) | GSM2595563 |
| Treatment | SRR5494755 | Nine days water deprivation (4-week-old Col-0) | GSM2595566 |

Additional sample information is provided in **sample_metadata.csv**.

---

# Software Requirements

| Software     |              Version | Purpose                                            |
| ------------ | -------------------: | -------------------------------------------------- |
| SRA Toolkit  | Version not recorded | Retrieval of sequencing data from NCBI SRA         |
| FastQC       |               0.11.9 | Raw sequencing quality assessment                  |
| Trim Galore! |               0.6.11 | Adapter and quality trimming                       |
| Cutadapt     |                  3.5 | Adapter removal                                    |
| Bowtie2      |                2.4.4 | Read alignment engine used by Bismark              |
| Bismark      |              0.25.1 | Bisulfite-aware alignment and methylation analysis |
| SAMtools     |                 1.13 | BAM/SAM processing                                 |
| R            |                4.5.2 | Downstream analysis and visualization              |
| Bioconductor |                 3.22 | R-based genomic analysis framework                 |
| data.table   |             1.18.2.1 | Methylation data processing                        |
| ggplot2      |                4.0.2 | Data visualization                                 |
| dplyr        |                1.1.4 | Data manipulation                                  |
| readr        |                2.2.0 | Data import                                        |
| edgeR        |                4.8.2 | Statistical analysis                               |

The complete software environment is available in **environment.yml**.

---

# Installation

Clone the repository:

```bash
git clone https://github.com/<username>/WGBS_Pipeline.git
cd WGBS_Pipeline
```

Create the Conda environment:

```bash
conda env create -f environment.yml
conda activate wgbs_pipeline
```

---

# Pipeline Execution

### Step 1. Download sequencing data

```bash
bash shell_scripts/01_data_download.sh
```

### Step 2. Perform quality control and preprocessing

```bash
bash shell_scripts/02_preprocessing.sh
```

### Step 3. Perform bisulfite alignment and methylation calling

```bash
bash shell_scripts/03_alignment_methylation.sh
```

### Step 4. Generate reports and organize outputs

```bash
bash shell_scripts/04_postprocessing_reporting.sh
```

### Step 5. Process methylation data

```bash
Rscript R_scripts/01_data_processing.R
```

### Step 6. Generate visualizations

```bash
Rscript R_scripts/02_visualization.R
```

### Step 7. Perform differential methylation analysis

```bash
Rscript R_scripts/03_differential_methylation.R
```

### Step 8. Perform quality assessment

```bash
Rscript R_scripts/04_quality_assessment.R
```

---

# Expected Outputs

The pipeline generates:

- FastQC reports
- MultiQC reports
- Trimmed FASTQ files
- BAM alignment files
- Deduplicated BAM files
- Cytosine methylation reports
- Coverage (.cov.gz) files
- BedGraph files
- M-bias reports
- Genome-wide methylation summaries
- Differential methylation tables
- Publication-quality figures

Representative outputs are available in the **example_results/** directory.

---

# Reproducibility

The workflow has been tested on **Ubuntu Linux (WSL2)** using the software versions listed in `environment.yml`. All scripts are modular and can be executed independently or as part of the complete pipeline.

---

# Citation

If you use this workflow in your research, please cite the associated book chapter:

**Aswini S., et al. Computational Approaches in Epigenomics. Methods in Molecular Biology.**

Please also cite the original software tools and the public datasets used in the case study.

---

# Contact

For questions regarding the workflow or the associated book chapter, please contact the corresponding author.