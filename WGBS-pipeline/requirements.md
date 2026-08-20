# System Requirements

## Operating System

The WGBS workflow was performed in a Linux-based environment using WSL (Ubuntu).

> The exact Ubuntu version used for the original analysis should be reported if available.

## Hardware Requirements

| Component | Recommended |
|-----------|-------------|
| CPU | Quad-core or higher |
| RAM | Minimum 16 GB (32 GB recommended) |
| Storage | ≥100 GB free space |

## Software Requirements

| Software | Version |
|----------|---------|
| SRA Toolkit | Version not recorded |
| FastQC | 0.11.9 |
| Trim Galore! | 0.6.11 |
| Cutadapt | 3.5 |
| Bowtie2 | 2.4.4 |
| Bismark | 0.25.1* |
| SAMtools | 1.13 |
| R | 4.5.2 |

\* Bismark v0.25.1 was confirmed from the methylation-extraction output. The corresponding Bismark alignment version should be reported as v0.25.1 only if confirmed from the alignment/log files.

## R Environment

| Component | Version |
|-----------|---------|
| R | 4.5.2 |
| Bioconductor | 3.22 |
| data.table | 1.18.2.1 |
| ggplot2 | 4.0.2 |
| dplyr | 1.1.4 |
| readr | 2.2.0 |
| edgeR | 4.8.2 |

## R Packages Used

- data.table
- ggplot2
- dplyr
- readr

`edgeR` was available in the R environment but should only be listed as a workflow dependency if it is explicitly used in the analysis scripts.

## Installation

The analysis environment can be recreated using the supplied `environment.yml` file:

```bash
conda env create -f environment.yml
conda activate wgbs_pipeline