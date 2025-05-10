# DNA Methylation and Intron Polyadenylation (IPA) QTM Analysis

This project aims to explore the relationship between **DNA methylation** and **intron polyadenylation (IPA)** across human cancers. Specifically, it identifies **methylation sites (CpGs)** that are quantitatively associated with changes in IPA usage, also known as **IPA-QTMs (intron polyadenylation quantitative trait methylations)**.

## Project Structure

.
├── eQTM/ # Scripts to compute methylation sites associated with gene expression (eQTMs)
├── spQTM/ # Scripts to compute methylation sites associated with splicing (spQTMs)
├── figures/ # Final figures and summary outputs of the current analysis
├── *.r # Core pipeline scripts for IPA-QTM analysis
└── README.md # Project documentatio

## Description

- **IPA-QTM analysis** identifies CpG sites whose methylation levels are associated with changes in intronic PAS (polyadenylation site) usage.
- **eQTM analysis** (in `eQTM/`) finds CpGs linked to gene expression variation.
- **spQTM analysis** (in `spQTM/`) finds CpGs linked to alternative splicing variation.
- **Figures and summary files** are stored in the `figures/` folder.

## Methods Overview

1. **Prepare data**: PDUI matrix, methylation matrix, and sample metadata.
2. **Run IPA-QTM analysis**: Scripts in the root directory implement the regression framework for identifying IPA-associated CpGs.
3. **Run eQTM and spQTM analysis**: Use respective directories to run expression- or splicing-based QTM discovery.
4. **Visualize results**: Boxplots, significance annotations, and correlation heatmaps are saved in the `figures/` directory.

## Software Requirements

Ensure the following R packages are installed:

```r
install.packages(c("tidyverse", "ggpubr"))
BiocManager::install(c("limma", "minfi", "DESeq2", "sva", "biomaRt"))
n
