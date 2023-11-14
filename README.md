# A Pancreatic Cancer Organoid Platform Identifies an Inhibitor Specific to Mutant KRAS

This repository contains the R scripts necessary to perform the analysis in our
manuscript Xiaohua Duan et al. “A Pancreatic Cancer Organoid Platform Identifies an
Inhibitor Specific to Mutant KRAS.”, as described in the supplementary methods and
main text.

### Input Data

The single cell RNA-seq data were generated with the 10X kit and pre-processed
using the 10X cellranger pipeline. The raw data are available in the GEO
database with accession#
[GSE207352](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?&acc=GSE207352).

### Requirements

The following R packages were used:
- Seurat (v3.1.0)
- ggplot2
- pheatmap
- dplyr
- tibble

The following python package was used:
- pyscenic

