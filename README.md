# strpip: Data Structure Pipelines for Genomic Analysis
[![R-CMD-check](https://github.com/hungms/strpip/workflows/R-CMD-check/badge.svg)](https://github.com/hungms/strpip/actions)
[![pkgdown](https://github.com/hungms/strpip/workflows/pkgdown/badge.svg)](https://github.com/hungms/strpip/actions)

A comprehensive R package for genomic data structure pipelines, providing tools for data format conversions, gene symbol mapping, gene set operations, expression data processing, and functional annotation integration.

## Documentation

For detailed documentation and examples, please visit our latest [documentation](https://hungms.github.io/strpip/).


## Features

### Gene Symbol Mapping
- Convert between human and mouse gene symbols
- Map between different gene ID types

### Specialized Gene Set Operations
- Identify X/Y chromosome genes
- Find mitochondrial genes
- Filter structural genes
- Perform set operations on gene lists

### Expression Processing
- Summarize gene expression data
- Filter and normalize expression matrices
- Handle duplicate gene entries

### Functional Annotation
- Integrate pathway annotations
- Access OmniPath resources
- Perform enrichment analyses

### Data Format Operations
- Read and write GCT (Gene Cluster Text) files
- Read and write GMT (Gene Matrix Transposed) files
- Convert between different genomic data formats

## Installation

You can install the development version of strpip from GitHub:

```r
# install.packages("remotes")
remotes::install_github("hungms/strpip")
```

## Dependencies

The package requires the following dependencies, which can be installed with:

```r
# Install CRAN packages
install.packages(c("dplyr", "tidyr", "tibble", "stringr", "magrittr"))

# Install Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("biomaRt", "OmnipathR"))
```

These packages provide:
- biomaRt: Access to Ensembl BioMart databases
- dplyr, tidyr, tibble, stringr: Data manipulation tools
- magrittr: For pipe operations
- OmnipathR: Access to protein interaction networks and pathway data



