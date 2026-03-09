# strpip: Data Structure Pipelines for Genomic Analysis

A comprehensive R package for genomic data structure pipelines,
providing tools for:

- Data format conversions (GCT, GMT, TSV)

- Gene symbol mapping between species

- Gene set operations and filtering

- Gene expression data processing

- Functional annotation integration

- Specialized gene set identification (XY, MT, BCR, TCR, etc.)

## Main functions

- Data Format Operations:

  [`read_gmt()`](https://hungms.github.io/strpip/reference/read_gmt.md),
  [`read_gct()`](https://hungms.github.io/strpip/reference/read_gct.md),
  [`write_gmt()`](https://hungms.github.io/strpip/reference/write_gmt.md),
  [`write_gct()`](https://hungms.github.io/strpip/reference/write_gct.md)

- Gene Symbol Mapping:

  `convert_mouse_to_human()`, `convert_human_to_mouse()`

- Gene Set Operations:

  [`get_xy_genes()`](https://hungms.github.io/strpip/reference/get_xy_genes.md),
  [`get_mt_genes()`](https://hungms.github.io/strpip/reference/get_mt_genes.md),
  [`get_str_genes()`](https://hungms.github.io/strpip/reference/get_str_genes.md)

- Expression Processing:

  [`summarize_genes()`](https://hungms.github.io/strpip/reference/summarize_genes.md)

- Annotation:

  [`run_annotation()`](https://hungms.github.io/strpip/reference/run_annotation.md)

## Dependencies

The package depends on:

- biomaRt: For accessing Ensembl BioMart

- dplyr, tidyr, tibble, stringr: For data manipulation

- magrittr: For pipe operations

- OmniPathR: For pathway data

## See also

Useful links:

- <https://hungms.github.io/deseq2pip/>

- <https://hungms.github.io/strpip/>

- Report bugs at <https://github.com/hungms/strpip/issues>

## Author

**Maintainer**: Hung M <hungm@example.com>
