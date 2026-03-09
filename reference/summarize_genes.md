# Summarize Gene Expression

This function aggregates expression values from multiple isoforms of the
same gene into a single gene-level expression value.

## Usage

``` r
summarize_genes(input, gene_sym_vec, normalized = FALSE)
```

## Arguments

- input:

  Gene expression data with isoforms/peaks in rows and samples in
  columns. Can be a data frame, data.table, or matrix.

- gene_sym_vec:

  Vector of gene symbols corresponding to each isoform/peak

- normalized:

  Logical. If TRUE, calculates the mean expression across
  isoforms/peaks. If FALSE, calculates the sum of expression across
  isoforms/peaks. Default is FALSE

## Value

A data frame, data.table, or matrix with unique gene names as row names.
Output type matches the class of the input.
