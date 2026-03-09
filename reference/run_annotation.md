# Add Gene Annotations

Adds functional annotations to a data frame containing gene symbols
using pre-computed OmniPath annotations. The function supports both
human and mouse gene symbols.

## Usage

``` r
run_annotation(
  df,
  gene_column = "gene",
  org.from = "human",
  org.to = "mouse",
  release = "105",
  one.to.many = FALSE
)
```

## Arguments

- df:

  Data frame containing gene symbols

- gene_column:

  Name of the column containing gene symbols. Default is "gene"

- org.from:

  Organism type. Either "human" or "mouse". Default is "human"

- org.to:

  Target organism. Either "human" or "mouse". Default is "mouse"

- release:

  Ensembl release version to use for annotations. Default is "105"

- one.to.many:

  Logical. Default is FALSE to return only unique mappings (one/many to
  one). If TRUE, returns all possible mouse gene mappings (one to many),
  including cases where one human gene maps to multiple mouse genes.

## Value

A data frame with additional columns containing functional
annotations: - Pathway information - Protein-protein interactions - Gene
regulatory relationships - Other functional annotations from OmniPath
