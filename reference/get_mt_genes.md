# Get genes from mitochondrial chromosome

Retrieves all gene symbols located on the mitochondrial chromosome for a
specified organism.

## Usage

``` r
get_mt_genes(org, ...)
```

## Arguments

- org:

  Organism to query. Either "human" or "mouse".

- ...:

  Additional arguments passed to `import_biomart_human` or
  `import_biomart_mouse`.

## Value

A character vector containing mitochondrial gene symbols.

## Examples

``` r
if (FALSE) { # \dontrun{
mt_genes <- get_mt_genes("human")
head(mt_genes)
} # }
```
