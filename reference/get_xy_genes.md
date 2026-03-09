# Get genes from X and Y chromosomes

Retrieves all gene symbols located on the X and Y chromosomes for a
specified organism.

## Usage

``` r
get_xy_genes(org, ...)
```

## Arguments

- org:

  Organism to query. Either "human" or "mouse".

- ...:

  Additional arguments passed to `import_biomart_human` or
  `import_biomart_mouse`.

## Value

A character vector containing gene symbols from X and Y chromosomes.

## Examples

``` r
if (FALSE) { # \dontrun{
xy_genes <- get_xy_genes("human")
head(xy_genes)
} # }
```
