# Get cell cycle genes

Retrieves gene symbols from cell cycle gene sets for a specified
organism.

## Usage

``` r
get_cc_genes(org)
```

## Arguments

- org:

  Organism to query. Either "human" or "mouse".

- ...:

  Additional arguments passed to `import_biomart_human` or
  `import_biomart_mouse`.

## Value

A character vector containing cell cycle gene symbols.

## Examples

``` r
if (FALSE) { # \dontrun{
cc_genes <- get_cc_genes("human")
head(cc_genes)
} # }
```
