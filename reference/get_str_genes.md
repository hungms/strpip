# Get genes matching specific patterns

Retrieves gene symbols matching predefined patterns for specific gene
families (BCR, TCR, MHC, HB, RB, MT) for a specified organism.

## Usage

``` r
get_str_genes(org, str, ...)
```

## Arguments

- org:

  Organism to query. Either "human" or "mouse".

- str:

  Character vector of patterns to search for. Valid options are: "bcr",
  "tcr", "mhc", "hb", "rb", "mt".

## Value

A character vector containing matching gene symbols.

## Examples

``` r
if (FALSE) { # \dontrun{
rb_genes <- get_str_genes("human", "rb")
head(rb_genes)
} # }
```
