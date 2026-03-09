# Convert genes between species efficiently

Convert genes between species efficiently

## Usage

``` r
convert_genes(genes, org.from = "human", org.to = "mouse", one.to.many = FALSE)
```

## Arguments

- genes:

  Vector of gene symbols to convert

- org.from:

  Source organism ("human" or "mouse")

- org.to:

  Target organism ("human" or "mouse")

- one.to.many:

  If TRUE, returns all possible mappings

## Value

Vector of converted gene symbols
