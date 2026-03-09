# Convert expression matrix or data frame between species efficiently

Convert expression matrix or data frame between species efficiently

## Usage

``` r
convert_exprs(
  exprs,
  org.from = "human",
  org.to = "mouse",
  many.to.one = TRUE,
  normalized = FALSE
)
```

## Arguments

- exprs:

  Expression matrix, data frame, or data.table

- org.from:

  Source organism ("human" or "mouse")

- org.to:

  Target organism ("human" or "mouse")

- many.to.one:

  If TRUE, aggregates multiple mappings

- normalized:

  If TRUE, uses mean; if FALSE, uses sum

## Value

Converted expression data with same type as input
