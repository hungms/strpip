# Convert data frame between species efficiently

Convert data frame between species efficiently

## Usage

``` r
convert_df(
  df,
  gene_column,
  org.from = "human",
  org.to = "mouse",
  one.to.many = FALSE
)
```

## Arguments

- df:

  Data frame containing gene symbols

- gene_column:

  Column name containing gene symbols

- org.from:

  Source organism ("human" or "mouse")

- org.to:

  Target organism ("human" or "mouse")

- one.to.many:

  If TRUE, returns all possible mappings

## Value

Data frame with converted gene symbols
