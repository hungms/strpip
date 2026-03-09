# Write GMT File

Saves a data frame or list as a Gene Matrix Transposed (GMT) file. The
output file will be formatted according to the GMT specification.

## Usage

``` r
write_gmt(input, file)
```

## Arguments

- input:

  A data frame where columns are gene sets and values are genes, or a
  named list where each element is a character vector of genes

- file:

  Output file path

## Value

No return value. Creates a GMT file at the specified location.

## Examples

``` r
if (FALSE) { # \dontrun{
# Using a data frame
write_gmt(gene_sets_df, "my_genesets.gmt")

# Using a list
gene_sets_list <- list(
  pathway1 = c("gene1", "gene2", "gene3"),
  pathway2 = c("gene2", "gene4", "gene5")
)
write_gmt(gene_sets_list, "my_genesets.gmt")
} # }
```
