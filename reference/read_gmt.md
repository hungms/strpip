# Read GMT File

Reads a Gene Matrix Transposed (GMT) file and converts it to a data
frame. GMT files are commonly used for gene set collections.

## Usage

``` r
read_gmt(gmt)
```

## Arguments

- gmt:

  Path to the GMT file

## Value

A data frame where: - Rows are genes - Columns are gene sets - Values
indicate gene set membership

## Examples

``` r
if (FALSE) { # \dontrun{
gmt_data <- read_gmt("path/to/genesets.gmt")
head(gmt_data)
} # }
```
