# Read GCT File

Reads a Gene Cluster Text (GCT) file and converts it to a data frame.
GCT files are commonly used for gene expression data.

## Usage

``` r
read_gct(gct)
```

## Arguments

- gct:

  Path to the GCT file

## Value

A data frame where: - Rows are genes - Columns are samples - Values are
expression measurements

## Examples

``` r
if (FALSE) { # \dontrun{
gct_data <- read_gct("path/to/expression.gct")
head(gct_data)
} # }
```
