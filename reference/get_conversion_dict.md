# Convert mouse gene symbols

Converts a vector of mouse (MGI) gene symbols to their human (HGNC)
equivalents. By default, returns a data frame with all possible
mappings.

## Usage

``` r
get_conversion_dict(org.from, org.to)
```

## Arguments

- genes:

  A vector of mouse gene symbols to convert

## Value

A data frame with all possible mappings.

## Examples

``` r
if (FALSE) { # \dontrun{
mouse_genes <- c("Trp53", "Cd4", "Cd8a")
human_genes <- convert_mouse_to_human(mouse_genes)
print(human_genes)
} # }
```
