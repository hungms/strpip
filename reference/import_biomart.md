# Import BioMart data from Ensembl archive

Retrieves archived ensembl biomarts from web and creates a mapping
dictionary between human and mouse genes.

## Usage

``` r
import_biomart(
  host = "https://dec2021.archive.ensembl.org",
  local = TRUE,
  release = "105"
)
```

## Arguments

- host:

  URL to retrieve archived ensembl biomarts. Default is the December
  2021 archive.

## Value

A data frame containing gene mappings and chromosome information

## Examples

``` r
if (FALSE) { # \dontrun{
biodict <- import_biomart()
head(biodict)
} # }
```
