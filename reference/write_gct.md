# Write GCT File

Saves a data frame as a Gene Cluster Text (GCT) file. The output file
will be formatted according to the GCT specification.

## Usage

``` r
write_gct(df, file)
```

## Arguments

- df:

  Data frame to save

- save_name:

  Name of the output file

- save_dir:

  Directory to save the file in. Defaults to current working directory.

## Value

No return value. Creates a GCT file at the specified location.

## Examples

``` r
if (FALSE) { # \dontrun{
write_gct(expression_df, "my_expression.gct", "output/")
} # }
```
