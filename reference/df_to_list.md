# Convert Data Frame to List

Converts a data frame into a list of vectors, removing empty strings.

## Usage

``` r
df_to_list(df)
```

## Arguments

- df:

  A data frame to convert

## Value

A list of vectors where each element corresponds to a column from the
input data frame.

## Examples

``` r
if (FALSE) { # \dontrun{
my_df <- data.frame(c("A", "B", "C"), c("D", "E"), c("F", "G", "H", "I"))
result <- df_to_list(my_df)
print(result)
} # }
```
