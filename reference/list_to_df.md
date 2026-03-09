# Convert List to Data Frame

Converts a list of vectors into a data frame, padding shorter vectors
with empty strings to match the length of the longest vector.

## Usage

``` r
list_to_df(list)
```

## Arguments

- list:

  A list of vectors to convert

## Value

A data frame where each column corresponds to a vector from the input
list. Shorter vectors are padded with empty strings to match the length
of the longest vector.

## Examples

``` r
if (FALSE) { # \dontrun{
my_list <- list(c("A", "B", "C"), c("D", "E"), c("F", "G", "H", "I"))
result <- list_to_df(my_list)
print(result)
} # }
```
