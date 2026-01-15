# Remove leading X

A convenience function to remove a leading capital X. Is case sensitive.

## Usage

``` r
remove_x(x)
```

## Arguments

- x:

  `character` or `string`.

## Value

Returns `character` or `string` with leading X removed.

## Examples

``` r
df <- data.frame('X1'=c(1,2))

remove_x(colnames(df))
#> [1] "1"
```
