# Insert cRAP numbers into a character vector

This function takes a character vector where each element is roughly in
the form of a UniProt header e.g. `"sp|XXXXXX|YYYY_YYYY Text goes here"`
and substitutes a cRAP number in place of the first `|` symbol, e.g.
`"sp|cRAP001|XXXXXX|YYYY_YYYY Text goes here"`.

## Usage

``` r
sub_crap(x, start = 1, width = 3)
```

## Arguments

- x:

  `character vector`, each element must have two `|` symbols with some
  text in between e.g. `|sometext|`.

- start:

  `numeric`, the number to increment from, default is `1`

- width:

  `numeric`, how many digits the cRAP number should be, default is `3`

## Value

Returns a `character vector` the same length as x.

## Examples

``` r
# basic use
sub_crap(c("|sometext|", "|moretext|"))
#> [1] "|cRAP001|sometext|" "|cRAP002|moretext|"

# start from a different number
sub_crap(c("|sometext|", "|moretext|"), start = 88)
#> [1] "|cRAP088|sometext|" "|cRAP089|moretext|"

# increase number width
sub_crap(c("|sometext|", "|moretext|"), start = 1111, width = 4)
#> [1] "|cRAP1111|sometext|" "|cRAP1112|moretext|"
```
