# Make Elements of a Character Vector Unique (Numbering From 1)

A variant of
[`base::make.unique()`](https://rdrr.io/r/base/make.unique.html) that
appends numeric suffixes to elements of a character vector. Unlike
[`make.unique()`](https://rdrr.io/r/base/make.unique.html), which leaves
the first occurrence unchanged, this function can number all members of
a duplicated group starting from 1. Optionally, it can number singletons
as well.

## Usage

``` r
make_unique_all(x, sep = ".", always = TRUE)
```

## Arguments

- x:

  A character vector.

- sep:

  A character string to separate the original element from the suffix.
  Defaults to `"."`.

- always:

  Logical. If `TRUE`, all elements (including singletons) get a suffix
  (e.g. `"B"` becomes `"B.1"`). If `FALSE` (default), singletons remain
  unchanged.

## Value

A character vector of the same length as `x`, with duplicates
disambiguated by appending a numeric suffix.

## See also

[`base::make.unique()`](https://rdrr.io/r/base/make.unique.html)

## Examples

``` r
# Mimics make.unique() but numbers first duplicates as well
make_unique_all(c("A", "A", "B"))
#> [1] "A.1" "A.2" "B.1"
# [1] "A.1" "A.2" "B"

# Multiple duplicate groups
make_unique_all(c("A", "A", "B", "B", "C"))
#> [1] "A.1" "A.2" "B.1" "B.2" "C.1"
# [1] "A.1" "A.2" "B.1" "B.2" "C"

# Force suffixes on all elements
make_unique_all(c("A", "A", "B"), always = TRUE)
#> [1] "A.1" "A.2" "B.1"
# [1] "A.1" "A.2" "B.1"

# Custom separator
make_unique_all(c("X", "X", "Y"), sep = "_")
#> [1] "X_1" "X_2" "Y_1"
# [1] "X_1" "X_2" "Y"
```
