# Remove duplicated full stops

A convenience function to remove any duplicated full stops (aka periods)
from the elements of a vector, e.g. will change `..` or `...` or `....`
etc. to just `.` Mainly used to fix column names.

## Usage

``` r
remove_dots(x)
```

## Arguments

- x:

  `character` or `string`. Contains duplicate full stops to be removed.

## Value

Returns `character` or `string` with duplicate full stops removed.

## Examples

``` r
df <- data.frame(
  column...name = c(1, 2, 3)
)

colnames(df) <- remove_dots(colnames(df))
```
