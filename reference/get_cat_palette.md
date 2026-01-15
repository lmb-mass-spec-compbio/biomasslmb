# Generate a colour-blind friendly palette for categorical colour encoding

For a given number of categories, identify a suitable palette of colours
which are colour-blind friendly. Palettes are derived from
<http://mkweb.bcgsc.ca/colorblind/palettes.mhtml#page-container>

## Usage

``` r
get_cat_palette(n)
```

## Arguments

- n:

  `numeric`. The number of colours required.

## Value

Returns a `character` with the Hex codes for the colour palette.
