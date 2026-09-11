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

## Examples

``` r
# a colourblind-friendly categorical palette
get_cat_palette(4)
#> [1] "#2271B2" "#d55e00" "#359B73" "#e69f00"

library(ggplot2)
ggplot(mtcars, aes(wt, mpg, colour = factor(cyl))) +
  geom_point() +
  scale_colour_manual(values = get_cat_palette(3)) +
  theme_biomasslmb()
```
