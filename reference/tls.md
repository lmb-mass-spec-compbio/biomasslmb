# Total least squares regression (for use with geom_smooth)

Total least squares regression (for use with geom_smooth)

## Usage

``` r
tls(formula, data, ...)
```

## Arguments

- formula:

  Formula of the form y ~ x.

- data:

  A data frame containing columns for x and y.

- ...:

  Ignored.

## Examples

``` r
# both axes are noisy estimates, so ordinary least squares would give a
# slope that depends on which variable is on the x axis
x <- rnorm(100)
df <- data.frame(x = x + rnorm(100, sd = 0.3),
                 y = x + rnorm(100, sd = 0.3))

unlist(tls(y ~ x, df))
#>   intercept       slope 
#> -0.02784179  0.95054451 

library(ggplot2)
ggplot(df, aes(x, y)) +
  geom_point() +
  geom_smooth(method = tls, formula = y ~ x, se = FALSE) +
  theme_biomasslmb()
```
