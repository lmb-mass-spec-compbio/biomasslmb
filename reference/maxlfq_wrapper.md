# Wrapper for iq::maxLFQ to use with QFeatures::aggregateFeatures

This function wraps `iq::maxLFQ()` so it can be passed to
[`QFeatures::aggregateFeatures()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html)
as the `fun` argument. It expects a numeric matrix with peptides as rows
and samples as columns, and returns a numeric vector of protein-level
intensities.

## Usage

``` r
maxlfq_wrapper(mat, ...)
```

## Arguments

- mat:

  A numeric matrix of peptide intensities (rows = peptides, columns =
  samples). Missing values should be `NA`.

- ...:

  Additional arguments passed to `iq::maxLFQ()`.

## Value

A numeric vector of length equal to the number of columns in `mat`,
containing protein-level abundance estimates.

## Details

The MaxLFQ algorithm computes protein abundances from peptide
intensities by using peptide ratios across samples and solving a
least-squares system. This wrapper extracts the `estimate` component
returned by `iq::maxLFQ()` and coerces it to a numeric vector, suitable
for
[`aggregateFeatures()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html).

## Examples

``` r
if (FALSE) { # \dontrun{
library(QFeatures)
library(iq)

# Example peptide-level matrix (log2 intensities)
mat <- matrix(rnorm(20), nrow = 5, ncol = 4)
colnames(mat) <- paste0("Sample", 1:4)
rownames(mat) <- paste0("Pep", 1:5)

protein_estimates <- maxlfq_wrapper(mat)
} # }
```
