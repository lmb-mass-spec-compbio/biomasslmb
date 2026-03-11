# Compatibility wrapper for QFeatures long-format extraction

Dispatches to the appropriate QFeatures long-format function depending
on the Bioconductor version and object class.

## Usage

``` r
qfeatures_long(object, ...)
```

## Arguments

- object:

  A QFeatures object (or object class supported by longForm/longFormat)

- ...:

  Additional arguments passed to the underlying long-format function,
  typically including `colvars` and `rowvars`.

## Value

A data structure identical to the output of the underlying long-format
function (usually a `LongTable` or `DataFrame`).

## Details

Specifically, this function:

- Uses
  [`longForm()`](https://rdrr.io/pkg/BiocGenerics/man/longForm.html)
  (BiocGenerics generic) if a method exists for the input object,

- Falls back to `QFeatures::longFormat()` otherwise.

This allows package or script code to remain compatible across QFeatures
releases without explicitly checking Bioconductor versions.

The function first checks whether
[`longForm()`](https://rdrr.io/pkg/BiocGenerics/man/longForm.html) is a
generic and whether there is a registered method for the input object
class. If so, it calls
[`longForm()`](https://rdrr.io/pkg/BiocGenerics/man/longForm.html).
Otherwise, it falls back to the legacy `longFormat()` function. An error
is thrown if neither is available.

## See also

`longFormat`,
[`longForm`](https://rdrr.io/pkg/BiocGenerics/man/longForm.html),
[`isGeneric`](https://rdrr.io/r/methods/GenericFunctions.html),
[`findMethods`](https://rdrr.io/r/methods/findMethods.html)
