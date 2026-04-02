# Filter a QFeatures assay to samples with complete group annotations

Subsets a QFeatures object to the samples belonging to a specified
assay, then removes any samples that have NA in one or more of the given
colData columns. Both the assay matrix and colData are filtered
consistently.

## Usage

``` r
filter_complete_groups(obj, i, group_cols)
```

## Arguments

- obj:

  A QFeatures object.

- i:

  Integer or character. Index or name of the assay to filter.

- group_cols:

  Character vector. Column names from colData(obj) to check for NA
  values. Samples with NA in any of these columns are removed.

## Value

A list with the following elements:

- mat:

  Numeric matrix (features x samples) with NA-group samples removed.

- cd_assay:

  Data frame of colData rows corresponding to the retained samples, in
  the same column order as `mat`.

- removed:

  Character vector of sample IDs that were removed due to NA in a group
  column. Zero-length if none were removed.

## Details

This is a utility used internally by
[`mnar_score()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/mnar_score.md)
and
[`mnar_global_score()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/mnar_global_score.md),
but is exported for use in any context where you need a clean sample set
with complete group annotations — for example, before building a design
matrix or running differential abundance testing.

## See also

[`mnar_score`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/mnar_score.md),
[`mnar_global_score`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/mnar_global_score.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Remove the internal standard sample (NA condition) before modelling
filtered <- filter_complete_groups(obj, "psms_norm", c("Condition", "timepoint"))
mat      <- filtered$mat        # 16 samples, IS removed
cd       <- filtered$cd_assay   # matching colData

# Check which samples were dropped
filtered$removed
} # }
```
