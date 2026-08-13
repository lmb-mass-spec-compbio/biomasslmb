# Summarise per-feature condition missingness scores into a dataset-level index

Aggregates the per-feature Tjur R^2 scores produced by
[`condition_miss_score()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/condition_miss_score.md)
into a single value in 0-1 representing how strongly experimental
condition predicts missingness across the dataset.

## Usage

``` r
condition_miss_index(
  summary_df,
  weight_by = c("miss_frac", "equal"),
  coverage_penalty = TRUE
)
```

## Arguments

- summary_df:

  A data.frame as returned in the `$summary` element of
  [`condition_miss_score()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/condition_miss_score.md).
  Must contain columns: `condition_miss_score`, `miss_frac`, and
  `condition_miss_class`.

- weight_by:

  Character. How to weight features when computing the mean score. One
  of:

  "miss_frac"

  :   (Default) Weight each feature by its missingness fraction.
      Features with more missing values contribute more – they represent
      a larger share of the imputation/modelling challenge.

  "equal"

  :   Unweighted mean across all informative features.

- coverage_penalty:

  Logical. If TRUE (default), scales the weighted mean score by the
  fraction of features that are informative (i.e. have at least
  `min_observed` observed and `min_missing` missing values as set in
  [`condition_miss_score()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/condition_miss_score.md)).
  When most features are fully observed and therefore uninformative, the
  dataset-level index is attenuated accordingly. If FALSE, the index
  reflects only the informative features without penalising for their
  proportion.

## Value

A list with the following elements:

- index:

  Numeric in 0-1. The dataset-level condition missingness index. Values
  near 0 indicate missingness that is not predictable from condition;
  values near 1 indicate missingness that is strongly structured by
  experimental condition.

- weighted_mean_score:

  Numeric. The weighted mean Tjur R^2 across informative features,
  before any coverage penalty is applied.

- coverage:

  Numeric. Fraction of features that were informative (had scoreable
  missingness).

- n_informative:

  Integer. Number of informative features.

- n_total:

  Integer. Total number of features.

- weight_by:

  Character. The weighting scheme used.

- coverage_penalty:

  Logical. Whether the coverage penalty was applied.

## Details

Unlike
[`global_condition_miss_score()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/global_condition_miss_score.md),
which fits a pooled logistic regression and can miss condition effects
that cancel across features (e.g. peptide A missing in condition X and
peptide B missing in condition Y), this function aggregates absolute
per-feature scores and correctly identifies such opposing patterns as
condition-structured missingness.

The index is computed as follows:

1.  Uninformative features (NA score) are excluded.

2.  The weighted mean Tjur R^2 is computed across informative features,
    using `miss_frac` or equal weights depending on `weight_by`.

3.  If `coverage_penalty = TRUE`, the result is multiplied by the
    fraction of all features that were informative.

The coverage penalty ensures that a dataset where only 5\\ have any
missingness does not score the same as one where 80\\ are missing in a
condition-structured way, even if the informative features have
identical per-feature scores.

When `coverage_penalty = FALSE` the index reflects the strength of
condition-missingness association among features that *do* have missing
values, which may be preferable when missingness prevalence is not
itself of interest.

## See also

[`condition_miss_score`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/condition_miss_score.md),
[`global_condition_miss_score`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/global_condition_miss_score.md)

## Examples

``` r
if (FALSE) { # \dontrun{
result <- condition_miss_score(obj, i = "peptides", group_cols = "condition")

# Default: miss_frac weighted, with coverage penalty
idx <- condition_miss_index(result$summary)
cat("Condition missingness index:", round(idx$index, 3), "\n")

# Strength among missing features only, no coverage penalty
idx2 <- condition_miss_index(result$summary,
                             weight_by        = "equal",
                             coverage_penalty = FALSE)
} # }
```
