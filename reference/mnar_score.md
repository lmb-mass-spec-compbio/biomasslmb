# Quantify the degree of MNAR vs MAR missingness in a QFeatures assay

For each feature, fits a logistic regression of binary missingness
against one or more experimental group variables from colData. Tjur's R²
(discrimination coefficient) of this model is used as a per-feature MNAR
score: high values indicate that missingness is structured by
experimental condition (MNAR-like), while values near zero indicate
missingness is unrelated to condition (MAR-like).

## Usage

``` r
mnar_score(
  obj,
  i,
  group_cols,
  min_observed = 2,
  min_missing = 1,
  store_results = TRUE
)
```

## Arguments

- obj:

  A QFeatures object.

- i:

  Integer or character. Index or name of the assay to analyse.

- group_cols:

  Character vector. One or more column names from colData(obj) that
  define experimental groups (e.g. c("condition", "treatment")).
  Multiple columns are combined into an interaction term.

- min_observed:

  Integer. Minimum number of observed (non-missing) values required for
  a feature to be included. Features below this threshold are returned
  as NA. Default: 2.

- min_missing:

  Integer. Minimum number of missing values required for a feature to be
  included. Features with no missingness are uninformative and returned
  as NA. Default: 1.

- store_results:

  Logical. If TRUE, MNAR scores and classifications are stored in
  rowData of the assay within the returned QFeatures object. Default:
  TRUE.

## Value

A list with the following elements:

- scores:

  Named numeric vector of per-feature MNAR scores (Tjur's R²).

- summary:

  A data.frame with per-feature details: n_observed, n_missing,
  miss_frac, mnar_score, mnar_class.

- global:

  A named numeric vector with dataset-level summaries: mean_mnar_score,
  median_mnar_score, prop_mnar, prop_mar, prop_uninformative.

- obj:

  The input QFeatures object, optionally with results written into
  rowData if store_results = TRUE.

## Details

The MNAR score for each feature is Tjur's R² (discrimination
coefficient) from a logistic regression:

missing_indicator ~ group

where missing_indicator is 1 if the value is missing and 0 if observed,
and group is a factor combining the specified colData columns. Tjur's R²
is computed as the difference in mean predicted probabilities between
missing and observed values:

Tjur's R² = mean(P(missing \| truly missing)) - mean(P(missing \| truly
observed))

It ranges from 0 (no discrimination — MAR) to 1 (perfect discrimination
— MNAR). Logistic regression is used in preference to a linear model as
it is the correct model for a binary outcome.

If multiple group_cols are provided they are combined into a single
interaction factor (e.g. condition:timepoint), so each unique
combination of levels is treated as a distinct group.

MNAR classification thresholds (mnar_class):

- "MNAR" : Tjur's R² \>= 0.5

- "mixed" : Tjur's R² in \[0.2, 0.5)

- "MAR" : Tjur's R² \< 0.2

- "uninformative" : excluded due to min_observed / min_missing filters

These thresholds are heuristic — inspect the score distribution before
applying fixed cutoffs to your dataset.

## See also

[`mnar_global_score`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/mnar_global_score.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(QFeatures)

# Single grouping variable
result <- mnar_score(obj, i = "peptides", group_cols = "condition")

# Multiple variables combined as an interaction
result <- mnar_score(obj, i = 1,
                     group_cols = c("condition", "timepoint"))

# Inspect per-feature results
hist(result$scores, breaks = 40,
     xlab = "MNAR score (Tjur's R²)", main = "Missingness structure")
print(result$global)
head(result$summary)
} # }
```
