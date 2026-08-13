# Score features by how well missingness can be predicted by experimental condition

For each feature, fits a logistic regression of binary missingness
against one or more experimental group variables from colData. Tjur's
R^2 (discrimination coefficient) of this model is used as a per-feature
score: high values indicate that missingness is structured by
experimental condition (condition-structured), while values near zero
indicate missingness is unrelated to condition (condition-independent).

## Usage

``` r
condition_miss_score(
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

  Logical. If TRUE, scores and classifications are stored in rowData of
  the assay within the returned QFeatures object. Default: TRUE.

## Value

A list with the following elements:

- scores:

  Named numeric vector of per-feature condition missingness scores
  (Tjur's R^2).

- summary:

  A data.frame with per-feature details: n_observed, n_missing,
  miss_frac, condition_miss_score, condition_miss_class.

- global:

  A named numeric vector with dataset-level summaries:
  mean_condition_miss_score, median_condition_miss_score,
  prop_condition_structured, prop_condition_independent,
  prop_uninformative.

- obj:

  The input QFeatures object, optionally with results written into
  rowData if store_results = TRUE.

## Details

Note: in proteomics, all missingness is partly intensity-dependent –
low- abundance features are more likely to be missing regardless of
condition. This function measures the additional, condition-specific
component: whether missingness is systematically higher in some
conditions than others for a given feature. A high score does not imply
the missingness is purely non-random; a low score does not imply it is
purely random.

The score for each feature is Tjur's R^2 (discrimination coefficient)
from a logistic regression:

missing_indicator ~ group

where missing_indicator is 1 if the value is missing and 0 if observed,
and group is a factor combining the specified colData columns. Tjur's
R^2 is computed as the difference in mean predicted probabilities
between missing and observed values:

Tjur's R^2 = mean(P(missing \| truly missing)) - mean(P(missing \| truly
observed))

It ranges from 0 (no discrimination – condition-independent) to 1
(perfect discrimination – condition-structured). Logistic regression is
used in preference to a linear model as it is the correct model for a
binary outcome.

If multiple group_cols are provided they are combined into a single
interaction factor (e.g. condition:timepoint), so each unique
combination of levels is treated as a distinct group.

Classification thresholds (condition_miss_class):

- "condition_structured" : Tjur's R^2 \>= 0.5

- "mixed" : Tjur's R^2 in \[0.2, 0.5)

- "condition_independent": Tjur's R^2 \< 0.2

- "uninformative" : excluded due to min_observed / min_missing filters

These thresholds are heuristic – inspect the score distribution before
applying fixed cutoffs to your dataset.

## See also

[`global_condition_miss_score`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/global_condition_miss_score.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(QFeatures)

# Single grouping variable
result <- condition_miss_score(obj, i = "peptides", group_cols = "condition")

# Multiple variables combined as an interaction
result <- condition_miss_score(obj, i = 1,
                               group_cols = c("condition", "timepoint"))

# Inspect per-feature results
hist(result$scores, breaks = 40,
     xlab = "Condition missingness score (Tjur's R^2)",
     main = "Missingness structure by condition")
print(result$global)
head(result$summary)
} # }
```
