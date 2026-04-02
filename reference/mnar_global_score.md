# Compute a single dataset-level MNAR index using logistic regression

Summarises missingness structure across the whole dataset into a single
value by fitting two pooled logistic regressions across all (feature x
sample) observations simultaneously:

## Usage

``` r
mnar_global_score(
  obj,
  i,
  group_cols,
  min_observed = 2,
  min_missing = 1,
  subset_features = TRUE,
  log_transform = FALSE
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
  a feature to be included in the model. Default: 2.

- min_missing:

  Integer. Minimum number of missing values required for a feature to be
  included in the model. Default: 1.

- subset_features:

  Logical. If TRUE (default), features with fewer than `min_observed`
  observed values or fewer than `min_missing` missing values are
  excluded, as they are uninformative and would dilute the intensity and
  group effect estimates.

- log_transform:

  Logical. If TRUE, apply log2 to the assay intensities before computing
  the per-feature mean intensity covariate. Set to TRUE if your data has
  not already been log-transformed. Default: FALSE.

## Value

A list with the following elements:

- tjur_intensity_only:

  Numeric. Tjur's R² of the intensity-only model. The discrimination
  between missing and observed explained by feature abundance alone —
  the universal intensity-dependent baseline present in all proteomics
  data.

- tjur_incremental:

  Numeric. The increase in Tjur's R² when condition is added on top of
  intensity. This is the primary quantity of interest:
  condition-specific missingness beyond what intensity already explains.
  Near zero means missingness is purely intensity-driven; positive
  values indicate condition-specific depletion of particular features.

- tjur_condition_fraction:

  Numeric. The fraction of total explainable missingness (tjur_full)
  that is condition-specific rather than intensity-driven. Ranges from 0
  (all intensity-driven) to 1 (all condition-specific).

- tjur_full:

  Numeric. Tjur's R² of the full model (intensity + condition). Equal to
  tjur_intensity_only + tjur_incremental.

- lrt_chi2:

  Numeric. Chi-squared statistic from the likelihood ratio test of the
  full vs intensity-only model.

- lrt_df:

  Integer. Degrees of freedom for the LRT.

- lrt_pvalue:

  Numeric. P-value for the LRT — tests whether condition adds
  significant explanatory power over intensity alone.

- model_intensity:

  The fitted intensity-only `glm` object.

- model_full:

  The fitted full `glm` object.

- interpretation:

  Character string. Human-readable summary of all key quantities.

- n_features:

  Integer. Number of features included in the models.

- n_obs_total:

  Integer. Total number of (feature x sample) rows.

## Details

mod_intensity : missing ~ intensity mod_full : missing ~ intensity +
group

The incremental Tjur's R² (full - intensity) is the primary quantity of
interest: it measures condition-specific missingness beyond what feature
abundance alone would predict.

Tjur's R² (discrimination coefficient) is computed as:

mean(P(missing \| truly missing)) - mean(P(missing \| truly observed))

It ranges from 0 (no discrimination) to 1 (perfect discrimination) and
is a natural measure of how well the model separates missing from
observed.

Pooled logistic regression (no random effects) is used in preference to
a linear mixed model. A per-feature random intercept would absorb the
very variance that intensity is meant to explain (low-abundance features
have both low intensity and high baseline missingness), suppressing the
intensity fixed effect. The pooled approach is analogous to the
detection probability curve (DPC) fitted by the limpa package.

The intensity covariate is the per-feature mean observed intensity,
scaled to mean 0 and sd 1. If your data has not been log-transformed
upstream, set `log_transform = TRUE` so that the sigmoidal
intensity-missingness relationship is approximately linear on the log
scale, consistent with how DPC curves are typically modelled.

The decomposition is:

tjur_full = tjur_intensity_only + tjur_incremental

A large tjur_intensity_only with small tjur_incremental means
missingness is driven by abundance uniformly across conditions — the
typical baseline. A substantial tjur_incremental means certain features
are specifically depleted in particular conditions beyond what their
overall abundance would predict.

## References

Tjur T (2009). Coefficients of determination in logistic regression
models — a new proposal: the coefficient of discrimination. The American
Statistician, 63(4), 366-372.

## See also

[`mnar_score`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/mnar_score.md)

## Examples

``` r
if (FALSE) { # \dontrun{
idx <- mnar_global_score(obj, i = "peptides", group_cols = "condition")

# Primary quantities
cat("Intensity-only Tjur R²:  ", round(idx$tjur_intensity_only,  3), "\n")
cat("Incremental (condition): ", round(idx$tjur_incremental,     3), "\n")
cat("Condition fraction:      ", round(idx$tjur_condition_fraction, 3), "\n")
cat("LRT p-value:             ", format.pval(idx$lrt_pvalue),        "\n")
cat(idx$interpretation, "\n")
} # }
```
