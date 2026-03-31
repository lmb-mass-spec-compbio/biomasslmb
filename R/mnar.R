#' Filter a QFeatures assay to samples with complete group annotations
#'
#' Subsets a QFeatures object to the samples belonging to a specified assay,
#' then removes any samples that have NA in one or more of the given colData
#' columns. Both the assay matrix and colData are filtered consistently.
#'
#' This is a utility used internally by \code{mnar_score()} and
#' \code{mnar_global_score()}, but is exported for use in any context where
#' you need a clean sample set with complete group annotations — for example,
#' before building a design matrix or running differential abundance testing.
#'
#' @param obj A QFeatures object.
#' @param i Integer or character. Index or name of the assay to filter.
#' @param group_cols Character vector. Column names from colData(obj) to check
#'   for NA values. Samples with NA in any of these columns are removed.
#'
#' @return A list with the following elements:
#'   \describe{
#'     \item{mat}{Numeric matrix (features x samples) with NA-group samples
#'       removed.}
#'     \item{cd_assay}{Data frame of colData rows corresponding to the
#'       retained samples, in the same column order as \code{mat}.}
#'     \item{removed}{Character vector of sample IDs that were removed due
#'       to NA in a group column. Zero-length if none were removed.}
#'   }
#'
#' @examples
#' \dontrun{
#' # Remove the internal standard sample (NA condition) before modelling
#' filtered <- filter_complete_groups(obj, "psms_norm", c("Condition", "timepoint"))
#' mat      <- filtered$mat        # 16 samples, IS removed
#' cd       <- filtered$cd_assay   # matching colData
#'
#' # Check which samples were dropped
#' filtered$removed
#' }
#'
#' @seealso \code{\link{mnar_score}}, \code{\link{mnar_global_score}}
#' @importFrom SummarizedExperiment assay colData rowData
#' @export
filter_complete_groups <- function(obj, i, group_cols) {

  # ── Input validation ──────────────────────────────────────────────────────────

  check_q(obj)
  check_se_exists(obj, i)

  cd           <- as.data.frame(SummarizedExperiment::colData(obj))
  missing_cols <- setdiff(group_cols, colnames(cd))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "Column(s) not found in colData: %s\nAvailable columns: %s",
      paste(missing_cols, collapse = ", "),
      paste(colnames(cd), collapse = ", ")
    ))
  }

  # ── Subset colData to this assay's samples ────────────────────────────────────

  mat        <- SummarizedExperiment::assay(obj[[i]])
  sample_ids <- colnames(mat)

  if (!all(sample_ids %in% rownames(cd))) {
    stop(
      "Some sample IDs in the assay are not present in colData. ",
      "This should not happen with a well-formed QFeatures object."
    )
  }

  cd_assay <- cd[sample_ids, , drop = FALSE]

  # ── Identify and remove samples with NA in any group column ───────────────────

  na_mask <- rowSums(is.na(cd_assay[, group_cols, drop = FALSE])) > 0
  removed <- rownames(cd_assay)[na_mask]

  if (length(removed) > 0) {
    message(sprintf(
      "Dropping %d sample(s) with NA in group column(s) [%s]: %s",
      length(removed),
      paste(group_cols, collapse = ", "),
      paste(removed, collapse = ", ")
    ))
    cd_assay <- cd_assay[!na_mask, , drop = FALSE]
    mat      <- mat[, rownames(cd_assay), drop = FALSE]
  }

  list(
    mat      = mat,
    cd_assay = cd_assay,
    removed  = removed
  )
}


#' Quantify the degree of MNAR vs MAR missingness in a QFeatures assay
#'
#' For each protein/feature, fits a linear model of binary missingness
#' against one or more experimental group variables from colData. The R²
#' of this model is used as a per-feature MNAR score: high R² indicates
#' that missingness is structured by experimental condition (MNAR-like),
#' while R² near zero indicates missingness is unrelated to condition (MAR-like).
#'
#' @param obj A QFeatures object.
#' @param i Integer or character. Index or name of the assay to analyse.
#' @param group_cols Character vector. One or more column names from colData(obj)
#'   that define experimental groups (e.g. c("condition", "treatment")).
#'   Multiple columns are combined into an interaction term.
#' @param min_observed Integer. Minimum number of observed (non-missing) values
#'   required for a feature to be included. Features below this threshold are
#'   returned as NA. Default: 2.
#' @param min_missing Integer. Minimum number of missing values required for a
#'   feature to be included. Features with no missingness are uninformative and
#'   returned as NA. Default: 1.
#' @param store_results Logical. If TRUE, MNAR scores and classifications are
#'   stored in rowData of the assay within the returned QFeatures object.
#'   Default: TRUE.
#'
#' @return A list with the following elements:
#'   \describe{
#'     \item{scores}{Named numeric vector of per-feature MNAR scores (R²).}
#'     \item{summary}{A data.frame with per-feature details: n_observed,
#'       n_missing, miss_frac, mnar_score, mnar_class.}
#'     \item{global}{A named numeric vector with dataset-level summaries:
#'       mean_mnar_score, median_mnar_score, prop_mnar, prop_mar,
#'       prop_uninformative.}
#'     \item{obj}{The input QFeatures object, optionally with results written
#'       into rowData if store_results = TRUE.}
#'   }
#'
#' @details
#' The MNAR score for each feature is the R² from a linear model:
#'
#'   missing_indicator ~ group
#'
#' where missing_indicator is 1 if the value is missing and 0 if observed,
#' and group is a factor combining the specified colData columns.
#'
#' If multiple group_cols are provided they are combined into a single
#' interaction factor (e.g. condition:timepoint), so each unique combination
#' of levels is treated as a distinct group.
#'
#' MNAR classification thresholds (mnar_class):
#'   - "MNAR"          : R² >= 0.5
#'   - "mixed"         : R² in [0.2, 0.5)
#'   - "MAR"           : R² <  0.2
#'   - "uninformative" : excluded due to min_observed / min_missing filters
#'
#' These thresholds are heuristic — inspect the score distribution before
#' applying fixed cutoffs to your dataset.
#'
#' To summarise the dataset to a single value, pass the output of this
#' function to \code{mnar_global_score()}.
#'
#' @examples
#' \dontrun{
#' library(QFeatures)
#'
#' # Single grouping variable
#' result <- mnar_score(obj, i = "peptides", group_cols = "condition")
#'
#' # Multiple variables combined as an interaction
#' result <- mnar_score(obj, i = 1,
#'                      group_cols = c("condition", "timepoint"))
#'
#' # Inspect per-feature results
#' hist(result$scores, breaks = 40,
#'      xlab = "MNAR score (R²)", main = "Missingness structure")
#' print(result$global)
#' head(result$summary)
#'
#' # Summarise to a single dataset-level value via mixed model
#' idx <- mnar_global_score(result,
#'                          i          = "peptides",
#'                          group_cols = "condition")
#' cat("Global MNAR index:", round(idx$r2_incremental, 3), "\n")
#' cat(idx$interpretation, "\n")
#' }
#'
#' @seealso \code{\link{mnar_global_score}}
#' @export
mnar_score <- function(obj,
                       i,
                       group_cols,
                       min_observed  = 2,
                       min_missing   = 1,
                       store_results = TRUE) {

  # ── Input validation ──────────────────────────────────────────────────────────

  check_q(obj)
  check_se_exists(obj, i)

  cd           <- as.data.frame(SummarizedExperiment::colData(obj))
  missing_cols <- setdiff(group_cols, colnames(cd))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "Column(s) not found in colData: %s\nAvailable columns: %s",
      paste(missing_cols, collapse = ", "),
      paste(colnames(cd), collapse = ", ")
    ))
  }

  if (!is.numeric(min_observed) || min_observed < 1) {
    stop("`min_observed` must be a positive integer.")
  }
  if (!is.numeric(min_missing) || min_missing < 1) {
    stop("`min_missing` must be a positive integer.")
  }

  # ── Extract assay and construct group factor ──────────────────────────────────

  se         <- obj[[i]]
  mat        <- SummarizedExperiment::assay(se)
  n_features <- nrow(mat)
  n_samples  <- ncol(mat)

  assay_name <- names(obj)[if (is.numeric(i)) i else which(names(obj) == i)]
  message(sprintf("Analysing assay '%s': %d features x %d samples",
                  assay_name, n_features, n_samples))

  # Filter to assay samples with complete group annotations
  filtered  <- filter_complete_groups(obj, i, group_cols)
  mat       <- filtered$mat
  cd_assay  <- filtered$cd_assay
  n_samples <- ncol(mat)

  if (length(group_cols) == 1) {
    group <- factor(cd_assay[[group_cols]])
    message(sprintf("Group variable: %s (%d levels: %s)",
                    group_cols, nlevels(group),
                    paste(levels(group), collapse = ", ")))
  } else {
    group_df <- cd_assay[, group_cols, drop = FALSE]
    group    <- factor(do.call(interaction,
                               c(group_df, list(sep = ":", drop = TRUE))))
    message(sprintf("Group interaction: %s (%d unique combinations)",
                    paste(group_cols, collapse = " x "), nlevels(group)))
  }

  if (nlevels(group) < 2) {
    stop("At least 2 group levels are required to fit the missingness model. ",
         "Check that your group_cols have sufficient variation.")
  }

  # ── Per-feature MNAR scoring ──────────────────────────────────────────────────

  miss_mat <- is.na(mat)

  scores <- vapply(seq_len(n_features), function(i) {
    m      <- miss_mat[i, ]
    n_miss <- sum(m)
    n_obs  <- sum(!m)
    if (n_obs < min_observed || n_miss < min_missing) return(NA_real_)

    # Perfect fit short-circuit: if missingness is fully determined by group
    # (e.g. always missing in condition A, never in condition B), the lm will
    # produce a perfect fit warning. R² = 1 is correct in this case — it is
    # the strongest possible MNAR signal — so we return it directly.
    per_group_var <- tapply(as.numeric(m), group, stats::var)
    if (all(per_group_var == 0, na.rm = TRUE)) return(1.0)

    df  <- data.frame(missing = as.numeric(m), group = group)
    mod <- stats::lm(missing ~ group, data = df)
    summary(mod)$r.squared
  }, numeric(1))

  names(scores) <- rownames(mat)

  # ── Per-feature summary table ─────────────────────────────────────────────────

  n_obs_vec  <- rowSums(!miss_mat)
  n_miss_vec <- rowSums(miss_mat)

  mnar_class <- dplyr::case_when(
    is.na(scores)  ~ "uninformative",
    scores >= 0.5  ~ "MNAR",
    scores >= 0.2  ~ "mixed",
    TRUE           ~ "MAR"
  )

  summary_df <- data.frame(
    feature    = rownames(mat),
    n_observed = n_obs_vec,
    n_missing  = n_miss_vec,
    miss_frac  = round(n_miss_vec / n_samples, 3),
    mnar_score = round(scores, 4),
    mnar_class = mnar_class,
    row.names  = rownames(mat),
    stringsAsFactors = FALSE
  )

  # ── Dataset-level summary ─────────────────────────────────────────────────────

  n_total  <- nrow(summary_df)
  n_inform <- sum(mnar_class != "uninformative")

  global <- c(
    n_features_total         = n_total,
    n_features_informative   = n_inform,
    n_features_uninformative = n_total - n_inform,
    mean_mnar_score          = round(mean(scores, na.rm = TRUE), 4),
    median_mnar_score        = round(stats::median(scores, na.rm = TRUE), 4),
    prop_mnar                = round(mean(mnar_class == "MNAR"), 4),
    prop_mixed               = round(mean(mnar_class == "mixed"), 4),
    prop_mar                 = round(mean(mnar_class == "MAR"), 4),
    prop_uninformative       = round(mean(mnar_class == "uninformative"), 4)
  )

  message(sprintf(
    "Results: %d informative features | mean MNAR score = %.3f | MNAR: %.1f%% | MAR: %.1f%%",
    n_inform, global["mean_mnar_score"],
    global["prop_mnar"] * 100, global["prop_mar"] * 100
  ))

  # ── Optionally write results into rowData ─────────────────────────────────────

  if (store_results) {
    rd            <- as.data.frame(SummarizedExperiment::rowData(se))
    rd$mnar_score <- scores[rownames(rd)]
    rd$mnar_class <- mnar_class[rownames(rd)]
    SummarizedExperiment::rowData(obj[[i]]) <- rd
  }

  list(
    scores  = scores,
    summary = summary_df,
    global  = global,
    obj     = obj
  )
}


#' Compute a single dataset-level MNAR index using a linear mixed model
#'
#' Summarises missingness structure across the whole dataset into a single
#' value by fitting one mixed model across all features simultaneously:
#'
#'   missing ~ group + (1 | feature)
#'
#' Two mixed models are fitted and compared to decompose missingness into
#' intensity-driven and condition-specific components:
#'
#'   mod_intensity : missing ~ intensity + (1 | feature)
#'   mod_full      : missing ~ intensity + group + (1 | feature)
#'
#' The incremental marginal R² (full - intensity) is the primary quantity of
#' interest: it measures condition-specific missingness beyond what feature
#' abundance alone would predict.
#'
#' @param mnar_result A list returned by \code{mnar_score()}.
#' @param i Integer or character. Must match the \code{i} used
#'   in the \code{mnar_score()} call.
#' @param group_cols Character vector. Must match the \code{group_cols} used
#'   in the \code{mnar_score()} call.
#' @param subset_features Logical. If TRUE (default), only features with a
#'   non-NA per-feature score are included in the mixed model. Features
#'   excluded by the min_observed / min_missing filters in \code{mnar_score()}
#'   are uninformative and would dilute the fixed effect estimate.
#'
#' @return A list with the following elements:
#'   \describe{
#'     \item{r2_intensity_only}{Numeric. Marginal R² of the intensity-only
#'       model. The proportion of missingness variance explained by feature
#'       abundance alone — the universal intensity-dependent baseline present
#'       in all proteomics data.}
#'     \item{r2_incremental}{Numeric. The increase in marginal R² when
#'       condition is added on top of intensity. This is the primary quantity
#'       of interest: condition-specific missingness beyond what intensity
#'       already explains. Near zero means missingness is purely
#'       intensity-driven; positive values indicate condition-specific
#'       depletion of particular features.}
#'     \item{r2_condition_fraction}{Numeric. The fraction of total explainable
#'       missingness (r2_full_marginal) that is condition-specific rather than
#'       intensity-driven. Ranges from 0 (all intensity-driven) to 1 (all
#'       condition-specific).}
#'     \item{r2_full_marginal}{Numeric. Marginal R² of the full model
#'       (intensity + condition). Equal to r2_intensity_only + r2_incremental.}
#'     \item{r2_full_conditional}{Numeric. Conditional R² of the full model.
#'       NA if the model is singular.}
#'     \item{lrt_chi2}{Numeric. Chi-squared statistic from the likelihood
#'       ratio test of the full vs intensity-only model.}
#'     \item{lrt_df}{Integer. Degrees of freedom for the LRT.}
#'     \item{lrt_pvalue}{Numeric. P-value for the LRT — tests whether
#'       condition adds significant explanatory power over intensity alone.}
#'     \item{singular_intensity}{Logical. Whether the intensity-only model
#'       is singular.}
#'     \item{singular_full}{Logical. Whether the full model is singular.}
#'     \item{model_intensity}{The fitted intensity-only \code{lmerMod} object.}
#'     \item{model_full}{The fitted full \code{lmerMod} object.}
#'     \item{interpretation}{Character string. Human-readable summary of all
#'       key quantities.}
#'     \item{n_features}{Integer. Number of features included in the models.}
#'     \item{n_obs_total}{Integer. Total number of (feature x sample) rows.}
#'   }
#'
#' @details
#' The intensity covariate is the per-feature mean observed intensity, scaled
#' to mean 0 and sd 1. This captures each feature's typical abundance level,
#' which drives intensity-dependent missingness. It is a feature-level
#' (not sample-level) covariate, so the same value is used for all samples
#' of a given feature.
#'
#' The decomposition is:
#'
#'   Total explainable missingness = intensity effect + condition effect
#'   r2_full_marginal              = r2_intensity_only + r2_incremental
#'
#' A large r2_intensity_only with small r2_incremental means missingness is
#' driven by abundance uniformly across conditions — the typical baseline.
#' A substantial r2_incremental means certain features are specifically
#' depleted in particular conditions beyond what their overall abundance
#' would predict.
#'
#' Marginal R² values use the Nakagawa & Schielzeth (2013) method via
#' \code{performance::r2()}. Both models use ML (REML = FALSE) to allow
#' likelihood ratio testing and for valid marginal R² computation.
#'
#' Requires packages: lme4, performance.
#'
#' @references
#' Nakagawa S & Schielzeth H (2013). A general and simple method for
#' obtaining R² from generalised linear mixed-effects models.
#' Methods in Ecology and Evolution, 4(2), 133-142.
#'
#' @examples
#' \dontrun{
#' result <- mnar_score(obj, i = "peptides", group_cols = "condition")
#'
#' idx <- mnar_global_score(result,
#'                          i          = "peptides",
#'                          group_cols = "condition")
#'
#' # Primary quantities
#' cat("Intensity-only R²:       ", round(idx$r2_intensity_only,     3), "\n")
#' cat("Incremental (condition): ", round(idx$r2_incremental,        3), "\n")
#' cat("Condition fraction:      ", round(idx$r2_condition_fraction, 3), "\n")
#' cat("LRT p-value:             ", format.pval(idx$lrt_pvalue),         "\n")
#' cat(idx$interpretation, "\n")
#' }
#'
#' @seealso \code{\link{mnar_score}}
#' @export
mnar_global_score <- function(mnar_result,
                              i,
                              group_cols,
                              subset_features = TRUE) {

  # ── Check dependencies ────────────────────────────────────────────────────────

  if (!requireNamespace("lme4", quietly = TRUE)) {
    stop("Package 'lme4' is required. Install with: install.packages('lme4')")
  }
  if (!requireNamespace("performance", quietly = TRUE)) {
    stop("Package 'performance' is required. Install with: install.packages('performance')")
  }

  # ── Input validation ──────────────────────────────────────────────────────────

  if (!is.list(mnar_result) ||
      !all(c("scores", "summary", "obj") %in% names(mnar_result))) {
    stop(
      "`mnar_result` must be the list returned by mnar_score(), ",
      "containing $scores, $summary, and $obj."
    )
  }

  obj    <- mnar_result$obj
  scores <- mnar_result$scores

  # ── Reconstruct assay matrix and group factor ─────────────────────────────────
  # Re-derived from the QFeatures object to ensure consistency with the
  # original mnar_score() call.

  mat <- SummarizedExperiment::assay(obj[[i]])

  # Filter to assay samples with complete group annotations
  filtered <- filter_complete_groups(obj, i, group_cols)
  mat      <- filtered$mat
  cd_assay <- filtered$cd_assay

  if (length(group_cols) == 1) {
    group <- factor(cd_assay[[group_cols]])
  } else {
    group_df <- cd_assay[, group_cols, drop = FALSE]
    group    <- factor(do.call(interaction,
                               c(group_df, list(sep = ":", drop = TRUE))))
  }

  if (nlevels(group) < 2) {
    stop("At least 2 group levels are required.")
  }

  # ── Optionally subset to informative features ─────────────────────────────────

  if (subset_features) {
    keep <- !is.na(scores)
    if (sum(keep) == 0) {
      stop("No informative features found (all scores are NA). ",
           "Check min_observed / min_missing thresholds in mnar_score().")
    }
    mat_use <- mat[keep, , drop = FALSE]
    message(sprintf("Using %d / %d features with non-NA per-feature scores.",
                    sum(keep), length(keep)))
  } else {
    mat_use <- mat
  }

  n_features <- nrow(mat_use)
  n_samples  <- ncol(mat_use)

  # ── Build long-format data frame ──────────────────────────────────────────────
  #
  # One row per (feature, sample) combination:
  #   missing   : binary outcome (1 = missing, 0 = observed)
  #   group     : fixed effect — experimental condition
  #   intensity : fixed effect — per-feature mean observed intensity, used as
  #               a proxy for abundance to capture intensity-dependent missingness.
  #               Repeated for each sample since it is a feature-level summary.
  #   feature   : random effect — absorbs residual per-feature baseline missingness

  message(sprintf(
    "Building long-format data: %d features x %d samples = %d rows...",
    n_features, n_samples, n_features * n_samples
  ))

  miss_mat <- is.na(mat_use)

  # Per-feature mean observed intensity — proxy for true abundance level.
  # Scale to mean 0, sd 1 so its coefficient is on a comparable scale to group.
  feature_mean_intensity <- rowMeans(mat_use, na.rm = TRUE)
  feature_mean_intensity <- scale(feature_mean_intensity)[, 1]

  long_df <- data.frame(
    missing   = as.integer(miss_mat),
    group     = rep(group,                    each  = n_features),
    intensity = rep(feature_mean_intensity,   times = n_samples),
    feature   = rep(rownames(mat_use),        times = n_samples),
    stringsAsFactors = FALSE
  )
  long_df$feature <- factor(long_df$feature)

  # ── Fit two linear mixed models ───────────────────────────────────────────────
  #
  # Both use ML (REML = FALSE) for marginal R² calculation and to allow
  # likelihood ratio testing between models.
  #
  # mod_intensity : missing ~ intensity + (1 | feature)
  #   Baseline — captures intensity-dependent missingness only.
  #   Marginal R² reflects how much missingness is explained by abundance alone,
  #   the universal technical phenomenon present in all proteomics data.
  #
  # mod_full      : missing ~ intensity + group + (1 | feature)
  #   Full model — adds condition on top of intensity.
  #   Incremental marginal R² (full - intensity) reflects condition-specific
  #   missingness beyond what intensity already explains: the signal of interest.

  message("Fitting baseline model: missing ~ intensity + (1 | feature) ...")
  mod_intensity <- lme4::lmer(
    missing ~ intensity + (1 | feature),
    data = long_df,
    REML = FALSE
  )

  message("Fitting full model: missing ~ intensity + group + (1 | feature) ...")
  mod_full <- lme4::lmer(
    missing ~ intensity + group + (1 | feature),
    data = long_df,
    REML = FALSE
  )

  # ── Likelihood ratio test: does condition add explanatory power? ──────────────

  lrt        <- anova(mod_intensity, mod_full)
  lrt_chi2   <- lrt$Chisq[2]
  lrt_df     <- lrt$Df[2]
  lrt_pvalue <- lrt[["Pr(>Chisq)"]][2]

  # ── Extract R² from both models, handling singularity ────────────────────────
  #
  # Singularity (random effect variance = 0) affects the conditional R² but
  # not the marginal R², which is the quantity of interest. We suppress the
  # misleading performance warning and handle it explicitly.

  .extract_r2 <- function(mod) {
    singular <- lme4::isSingular(mod)
    r2_vals  <- suppressWarnings(performance::r2(mod))
    list(
      marginal    = as.numeric(r2_vals$R2_marginal),
      conditional = if (singular) NA_real_
                    else as.numeric(r2_vals$R2_conditional),
      singular    = singular
    )
  }

  r2_int       <- .extract_r2(mod_intensity)
  r2_full_vals <- .extract_r2(mod_full)

  r2_intensity_only  <- r2_int$marginal
  r2_full_marginal   <- r2_full_vals$marginal
  r2_incremental     <- r2_full_marginal - r2_intensity_only

  # Fraction of total explainable missingness that is condition-specific
  # (as opposed to intensity-driven). NA-safe.
  r2_condition_fraction <- if (r2_full_marginal > 0)
    r2_incremental / r2_full_marginal
  else
    NA_real_

  if (r2_int$singular || r2_full_vals$singular) {
    message(
      "Note: one or both models are singular (random effect variance = 0). ",
      "Marginal R² values are unaffected. Conditional R² is set to NA ",
      "where applicable."
    )
  }

  # ── Interpret ─────────────────────────────────────────────────────────────────

  interpretation <- sprintf(paste(
    "Intensity-only R²  = %.3f: proportion of missingness explained by",
    "feature abundance alone (universal intensity-dependent baseline).",
    "Incremental R²     = %.3f: additional missingness explained by",
    "experimental condition beyond intensity (condition-specific signal).",
    "Condition fraction = %.3f: %.1f%% of explainable missingness is",
    "condition-specific rather than purely intensity-driven.",
    "LRT p-value        = %s: %s that condition adds explanatory power."),
    r2_intensity_only,
    r2_incremental,
    r2_condition_fraction,
    r2_condition_fraction * 100,
    format.pval(lrt_pvalue, digits = 3),
    ifelse(lrt_pvalue < 0.05, "significant evidence", "no significant evidence")
  )

  message(sprintf(
    "Intensity-only R²: %.4f | Incremental (condition) R²: %.4f | Condition fraction: %.3f",
    r2_intensity_only, r2_incremental, r2_condition_fraction
  ))
  message(sprintf(
    "LRT: chi²(df=%d) = %.2f, p = %s",
    lrt_df, lrt_chi2, format.pval(lrt_pvalue, digits = 3)
  ))

  # ── Return ────────────────────────────────────────────────────────────────────

  invisible(list(
    # Core quantities
    r2_intensity_only     = r2_intensity_only,
    r2_incremental        = r2_incremental,
    r2_condition_fraction = r2_condition_fraction,

    # Full model marginal R² (= r2_intensity_only + r2_incremental)
    r2_full_marginal      = r2_full_marginal,

    # Conditional R² (intensity + condition + feature random effect)
    r2_full_conditional   = r2_full_vals$conditional,

    # Likelihood ratio test of condition term
    lrt_chi2              = lrt_chi2,
    lrt_df                = lrt_df,
    lrt_pvalue            = lrt_pvalue,

    # Singularity flags
    singular_intensity    = r2_int$singular,
    singular_full         = r2_full_vals$singular,

    # Models and metadata
    model_intensity       = mod_intensity,
    model_full            = mod_full,
    interpretation        = interpretation,
    n_features            = n_features,
    n_obs_total           = nrow(long_df)
  ))
}
