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
#' For each feature, fits a logistic regression of binary missingness against
#' one or more experimental group variables from colData. Tjur's R²
#' (discrimination coefficient) of this model is used as a per-feature MNAR
#' score: high values indicate that missingness is structured by experimental
#' condition (MNAR-like), while values near zero indicate missingness is
#' unrelated to condition (MAR-like).
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
#'     \item{scores}{Named numeric vector of per-feature MNAR scores
#'       (Tjur's R²).}
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
#' The MNAR score for each feature is Tjur's R² (discrimination coefficient)
#' from a logistic regression:
#'
#'   missing_indicator ~ group
#'
#' where missing_indicator is 1 if the value is missing and 0 if observed,
#' and group is a factor combining the specified colData columns. Tjur's R² is
#' computed as the difference in mean predicted probabilities between missing
#' and observed values:
#'
#'   Tjur's R² = mean(P(missing | truly missing)) - mean(P(missing | truly observed))
#'
#' It ranges from 0 (no discrimination — MAR) to 1 (perfect discrimination —
#' MNAR). Logistic regression is used in preference to a linear model as it
#' is the correct model for a binary outcome.
#'
#' If multiple group_cols are provided they are combined into a single
#' interaction factor (e.g. condition:timepoint), so each unique combination
#' of levels is treated as a distinct group.
#'
#' MNAR classification thresholds (mnar_class):
#'   - "MNAR"          : Tjur's R² >= 0.5
#'   - "mixed"         : Tjur's R² in [0.2, 0.5)
#'   - "MAR"           : Tjur's R² <  0.2
#'   - "uninformative" : excluded due to min_observed / min_missing filters
#'
#' These thresholds are heuristic — inspect the score distribution before
#' applying fixed cutoffs to your dataset.
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
#'      xlab = "MNAR score (Tjur's R²)", main = "Missingness structure")
#' print(result$global)
#' head(result$summary)
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

  miss_mat     <- is.na(mat)
  n_separation <- 0L

  scores <- vapply(seq_len(n_features), function(idx) {
    m      <- miss_mat[idx, ]
    n_miss <- sum(m)
    n_obs  <- sum(!m)
    if (n_obs < min_observed || n_miss < min_missing) return(NA_real_)

    # Perfect separation: missingness is fully determined by group
    # (e.g. always missing in condition A, never in condition B).
    # glm will not converge in this case; Tjur's R² = 1 is correct.
    per_group_var <- tapply(as.numeric(m), group, stats::var)
    if (all(per_group_var == 0, na.rm = TRUE)) return(1.0)

    df  <- data.frame(missing = as.integer(m), group = group)
    mod <- withCallingHandlers(
      stats::glm(missing ~ group, data = df, family = stats::binomial()),
      warning = function(w) {
        if (grepl("fitted probabilities numerically 0 or 1", conditionMessage(w),
                  fixed = TRUE)) {
          n_separation <<- n_separation + 1L
          invokeRestart("muffleWarning")
        }
      }
    )
    preds <- stats::fitted(mod)
    mean(preds[m]) - mean(preds[!m])  # Tjur's R²
  }, numeric(1))

  if (n_separation > 0L) {
    message(sprintf(
      "%d feature(s) showed near-complete separation (missingness almost fully determined by group). Tjur's R\u00b2 for these features will be close to 1 \u2014 this indicates strong MNAR structure and is expected.",
      n_separation
    ))
  }

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


#' Compute a single dataset-level MNAR index using logistic regression
#'
#' Summarises missingness structure across the whole dataset into a single
#' value by fitting two pooled logistic regressions across all
#' (feature x sample) observations simultaneously:
#'
#'   mod_intensity : missing ~ intensity
#'   mod_full      : missing ~ intensity + group
#'
#' The incremental Tjur's R² (full - intensity) is the primary quantity of
#' interest: it measures condition-specific missingness beyond what feature
#' abundance alone would predict.
#'
#' @param obj A QFeatures object.
#' @param i Integer or character. Index or name of the assay to analyse.
#' @param group_cols Character vector. One or more column names from colData(obj)
#'   that define experimental groups (e.g. c("condition", "treatment")).
#'   Multiple columns are combined into an interaction term.
#' @param min_observed Integer. Minimum number of observed (non-missing) values
#'   required for a feature to be included in the model. Default: 2.
#' @param min_missing Integer. Minimum number of missing values required for a
#'   feature to be included in the model. Default: 1.
#' @param subset_features Logical. If TRUE (default), features with fewer than
#'   \code{min_observed} observed values or fewer than \code{min_missing} missing
#'   values are excluded, as they are uninformative and would dilute the
#'   intensity and group effect estimates.
#' @param log_transform Logical. If TRUE, apply log2 to the assay intensities
#'   before computing the per-feature mean intensity covariate. Set to TRUE if
#'   your data has not already been log-transformed. Default: FALSE.
#'
#' @return A list with the following elements:
#'   \describe{
#'     \item{tjur_intensity_only}{Numeric. Tjur's R² of the intensity-only
#'       model. The discrimination between missing and observed explained by
#'       feature abundance alone — the universal intensity-dependent baseline
#'       present in all proteomics data.}
#'     \item{tjur_incremental}{Numeric. The increase in Tjur's R² when
#'       condition is added on top of intensity. This is the primary quantity
#'       of interest: condition-specific missingness beyond what intensity
#'       already explains. Near zero means missingness is purely
#'       intensity-driven; positive values indicate condition-specific
#'       depletion of particular features.}
#'     \item{tjur_condition_fraction}{Numeric. The fraction of total
#'       explainable missingness (tjur_full) that is condition-specific rather
#'       than intensity-driven. Ranges from 0 (all intensity-driven) to 1 (all
#'       condition-specific).}
#'     \item{tjur_full}{Numeric. Tjur's R² of the full model
#'       (intensity + condition). Equal to tjur_intensity_only +
#'       tjur_incremental.}
#'     \item{lrt_chi2}{Numeric. Chi-squared statistic from the likelihood
#'       ratio test of the full vs intensity-only model.}
#'     \item{lrt_df}{Integer. Degrees of freedom for the LRT.}
#'     \item{lrt_pvalue}{Numeric. P-value for the LRT — tests whether
#'       condition adds significant explanatory power over intensity alone.}
#'     \item{model_intensity}{The fitted intensity-only \code{glm} object.}
#'     \item{model_full}{The fitted full \code{glm} object.}
#'     \item{interpretation}{Character string. Human-readable summary of all
#'       key quantities.}
#'     \item{n_features}{Integer. Number of features included in the models.}
#'     \item{n_obs_total}{Integer. Total number of (feature x sample) rows.}
#'   }
#'
#' @details
#' Tjur's R² (discrimination coefficient) is computed as:
#'
#'   mean(P(missing | truly missing)) - mean(P(missing | truly observed))
#'
#' It ranges from 0 (no discrimination) to 1 (perfect discrimination) and is
#' a natural measure of how well the model separates missing from observed.
#'
#' Pooled logistic regression (no random effects) is used in preference to a
#' linear mixed model. A per-feature random intercept would absorb the very
#' variance that intensity is meant to explain (low-abundance features have
#' both low intensity and high baseline missingness), suppressing the intensity
#' fixed effect. The pooled approach is analogous to the detection probability
#' curve (DPC) fitted by the limpa package.
#'
#' The intensity covariate is the per-feature mean observed intensity, scaled
#' to mean 0 and sd 1. If your data has not been log-transformed upstream,
#' set \code{log_transform = TRUE} so that the sigmoidal intensity-missingness
#' relationship is approximately linear on the log scale, consistent with how
#' DPC curves are typically modelled.
#'
#' The decomposition is:
#'
#'   tjur_full = tjur_intensity_only + tjur_incremental
#'
#' A large tjur_intensity_only with small tjur_incremental means missingness
#' is driven by abundance uniformly across conditions — the typical baseline.
#' A substantial tjur_incremental means certain features are specifically
#' depleted in particular conditions beyond what their overall abundance
#' would predict.
#'
#' @references
#' Tjur T (2009). Coefficients of determination in logistic regression models —
#' a new proposal: the coefficient of discrimination.
#' The American Statistician, 63(4), 366-372.
#'
#' @examples
#' \dontrun{
#' idx <- mnar_global_score(obj, i = "peptides", group_cols = "condition")
#'
#' # Primary quantities
#' cat("Intensity-only Tjur R²:  ", round(idx$tjur_intensity_only,  3), "\n")
#' cat("Incremental (condition): ", round(idx$tjur_incremental,     3), "\n")
#' cat("Condition fraction:      ", round(idx$tjur_condition_fraction, 3), "\n")
#' cat("LRT p-value:             ", format.pval(idx$lrt_pvalue),        "\n")
#' cat(idx$interpretation, "\n")
#' }
#'
#' @seealso \code{\link{mnar_score}}
#' @export
mnar_global_score <- function(obj,
                              i,
                              group_cols,
                              min_observed    = 2,
                              min_missing     = 1,
                              subset_features = TRUE,
                              log_transform   = FALSE) {

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

  # ── Extract assay matrix and group factor ─────────────────────────────────────

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
    n_obs  <- rowSums(!is.na(mat))
    n_miss <- rowSums(is.na(mat))
    keep   <- n_obs >= min_observed & n_miss >= min_missing
    if (sum(keep) == 0) {
      stop("No informative features found. ",
           "Check min_observed / min_missing thresholds.")
    }
    mat_use <- mat[keep, , drop = FALSE]
    message(sprintf("Using %d / %d features meeting min_observed / min_missing thresholds.",
                    sum(keep), nrow(mat)))
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
  #   intensity : fixed effect — per-feature mean observed intensity, scaled to
  #               mean 0 sd 1. Repeated for each sample as it is a feature-level
  #               summary. Log2-transformed first if log_transform = TRUE.

  message(sprintf(
    "Building long-format data: %d features x %d samples = %d rows...",
    n_features, n_samples, n_features * n_samples
  ))

  miss_mat <- is.na(mat_use)

  feature_mean_intensity <- rowMeans(mat_use, na.rm = TRUE)
  if (log_transform) feature_mean_intensity <- log2(feature_mean_intensity)
  feature_mean_intensity <- scale(feature_mean_intensity)[, 1]

  long_df <- data.frame(
    missing   = as.integer(miss_mat),
    group     = rep(group,                   each  = n_features),
    intensity = rep(feature_mean_intensity,  times = n_samples),
    stringsAsFactors = FALSE
  )

  # ── Fit two pooled logistic regressions ──────────────────────────────────────
  #
  # Pooled (no random effects) to avoid the per-feature random intercept
  # absorbing the intensity signal it is meant to measure. This is analogous
  # to how DPC curves are fitted.
  #
  # mod_intensity : missing ~ intensity
  #   Baseline — how much missingness is explained by feature abundance alone.
  #
  # mod_full      : missing ~ intensity + group
  #   Adds condition. Incremental Tjur R² is the condition-specific signal.

  .catch_separation <- function(w) {
    if (grepl("fitted probabilities numerically 0 or 1", conditionMessage(w),
              fixed = TRUE)) {
      message("Note: some features show near-complete separation between missing ",
              "and observed values. Fitted probabilities close to 0 or 1 are ",
              "expected for strong MNAR features and do not affect Tjur's R\u00b2.")
      invokeRestart("muffleWarning")
    }
  }

  message("Fitting baseline model: missing ~ intensity ...")
  mod_intensity <- withCallingHandlers(
    stats::glm(missing ~ intensity, data = long_df, family = stats::binomial()),
    warning = .catch_separation
  )

  message("Fitting full model: missing ~ intensity + group ...")
  mod_full <- withCallingHandlers(
    stats::glm(missing ~ intensity + group, data = long_df, family = stats::binomial()),
    warning = .catch_separation
  )

  # ── Likelihood ratio test: does condition add explanatory power? ──────────────

  lrt        <- anova(mod_intensity, mod_full, test = "Chisq")
  lrt_chi2   <- lrt$Deviance[2]
  lrt_df     <- lrt$Df[2]
  lrt_pvalue <- lrt[["Pr(>Chi)"]][2]

  # ── Tjur's R² for both models ─────────────────────────────────────────────────
  #
  # Tjur's R² = mean(P(missing | truly missing)) - mean(P(missing | truly observed))
  # Ranges from 0 (no discrimination) to 1 (perfect discrimination).

  .tjur_r2 <- function(mod, missing_vec) {
    preds <- stats::fitted(mod)
    mean(preds[missing_vec == 1L]) - mean(preds[missing_vec == 0L])
  }

  tjur_intensity_only <- .tjur_r2(mod_intensity, long_df$missing)
  tjur_full           <- .tjur_r2(mod_full,      long_df$missing)
  tjur_incremental    <- tjur_full - tjur_intensity_only

  tjur_condition_fraction <- if (tjur_full > 0)
    tjur_incremental / tjur_full
  else
    NA_real_

  # ── Interpret ─────────────────────────────────────────────────────────────────

  interpretation <- sprintf(paste(
    "Intensity-only Tjur R²  = %.3f: discrimination between missing and",
    "observed explained by feature abundance alone (intensity-dependent baseline).",
    "Incremental Tjur R²     = %.3f: additional discrimination explained by",
    "experimental condition beyond intensity (condition-specific signal).",
    "Condition fraction       = %.3f: %.1f%% of explainable missingness is",
    "condition-specific rather than purely intensity-driven.",
    "LRT p-value              = %s: %s that condition adds explanatory power."),
    tjur_intensity_only,
    tjur_incremental,
    tjur_condition_fraction,
    tjur_condition_fraction * 100,
    format.pval(lrt_pvalue, digits = 3),
    ifelse(lrt_pvalue < 0.05, "significant evidence", "no significant evidence")
  )

  message(sprintf(
    "Intensity-only Tjur R²: %.4f | Incremental (condition) Tjur R²: %.4f | Condition fraction: %.3f",
    tjur_intensity_only, tjur_incremental, tjur_condition_fraction
  ))
  message(sprintf(
    "LRT: chi²(df=%d) = %.2f, p = %s",
    lrt_df, lrt_chi2, format.pval(lrt_pvalue, digits = 3)
  ))

  # ── Return ────────────────────────────────────────────────────────────────────

  invisible(list(
    tjur_intensity_only     = tjur_intensity_only,
    tjur_incremental        = tjur_incremental,
    tjur_condition_fraction = tjur_condition_fraction,
    tjur_full               = tjur_full,
    lrt_chi2                = lrt_chi2,
    lrt_df                  = lrt_df,
    lrt_pvalue              = lrt_pvalue,
    model_intensity         = mod_intensity,
    model_full              = mod_full,
    interpretation          = interpretation,
    n_features              = n_features,
    n_obs_total             = nrow(long_df)
  ))
}
