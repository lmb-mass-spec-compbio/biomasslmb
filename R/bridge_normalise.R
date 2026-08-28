#' Put TMT plexes onto a common scale using bridge (pooled reference) channels
#'
#' @description
#' An experiment too large for a single TMT plex is split over several plexes,
#' each of which is a separate MS run. A feature is therefore quantified from a
#' different set of PSMs in each plex, so its abundances are not directly
#' comparable between plexes and merging them confounds plex with biology.
#'
#' A bridge channel is a pooled sample, identical in every plex, which provides
#' a per-feature reference point. For each feature, this function takes the
#' difference between its abundance in a plex's bridge channel(s) and its
#' abundance across the plexes' bridge channels, and removes that difference
#' from every sample in the plex.
#'
#' @details
#' Unlike [center_normalise_to_ref()], which removes a single offset per
#' sample, the correction here is per feature *and* per plex, because that is
#' what a plex effect is: different plexes quantify a feature from different
#' PSMs. The function is agnostic to the feature level, so it can be applied to
#' a joined protein-level assay or to a joined peptide- or site-level assay.
#'
#' **Several bridge channels per plex.** These are collapsed to one reference
#' per plex with `fun`, ignoring missing values, before the plexes are compared.
#' Collapsing within a plex first matters when the plexes carry unequal numbers
#' of bridge channels, since pooling all bridge channels together would give a
#' plex with more of them more weight in defining the common scale.
#'
#' **Features without a usable bridge value in a plex.** These cannot be placed
#' on the common scale, and `on_missing` decides what happens to them. Note
#' that `on_missing='ignore'` does not assume that the plex has no offset for
#' the feature; because the reference is centred across plexes, it assumes the
#' plex's offset is the average one.
#'
#' After correction the bridge channels agree by construction, so their
#' agreement is not evidence that the correction worked. The evidence is what
#' happened to the other samples, e.g. with [plot_pca()] coloured by plex.
#'
#' @param obj `SummarizedExperiment` containing quantification for all plexes,
#' e.g. the assay produced by `QFeatures::joinAssays`.
#' @param plex_col `string`. Name of the `colData` column identifying the plex.
#' @param bridge_cols `logical` vector of length `ncol(obj)`, or `character`
#' vector of column names, identifying the bridge channels.
#' @param fun `function`. Used to collapse several bridge channels within a
#' plex to a single reference. Applied to the non-missing values only.
#' Default is `mean`.
#' @param reference `string`. Plex to align the others to. Default (`NULL`) is
#' to align every plex to the mean of the per-plex references, which preserves
#' the overall abundance scale.
#' @param on_missing `string`. What to do with a feature that has no usable
#' bridge value in a plex: `'na'` (default) sets that plex's values for the
#' feature to `NA`, `'drop'` removes the feature, and `'ignore'` leaves that
#' plex uncorrected.
#' @param min_bridge `numeric`. Minimum number of non-missing bridge values
#' required in a plex for the reference to be used. Default is 1.
#' @param on_log_scale `logical`. Input data is log-transformed, so the
#' correction is subtractive. If `FALSE`, the correction is a ratio.
#' @param verbose `logical`. Report the size of the correction and the number
#' of features affected by missing bridge values. Default is TRUE.
#' @return Returns a `SummarizedExperiment` with the plexes on a common scale.
#'
#' @importFrom SummarizedExperiment assay "assay<-" colData
#' @export
bridge_normalise <- function(obj,
                             plex_col,
                             bridge_cols,
                             fun = mean,
                             reference = NULL,
                             on_missing = c("na", "drop", "ignore"),
                             min_bridge = 1,
                             on_log_scale = TRUE,
                             verbose = TRUE){

  check_se(obj)
  check_colData_col(obj, plex_col)

  on_missing <- match.arg(on_missing)

  is_bridge <- resolve_bridge_cols(obj, bridge_cols)

  plex <- as.character(colData(obj)[[plex_col]])

  if(any(is.na(plex))){
    stop(sprintf(
      "column %s contains missing values, so %s samples cannot be assigned to a plex",
      plex_col, sum(is.na(plex))))
  }

  plexes <- if(is.factor(colData(obj)[[plex_col]])){
    intersect(levels(colData(obj)[[plex_col]]), plex)
  } else {
    unique(plex)
  }

  n_bridge_per_plex <- vapply(plexes, function(p) sum(is_bridge & plex == p),
                              numeric(1))

  if(any(n_bridge_per_plex == 0)){
    stop(sprintf("no bridge channel in plex(es): %s",
                 paste(plexes[n_bridge_per_plex == 0], collapse = ", ")))
  }

  if(length(unique(n_bridge_per_plex)) > 1){
    warning(sprintf(
      "unequal numbers of bridge channels per plex (%s), so the corrections are estimated with unequal precision",
      paste(sprintf("%s: %s", plexes, n_bridge_per_plex), collapse = "; ")))
  }

  if(!is.null(reference) && !reference %in% plexes){
    stop(sprintf("`reference` plex %s is not present in column %s. Available plexes are: %s",
                 reference, plex_col, paste(plexes, collapse = ", ")))
  }

  quant <- assay(obj)

  # Per-plex reference for each feature, and the number of bridge values it was
  # estimated from
  plex_ref <- matrix(NA_real_, nrow = nrow(quant), ncol = length(plexes),
                     dimnames = list(rownames(quant), plexes))
  n_bridge <- plex_ref
  has_data <- matrix(FALSE, nrow = nrow(quant), ncol = length(plexes),
                     dimnames = list(rownames(quant), plexes))

  for(p in plexes){
    bridge_quant <- quant[, is_bridge & plex == p, drop = FALSE]

    n_bridge[, p] <- rowSums(!is.na(bridge_quant))

    plex_ref[, p] <- apply(bridge_quant, 1, function(x){
      x <- x[!is.na(x)]
      if(length(x) < min_bridge) NA_real_ else fun(x)
    })

    has_data[, p] <- rowSums(!is.na(quant[, plex == p, drop = FALSE])) > 0
  }

  across_plex_ref <- if(is.null(reference)){
    rowMeans(plex_ref, na.rm = TRUE)
  } else {
    plex_ref[, reference]
  }

  # rowMeans of an all-NA row returns NaN
  across_plex_ref[is.nan(across_plex_ref)] <- NA_real_

  if(on_log_scale){
    correction <- plex_ref - across_plex_ref
  } else {
    correction <- plex_ref / across_plex_ref
  }

  # A plex is only a problem for a feature if the feature has data there
  uncorrectable <- is.na(correction) & has_data

  if(on_missing == "ignore"){
    correction[is.na(correction)] <- if(on_log_scale) 0 else 1
  }

  for(p in plexes){
    plex_cols <- plex == p
    if(on_log_scale){
      quant[, plex_cols] <- quant[, plex_cols, drop = FALSE] - correction[, p]
    } else {
      quant[, plex_cols] <- quant[, plex_cols, drop = FALSE] / correction[, p]
    }
  }

  assay(obj) <- quant

  if(on_missing == "drop"){
    obj <- obj[rowSums(uncorrectable) == 0, ]
  }

  if(verbose){
    report_bridge_normalise(correction, n_bridge, n_bridge_per_plex,
                            uncorrectable, on_missing, on_log_scale)
  }

  return(obj)
}


#' Resolve the bridge channels to a logical vector
#'
#' @noRd
resolve_bridge_cols <- function(obj, bridge_cols){

  if(is.logical(bridge_cols)){
    if(length(bridge_cols) != ncol(obj)){
      stop(sprintf("logical `bridge_cols` must have one value per sample (%s), not %s",
                   ncol(obj), length(bridge_cols)))
    }
    if(any(is.na(bridge_cols))){
      stop("`bridge_cols` contains missing values. A predicate such as `grepl('Pool', obj$ID)` avoids these")
    }
    is_bridge <- bridge_cols
  } else if(is.character(bridge_cols)){
    missing_cols <- setdiff(bridge_cols, colnames(obj))
    if(length(missing_cols) > 0){
      stop(sprintf("`bridge_cols` not found in the samples: %s",
                   paste(missing_cols, collapse = ", ")))
    }
    is_bridge <- colnames(obj) %in% bridge_cols
  } else {
    stop("`bridge_cols` must be a logical vector with one value per sample, or a character vector of sample names")
  }

  if(!any(is_bridge)){
    stop("`bridge_cols` identifies no bridge channels")
  }

  return(is_bridge)
}


#' Report what bridge_normalise did
#'
#' @noRd
report_bridge_normalise <- function(correction, n_bridge, n_bridge_per_plex,
                                    uncorrectable, on_missing, on_log_scale){

  plexes <- colnames(correction)

  centre <- if(on_log_scale) 0 else 1

  message(sprintf("%s features across %s plexes", nrow(correction), length(plexes)))

  for(p in plexes){
    message(sprintf(
      "plex %s: %s bridge channel(s), median absolute correction %.3f",
      p, n_bridge_per_plex[[p]],
      stats::median(abs(correction[, p] - centre), na.rm = TRUE)))
  }

  # Features whose reference came from fewer bridge channels than the plex has
  reduced <- vapply(plexes, function(p){
    sum(n_bridge[, p] > 0 & n_bridge[, p] < n_bridge_per_plex[[p]])
  }, numeric(1))

  if(sum(reduced) > 0){
    message(sprintf(
      "%s feature-plex combinations had a reference estimated from fewer bridge channels than the plex carries",
      sum(reduced)))
  }

  n_uncorrectable <- colSums(uncorrectable)

  if(sum(n_uncorrectable) > 0){
    action <- switch(on_missing,
                     na = "set to NA",
                     drop = "removed from the assay",
                     ignore = "left uncorrected")
    for(p in plexes[n_uncorrectable > 0]){
      message(sprintf(
        "plex %s: %s features have data but no usable bridge value => %s",
        p, n_uncorrectable[[p]], action))
    }
    if(on_missing == "drop"){
      message(sprintf("%s features removed in total",
                      sum(rowSums(uncorrectable) > 0)))
    }
  }
}
