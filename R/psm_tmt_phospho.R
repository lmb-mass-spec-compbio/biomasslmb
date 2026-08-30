#' PSM-level PD output for the phospho-enriched fraction of a TMTpro experiment
#'
#' @description Proteome Discoverer PSM-level output for the phospho-enriched
#' fraction of a real TMTpro experiment in mouse fibroblasts, comparing a drug
#' treatment against a vehicle control at two timepoints (four replicates of
#' each combination), subsetted to a random selection of proteins plus a sample
#' of contaminant and non-unique-master-protein PSMs for use in the
#' `PTM site quantification` vignette. Includes the `ptmRS` columns needed by
#' \code{\link{parse_PTM_scores_pd}}. Sample columns are named to match the
#' rownames of \code{\link{tmt_phospho_design}}.
#'
#' \code{\link{psm_tmt_phospho_total}} is the matched total (non-enriched)
#' fraction of the same labelled pool.
#'
#' @keywords datasets
"psm_tmt_phospho"

#' PSM-level PD output for the total fraction matching `psm_tmt_phospho`
#'
#' @description Proteome Discoverer PSM-level output for the total
#' (non-enriched) fraction of the experiment described in
#' \code{\link{psm_tmt_phospho}}, acquired from the same labelled pool and so
#' sharing its experimental design. Used to normalise the enriched fraction and
#' to express site abundance relative to protein abundance.
#'
#' @keywords datasets
"psm_tmt_phospho_total"

#' Experimental design for `psm_tmt_phospho` and `psm_tmt_phospho_total`
#'
#' @description A `data.frame` giving the Condition (`Control`/`Treated`),
#' Timepoint (`T1`/`T2`) and Replicate for each sample in
#' \code{\link{psm_tmt_phospho}} and \code{\link{psm_tmt_phospho_total}}, with
#' rownames matching the quantification column names.
#'
#' @keywords datasets
"tmt_phospho_design"
