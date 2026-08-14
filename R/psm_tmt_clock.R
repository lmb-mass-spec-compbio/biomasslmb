#' PSM-level PD output for a TMT12plex Control vs Mutant experiment
#'
#' @description Proteome Discoverer PSM-level output for a real TMT12plex
#' experiment comparing Control and Mutant samples (6 replicates each),
#' subsetted to ~600 randomly-selected proteins plus a sample of contaminant
#' and non-unique-master-protein PSMs for use in the
#' `TMT QC PSM-level quantification and summarisation to protein-level abundance`
#' vignette. Sample columns are named to match the rownames of
#' \code{\link{tmt_clock_design}}.
#'
#' @keywords datasets
"psm_tmt_clock"

#' Experimental design for `psm_tmt_clock`
#'
#' @description A `data.frame` giving the Condition (`Control`/`Mutant`) and
#' Replicate for each sample in \code{\link{psm_tmt_clock}}, with rownames
#' matching the quantification column names.
#'
#' @keywords datasets
"tmt_clock_design"
