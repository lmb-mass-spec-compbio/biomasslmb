#' PSM-level MaxQuant output for a TMT18plex IP vs control experiment
#'
#' @description MaxQuant `evidence.txt` PSM-level output for a real
#' TMT18plex immunoprecipitation experiment comparing a bait pulldown
#' (`IP`) with a control pulldown (`Control`), 6 replicates each, subsetted
#' to ~600 randomly-selected proteins plus a sample of contaminant and
#' non-unique-master-protein PSMs for use as the MaxQuant worked example in
#' the `Enrichment designs: IP, BioID and TurboID` vignette. Sample columns are named to match the rownames of
#' \code{\link{tmt_per2_mq_design}}. The original condition labels and the
#' `Raw.file` column (which identified the specific source project) have
#' been anonymised to generic values, and PSMs identifying the bait
#' construct itself have been excluded.
#'
#' @keywords datasets
"psm_tmt_per2_mq"

#' Experimental design for `psm_tmt_per2_mq`
#'
#' @description A `data.frame` giving the Condition (`IP`/`Control`) and
#' Replicate for each sample in \code{\link{psm_tmt_per2_mq}}, with
#' rownames matching the quantification column names.
#'
#' @keywords datasets
"tmt_per2_mq_design"

#' TMT data (MaxQuant input)
#'
#' @description `QFeatures` object containing the MaxQuant TMT IP-vs-control
#' data (\code{\link{psm_tmt_per2_mq}}) processed through PSM filtering and
#' `robustSummary` summarisation to protein level, as produced in Part A of
#' the `Enrichment designs: IP, BioID and TurboID` vignette.
#'
#' @keywords datasets
"tmt_qf_mq"
