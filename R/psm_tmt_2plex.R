#' PSM-level PD output for a two-plex TMTpro 18plex experiment
#'
#' @description Proteome Discoverer PSM-level output for a real experiment
#' spread over two TMTpro 18plex plexes, each carrying a pooled bridge
#' channel. The design is a 3 x 3 factorial: three cell lines (a wild type and
#' two knockouts) crossed with three drug treatment timepoints, three
#' replicates each. Subsetted to ~300 randomly-selected proteins quantified in
#' both plexes, plus a sample of contaminant and non-unique-master-protein
#' PSMs, for use in Part C of the
#' `TMT QC PSM-level quantification and summarisation to protein-level abundance`
#' vignette.
#'
#' The `Plex` column identifies which plex each PSM came from, and the
#' abundance columns are named by TMT tag, since the same tag denotes a
#' different sample in each plex. \code{\link{tmt_2plex_design}} maps tag and
#' plex to sample.
#'
#' @keywords datasets
"psm_tmt_2plex"

#' Experimental design for `psm_tmt_2plex`
#'
#' @description A `data.frame` with one row per labelled channel, giving the
#' Genotype, Time and Replicate for each sample in
#' \code{\link{psm_tmt_2plex}}. Channels that were not labelled in a given plex
#' are absent. `quantCols` gives the TMT tag and `runCol` the plex, as
#' required by `QFeatures::readQFeatures(runCol = 'Plex')`. The pooled bridge
#' channel present in both plexes has a `Genotype` of `Bridge`.
#'
#' @keywords datasets
"tmt_2plex_design"
