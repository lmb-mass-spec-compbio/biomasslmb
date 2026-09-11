#' PSM-level MaxQuant output for a TMT18plex whole-proteome factorial experiment
#'
#' @description MaxQuant `evidence.txt` PSM-level output for a real TMT18plex
#' whole-proteome experiment with a crossed design: three cell lines
#' (`Line_A`, `Line_B` and the double `Line_AB`) each treated with a vehicle
#' (`Control`) or a compound (`Treated`), in 3 replicates. Subsetted to ~1500
#' randomly-selected proteins plus samples of contaminant, decoy and
#' non-unique-master-protein PSMs, and used as the MaxQuant worked example in
#' the `TMT workflow: PSM QC and protein summarisation` vignette and as the
#' factorial worked example in
#' `Data exploration and statistical testing`. Sample columns are named to
#' match the rownames of \code{\link{tmt_factorial_design}}.
#'
#' The proteins were sampled at random rather than chosen for effect size, so
#' the number of differentially abundant proteins reported in the vignettes is
#' what a random slice of this experiment gives, not what the full experiment
#' gives.
#'
#' The design labels, the `Raw.file` column and PSMs identifying the tagged
#' constructs themselves have been anonymised or excluded, since those
#' identify the source project rather than describing the experiment. Columns
#' cross-referencing MaxQuant output tables that are not packaged here, and
#' columns constant across every retained PSM, have been dropped.
#'
#' @keywords datasets
"psm_tmt_factorial"

#' Experimental design for `psm_tmt_factorial`
#'
#' @description A `data.frame` giving the `Genotype` (`Line_A`, `Line_B` or
#' `Line_AB`), `Treatment` (`Control` or `Treated`) and `Replicate` for each
#' sample in \code{\link{psm_tmt_factorial}}, with rownames matching the
#' quantification column names. The two factors are fully crossed, with 3
#' replicates in each of the 6 combinations.
#'
#' @keywords datasets
"tmt_factorial_design"

#' Protein-level abundances for the TMT factorial experiment
#'
#' @description `QFeatures` object containing the protein-level abundances
#' produced from \code{\link{psm_tmt_factorial}} by the PSM filtering and
#' summarisation in the `TMT workflow: PSM QC and protein summarisation`
#' vignette, read back in by
#' `Data exploration and statistical testing`. Only the `protein` assay is
#' retained: the PSM-level assays it was built from are available as
#' \code{\link{psm_tmt_factorial}}, and keeping them here would make the
#' object an order of magnitude larger.
#'
#' @keywords datasets
"tmt_qf_factorial"
