#' TMT data
#'
#' @description `QFeatures` object holding the single-plex TMT dataset
#' (`psm_tmt_clock`) processed through PSM QC, filtering and summarisation to
#' protein level with `sum`, as produced in the
#' `TMT workflow: PSM QC and protein summarisation` vignette.
#'
#' @keywords datasets
"tmt_qf"

#' LFQ-DDA data
#'
#' @description `QFeatures` object holding the LFQ-DDA whole-proteome dataset
#' (`lfq_dda_pd_PeptideGroups.txt`) processed through peptide QC, filtering and
#' summarisation to protein level with `robustSummary`, as produced in the
#' `LFQ-DDA workflow: peptide QC and protein summarisation` vignette.
#'
#' @keywords datasets
"lfq_qf"

#' LFQ-DIA data
#'
#' @description `QFeatures` object holding the LFQ-DIA plasma dataset
#' (`monkeypox_plasma_proteomes.parquet`) processed through precursor QC,
#' filtering and summarisation to protein level with `robustSummary`, as
#' produced in the
#' `LFQ-DIA workflow: precursor QC and protein summarisation` vignette.
#'
#' @keywords datasets
"dia_qf"

#' LFQ-DDA TurboID pulldown data
#'
#' @description `QFeatures` object holding the TurboID proximity labelling
#' dataset (`lfq_dda_pd_turboid_PeptideGroups.txt`), three biotin-treated
#' samples against three untreated controls, processed through peptide QC and
#' summarisation to protein level with `robustSummary`, as produced in Part B
#' of the `Enrichment designs: IP, BioID and TurboID` vignette.
#'
#' The `protein` assay is masked so that every value rests on at least two
#' peptides, and so retains the presence/absence structure a pulldown produces.
#' `protein_imputed` and `protein_restricted` are the blanket and restricted
#' imputations of it, the latter filling only control samples with at most one
#' quantified replicate.
#'
#' @keywords datasets
"lfq_qf_turboid"
