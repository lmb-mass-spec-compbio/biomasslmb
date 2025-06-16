


#' Adds new feature describing the average reporter Signal/Noise ratio.
#'
#' @description PD column `Average.Reporter.SN` is NA when all tags have missing
#' values and imputes a value of zero for NA when there not all tags are missing.
#' So, with a 10-plex TMT experiment, if a single tag has a SN of 10 and all others
#' are NA, PD will report an average SN of 10/10 = 1.
#'
#' This function adds an average reporter SN column which ignores missing
#' values. In the example above, the average SN value reported would be 10.
#' Where all values are NA, the average SN will remain NA.
#'
#' The assumption is that the intensity values in the `exprs` matrix
#' are Signal/Noise ratios, as reported by PD by default
#'
#' @param obj `SummarizedExperiment`. Should contain PSMs-level TMT quantification
#' @param sn_col `string`. Name of output column containing the average signal:noise
#'
#' @return Returns an `MSnSet` with the average SN included as a new feature column.
#' @export
#' @importFrom SummarizedExperiment rowData assay
update_average_sn <- function(obj,
                              sn_col='Average.Reporter.SN'){

  check_se_psm(obj)

  rowData(obj)[[sn_col]] <- rowMeans(assay(obj), na.rm=TRUE)

  return(obj)
}



#' Filter a PSM-level summarizedExperiment to remove low quality PSMs
#'
#' @description Filter PSMs from TMT quantification to remove the following:
#'
#' 1. Missing values (`NA`) for all tags
#' 2. Interference/co-isolation above a set value (default=100, e.g no filtering)
#' 3. Signal:noise ratio below a set value (default=0, e.g no filtering)
#' 4. Quan.Info is not empty ('')
#' 5. PSM.Ambiguity is not 'Selected' or 'Unambiguous'
#'
#' @param obj `summarizedExperiment`. Should contain PSMs-level TMT quantification
#' @param inter_thresh `numeric`. Maximum allowed interference/co-isolation
#' @param sn_thresh `numeric`. Minimum allowed average signal:noise
#' @param master_protein_col `string`. Name of column containing master
#' proteins.
#' @param inter_col `string`. Name of column containing the interference value.
#' @param sn_col `string`. Name of column containing the signal:noise value.
#' @param from_PD `logical`. Input is from ProteomeDiscover. If set, will filter using
#' Quan.Info and PSM.Ambiguity columns. Default is TRUE
#' @param verbose `logical`. Default is TRUE, use verbose output messages.
#'
#' @return Returns an `summarizedExperiment` with the filtered PSMs.
#' @export
filter_TMT_PSMs <- function(obj,
                            inter_thresh=100,
                            sn_thresh=0,
                            master_protein_col='Master.Protein.Accessions',
                            inter_col='Isolation.Interference.in.Percent',
                            sn_col='Average.Reporter.SN',
                            from_PD=TRUE,
                            verbose=TRUE){

  check_se_psm(obj)

  if(verbose) message("Filtering PSMs...")
  if(verbose) message_parse(rowData(obj),
                            master_protein_col,
                            "Initial PSMs")

  if(from_PD){
    obj <- obj[rowData(obj)$Quan.Info=='',]
    if(verbose) message_parse(rowData(obj),
                              master_protein_col,
                              "PSMs with Quan.Info removed")

    obj <- obj[rowData(obj)$PSM.Ambiguity %in% c('Selected', 'Unambiguous'),]
    if(verbose) message_parse(rowData(obj),
                              master_protein_col,
                              "PSMs which are not selected or unambiguous removed")
  }

  obj <- obj[rowSums(is.finite(assay(obj)))>0,]
  if(verbose) message_parse(rowData(obj),
                            master_protein_col,
                            "Removing PSMs without quantification values")

  if(inter_thresh<100){
    obj <- obj[rowData(obj)[[inter_col]]<=inter_thresh,]
    if(verbose) message_parse(rowData(obj),
                              master_protein_col,
                              "Removing PSMs with high Co-isolation/interference")
  } else{
    if(verbose) message(sprintf('Not performing filtering by interference thresholding (`inter_thresh`=%s)',
                                inter_thresh))
  }

  if(sn_thresh>0){

    obj <- obj[rowData(obj)[[sn_col]]>=sn_thresh,]

    if(verbose) message_parse(rowData(obj),
                              master_protein_col,
                              "Removing PSMs with low average S:N ratio")
  } else{
    if(verbose) message(sprintf('Not performing filtering by average S:N ratio (`sn_thresh`=%s)',
                                sn_thresh))
  }

  return(obj)

}
