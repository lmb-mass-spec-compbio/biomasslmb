#' Read in data from DIA-NN
#'
#' @description This function reads in the report.tsv file from DIA-NN, performs
#' filtering using the Q-Value columns and then joins together the quantification
#' from the individual samples
#'
#' @param diann_report_infile `character`. Filepath to DIA-NN report.tsv file
#' @param global_protein_q `numeric`. Threshold for global protein FDR
#' @param run_protein_q `numeric`. Threshold for run-specific protein FDR
#' @param run_precursor_q `numeric`. Threshold for run-specific precursor FDR
#' @param return_sep_quant `logical`. Return the individual sample-level quantification as well
#' @export
readDIANNFilterQJoin <- function(diann_report_infile,
                                 global_protein_q = 0.01,
                                 run_protein_q = 0.01,
                                 run_precursor_q = 0.01,
                                 return_sep_quant=FALSE){

  # Read in the DIA-NN data at precursor level
  qf <- readQFeaturesFromDIANN(assayData = read.delim(diann_report_infile),
                               quantCols = "Precursor.Quantity",
                               runCol = "Run",
                               fnames = "Precursor.Id")

  # Precursor filtering using Q.Value columns. Consider adding (optional) PEP filtering here too
  qf <- qf %>%
    filterFeatures(
      VariableFilter(field='Lib.PG.Q.Value', value=global_protein_q, condition='<=')) %>%
    filterFeatures(
      VariableFilter(field='PG.Q.Value', value=run_protein_q, condition='<=')) %>%
    filterFeatures(
      VariableFilter(field='Q.Value', value=run_precursor_q, condition='<='))

  # Each Precursor.Id value occurs multiple times across the samples and thus
  # gets a numerical identifier added to the end.
  # We need to remove this before we join the assays
  n_samples <- length(qf)
  for (i in 1:n_samples){
    rownames(qf[[i]]) <- gsub("\\.\\d+", "", rownames(qf[[i]]))
  }



  # Now we can join the assays as the rownames for precursors should be consistent.

  qf<- joinAssays(x = qf,
                  i = 1:n_samples,
                  name = "peptides_fdr_cntrl")

  if(return_sep_quant){
    # return separate and joined quantification
    return(qf)
  } else{
    # Return just the FDR controlled peptide(precursor)-level abundances
    return(qf[,,'peptides_fdr_cntrl'])
  }
}
