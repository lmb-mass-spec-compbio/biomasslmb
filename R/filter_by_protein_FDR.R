#' Filter to remove peptides from proteins failing the protein FDR threshold
#'
#' @description This function filters the content of a summarised experiment
#' based on a separate file detailing the proteins passing the FDR confidence threshold
#'
#' This function is used to filter peptides to retain only those from proteins
#' passing the FDR threshold
#'
#' @param obj `SummarisedExperiment` containing peptide-level output from Proteome Discoverer.
#' @param protein_fdr_filename `string`. Filepath for protein-level output from Proteome Discoverer
#' @param protein_col_peptide `string`. Name of column containing master
#' proteins in the peptide-level `obj`
#' @param protein_col_protein `string`. Name of column containing master
#' proteins in the protein-level `protein_fdr_filename`
#' @param protein_FDR_col `string`. Name of column containing FDR information in
#' the protein-level `protein_fdr_filename`
#' @param retain_proteins `character vector`. Vector of protein accessions to always retain,
#' even if they do not pass the FDR threshold. Default is NULL.
#' @return Returns a `SummarisedExperiment` with the filtered Proteome Discoverer output.
#' @export
filter_by_protein_fdr <- function(obj,
                                  protein_fdr_filename,
                                  protein_col_peptide = "Master.Protein.Accessions",
                                  protein_col_protein = "Accession",
                                  protein_FDR_col = 'Protein.FDR.Confidence.Combined',
                                  retain_proteins = NULL){

  check_se(obj)

  proteinData <- read.delim(file = protein_fdr_filename)

  ## Extract master protein accessions from our Peptide-level data
  proteins_in_data <- obj %>%
    rowData() %>%
    as_tibble() %>%
    select(!!sym(protein_col_peptide))

  ## Extract protein accessions and corresponding confidence from search output file
  protein_search_output <- proteinData %>%
    select(!!sym(protein_col_protein), !!sym(protein_FDR_col)) %>%
    distinct()

  ## Combine data
  protein_fdr <- left_join(
    x = proteins_in_data,
    y = protein_search_output,
    by = join_by(
      !!sym(protein_col_peptide) == !!(protein_col_protein)))

  rowData(obj)$Protein.Confidence <-
    protein_fdr$Protein.FDR.Confidence.Combined


  obj <- obj[(!is.na(rowData(obj)$Protein.Confidence) | rowData(obj)[[protein_col_peptide]] %in% retain_proteins), ]

  obj <- obj[(rowData(obj)$Protein.Confidence == 'High' | rowData(obj)[[protein_col_peptide]] %in% retain_proteins),]

  message_parse(rowData(obj), protein_col_peptide,
                'Removing peptides from proteins not passing the FDR threshold')

  return(obj)
}
