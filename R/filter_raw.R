#' @noRd
remove_no_quant_assay <- function(obj, master_protein_col){
  # Remove features with no quantification according to the assay data
  obj <- obj[rowSums(is.finite(assay(obj)))>0,]
  message_parse(rowData(obj), master_protein_col, "features without quantification removed")
  return(obj)
}

#' @noRd
remove_non_unique_master_protein <- function(obj, master_protein_col){
  # Remove features without a unique master protein annotation
  obj <- obj[rowData(obj)[["Number.of.Protein.Groups"]] == 1, ]
  message_parse(rowData(obj), master_protein_col, "features with non-unique master proteins removed")
  return(obj)
}


#' @noRd
remove_no_master <- function(obj, master_protein_col){
  # remove features without a master protein
  missing_master_protein <- (is.na(rowData(obj)[[master_protein_col]]) |
                               rowData(obj)[[master_protein_col]] == '')

  obj <- obj[!missing_master_protein,]
  message_parse(rowData(obj), master_protein_col, "features without a master protein removed")

  return(obj)
}


#' @noRd
remove_contaminant <- function(obj,
                               contaminant_proteins,
                               filter_associated_contaminant,
                               master_protein_col,
                               protein_col=NULL){

  if (is.null(contaminant_proteins)) {
    stop("must supply the contaminant_proteins argument to filter contaminant proteins")
  }

  message(sprintf(
    "%s contaminant proteins supplied",
    length(contaminant_proteins)
  ))

  # Identify features for contaminant_proteins
  if(!is.null(protein_col)){
    contaminant_features <- rowData(obj)[[master_protein_col]] %in% contaminant_proteins |
      grepl("contaminant", rowData(obj)[[protein_col]])
  } else{
    contaminant_features <- rowData(obj)[[master_protein_col]] %in% contaminant_proteins
  }


  # identify associated contaminant proteins first if necessary
  if (filter_associated_contaminant) {

    if (is.null(protein_col)) {
      stop("must supply the protein_col argument to identify and remove associated contaminants")
    }

    associated_contaminant <- rowData(obj[contaminant_features,])[[protein_col]] %>%
      strsplit("; ") %>%
      unlist()

    associated_contaminant <- associated_contaminant[!grepl("cRAP", associated_contaminant)]

    message(sprintf(
      "%s proteins identified as 'contaminant associated'",
      length(associated_contaminant)
    ))
  }

  # then remove normal contaminant proteins
  obj <- obj[!contaminant_features,]
  message_parse(rowData(obj), master_protein_col, "contaminant features removed")

  # then remove associated contaminant proteins if necessary
  if (filter_associated_contaminant) {
    if (length(associated_contaminant) > 0) {
      # remove isoforms
      associated_contaminant_no_isoform <- unique(sapply(strsplit(associated_contaminant, "-"), "[[", 1))
      associated_contaminant_regex <- paste(associated_contaminant_no_isoform, collapse = "|")
      obj <- obj[!grepl(associated_contaminant_regex, rowData(obj)[[protein_col]]), ]

      message_parse(rowData(obj), master_protein_col, "associated contaminant features removed")
    }
  }

  return(obj)
}

#' Filter Proteome Discoverer DDA output
#'
#' @description This function filters the output .txt files (peptide groups or PSMs) from
#' Proteome Discoverer for DDA, based on various criteria:
#'
#' 1. Remove features without a master protein
#' 2. Remove features without a unique master protein  (i.e.
#'    Number.of.Protein.Groups == 1)
#' 3. Remove features matching a contaminant protein
#' 4. Remove features matching any protein associated with
#'    a contaminant protein (see below)
#' 5. Remove features without quantification values
#'
#' @details **Associated contaminant proteins** are proteins which have at least
#' one feature shared with a contaminant protein. It has been observed that the contaminant
#' fasta files often do not contain all possible contaminant proteins e.g. some features
#' can be assigned to a keratin which is not in the provided contaminant database.
#'
#' In the example below, using `filter_associated_contaminant = TRUE` will filter out f2 and f3 in
#' addition to f1, regardless of the value in the Master.Protein.Accession column.
#'
#' ```
#' feature  Protein.Accessions         Master.Protein.Accessions
#' f1       protein1, protein2, contaminant,  protein1,
#' f2       protein1, protein3         protein3,
#' f3       protein2                   protein2
#' ```
#' @param obj `SummarisedExperiment` containing output from Proteome Discoverer.
#' Use \code{\link[QFeatures]{readQFeatures}} to read in .txt file
#' @param master_protein_col `string`. Name of column containing master
#' proteins.
#' @param protein_col `string`. Name of column containing all protein
#' matches.
#' @param unique_master `logical`. Filter out features without a unique
#' master protein.
#' @param filter_contaminant `logical`. Filter out features which match a contaminant
#' protein.
#' @param contaminant_proteins `character vector`. The protein IDs form the contaminant proteins
#' @param crap_proteins `character vector`. Same as contaminant_proteins. Available for backwards compatibility. Default is NULL. If both contaminant_proteins and crap_proteins are set, an error is thrown.
#' @param filter_associated_contaminant `logical`. Filter out features which
#' match a contaminant associated protein.
#' @param remove_no_quant `logical`. Remove features with no quantification
#' @return Returns a `SummarisedExperiment` with the filtered Proteome Discoverer output.
#' @examples
#' \dontrun{
#'
#' #### PSMs.txt example ####
#' # load PD PSMs.txt output
#' tmt_qf <- readQFeatures(assayData = psm_tmt_total,
#'  quantCols = 36:45,
#'  name = "psms_raw")
#'
#' # extract the UniProt accessions from the contaminant FASTA headers
#' contaminant_accessions <- get_contaminant_fasta_accessions(contaminant_fasta_inf)
#'
#' # filter the PSMs
#' psm2 <- filter_features_pd_dda(
#'   obj = tmt_qf[['psms_raw']],
#'   master_protein_col = "Master.Protein.Accessions",
#'   protein_col = "Protein.Accessions",
#'   unique_master = TRUE,
#'   TMT = TRUE,
#'   filter_contaminant = TRUE,
#'   contaminant_proteins = contaminant_accessions,
#'   filter_associated_contaminant = TRUE
#' )
#'
#'
#' }
#' @export
filter_features_pd_dda <- function(obj,
                                   master_protein_col = "Master.Protein.Accessions",
                                   protein_col = "Protein.Accessions",
                                   unique_master = TRUE,
                                   filter_contaminant = TRUE,
                                   contaminant_proteins = NULL,
                                   crap_proteins = NULL,
                                   filter_associated_contaminant = TRUE,
                                   remove_no_quant = TRUE){

  # Check use of contaminants/crap arguments
  if(!is.null(crap_proteins)){
    if(!is.null(contaminant_proteins)){
      stop(paste0('Both crap_proteins and contaminant_proteins arguments were ',
                  'provided. They are the same, with the former only being provided ',
                  'for backwards compatibility. Please use only contaminant_proteins'))

    } else{
      contaminant_proteins <- crap_proteins
    }
  }

  # check arguments
  check_se(obj)

  if('Identifying.Node' %in% colnames(rowData(obj))){
    if(length(unique(rowData(obj)$Identifying.Node))>1){
      message('Results from multiple search engines detected. See documentation for
              update_peptide_assignments() and remove_redundant_psm_quant()
              functions to handle this')
    }
  }

  # print input summary
  message("Filtering data...")
  message_parse(rowData(obj), master_protein_col, "Input")

  # remove contaminant proteins
  if (filter_contaminant) {

    obj <- remove_contaminant(obj,
                       contaminant_proteins,
                       filter_associated_contaminant,
                       master_protein_col,
                       protein_col)

    if('Contaminant' %in% colnames(rowData(obj))){
      obj <- obj[rowData(obj)$Contaminant=='False',]
      message_parse(rowData(obj), master_protein_col, "PD-labelled 'Contaminants' removed")
    }
  }

  obj <- remove_no_master(obj, master_protein_col)

  # remove features with non-unique master proteins
  if (unique_master) {
        # Add option(s) here to handle cases where column does not exist
    obj <- remove_non_unique_master_protein(obj, master_protein_col)
  }

  # remove features with quantification warnings if necessary
  if (remove_no_quant) {
    obj <- remove_no_quant_assay(obj, master_protein_col)
  }

  return(obj)

}


#' Filter DIA-NN output
#'
#' @description This function filters the precursor .txt file from
#' DIA-NN, based on various criteria:
#'
#' 1. Remove features without a master protein (Protein.Group column)
#' 2. Remove features without a unique master protein
#' 3. Remove features matching a contaminant protein
#' 4. Remove features matching any protein associated with
#'    a contaminant protein (see below)
#' 5. Remove features without quantification values
#'
#' @details **Associated contaminant proteins** are proteins which have at least
#' one feature shared with a contaminant protein. It has been observed that the contaminant
#' fasta files often do not contain all possible contaminant proteins e.g. some features
#' can be assigned to a keratin which is not in the provided contaminant database.
#'
#' In the example below, using `filter_associated_contaminant = TRUE` will filter out f2 and f3 in
#' addition to f1, regardless of the value in the Master.Protein.Accession column.
#'
#' ```
#' feature  Protein.Accessions         Master.Protein.Accessions
#' f1       protein1, protein2, contaminant,  protein1,
#' f2       protein1, protein3         protein3,
#' f3       protein2                   protein2
#' ```
#' @param obj `SummarisedExperiment` containing output from Proteome Discoverer.
#' Use \code{\link[QFeatures]{readQFeatures}} to read in .txt file
#' @param master_protein_col `string`. Name of column containing master
#' proteins.
#' @param protein_col `string`. Name of column containing all protein
#' matches.
#' @param unique_master `logical`. Filter out features without a unique
#' master protein.
#' @param filter_contaminant `logical`. Filter out features which match a contaminant
#' protein.
#' @param contaminant_proteins `character vector`. The protein IDs form the contaminant proteins
#' @param filter_associated_contaminant `logical`. Filter out features which
#' match a contaminant associated protein.
#' @param remove_no_quant `logical`. Remove features with no quantification
#' @return Returns a `SummarisedExperiment` with the filtered Proteome Discoverer output.
#' @export
filter_features_diann <- function(obj,
                                 master_protein_col = "Protein.Group",
                                 protein_col = "Protein.Ids",
                                 unique_master = TRUE,
                                 filter_contaminant = TRUE,
                                 contaminant_proteins = NULL,
                                 filter_associated_contaminant = TRUE,
                                 remove_no_quant = TRUE){


  # check arguments
  check_se(obj)

  # print input summary
  message("Filtering data...")
  message_parse(rowData(obj), master_protein_col, "Input")

  # remove contaminant proteins
  if (filter_contaminant) {

    obj <- remove_contaminant(obj,
                       contaminant_proteins,
                       filter_associated_contaminant,
                       master_protein_col,
                       protein_col)

  }

  obj <- remove_no_master(obj, master_protein_col)

  # remove features with non-unique master proteins
  if (unique_master) {
    # First, add a column with the number of master proteins
    rowData(obj)[["Number.of.Protein.Groups"]] <- sapply(
      strsplit(rowData(obj)[[master_protein_col]], split=';'), length)

    obj <- remove_non_unique_master_protein(obj, master_protein_col)
  }

  # remove features with quantification warnings if necessary
  if (remove_no_quant) {
    obj <- remove_no_quant_assay(obj, master_protein_col)
  }

  return(obj)

}


#' Filter Spectronaut output
#'
#' @description This function filters the precursor .txt file from
#' Spectronaut, based on various criteria:
#'
#' 1. Remove features without a master protein (PG.ProteinAccessions column)
#' 2. Remove features without a unique master protein
#' 3. Remove features matching a contaminant protein
#' 5. Remove features without quantification values
#'
#' @param obj `SummarisedExperiment` containing output from Proteome Discoverer.
#' Use \code{\link[QFeatures]{readQFeatures}} to read in .txt file
#' @param master_protein_col `string`. Name of column containing master
#' proteins.
#' @param unique_master `logical`. Filter out features without a unique
#' master protein.
#' @param filter_contaminant `logical`. Filter out features which match a contaminant
#' protein.
#' @param contaminant_proteins `character vector`. The protein IDs form the contaminant proteins
#' @param remove_no_quant `logical`. Remove features with no quantification
#' @return Returns a `SummarisedExperiment` with the filtered Proteome Discoverer output.
#' @export
filter_features_sn <- function(obj,
                                 master_protein_col = "PG.ProteinAccessions",
                                 unique_master = TRUE,
                                 filter_contaminant = TRUE,
                                 contaminant_proteins = NULL,
                                 remove_no_quant = TRUE){


  # check arguments
  check_se(obj)

  # print input summary
  message("Filtering data...")
  message_parse(rowData(obj), master_protein_col, "Input")

  # remove contaminant proteins
  if (filter_contaminant) {

    obj <- remove_contaminant(obj,
                       contaminant_proteins,
                       filter_associated_contaminant=FALSE,
                       master_protein_col)

  }

  obj <- remove_no_master(obj, master_protein_col)

  # remove features with non-unique master proteins
  if (unique_master) {
    # First, add a column with the number of master proteins
    rowData(obj)[["Number.of.Protein.Groups"]] <- sapply(
      strsplit(rowData(obj)[[master_protein_col]], split=';'), length)

    obj <- remove_non_unique_master_protein(obj, master_protein_col)
  }

  # remove features with quantification warnings if necessary
  if (remove_no_quant) {
    obj <- remove_no_quant_assay(obj, master_protein_col)
  }

  return(obj)

}
