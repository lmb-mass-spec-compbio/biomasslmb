#' Update PSM level protein/master protein assignments using peptide-level assignments
#'
#' @description
#' When multiple search engines are used in Proteome Discoverer for TMT proteomics
#' the PSM to protein and master protein
#' assignments may not be consistent between the search engines. We can resolve this
#' using the peptide level output from Proteome Discoverer, which has consistent
#' assignments
#'
#' @details **Leucine/Isoleucine**
#' Leucine/Isoleucine are isobaric and only generate distinct fragments with higher
#' collision energies which are not typically used in proteomics analysis. Thus,
#' L and I are usually indistinguishable.
#'
#' The spectrum search engines deal with this differently in terms of how they report the
#' multiple possible sequence matches and their assignments.
#'
#' The peptide level Proteome Discoverer output selects one possible sequence and
#' uses the union of reported assignments from the PSM level for the peptide level assignments.
#' Thus, when we use the peptide level assignments, we will lose PSM level data
#' for the redundant sequences which were not selected. This is not an issue but will result
#' in the function output being shorter than the input.
#'
#'
#' @param obj `SummarisedExperiment` containing PSM-level output from Proteome Discoverer.
#' @param pep_inf `string`. Filepath for Peptide-level output
#' @param verbose `boolean`. Default is FALSE; Don't output tallies of PSMs/proteins.
#' @param master_protein_col `string`. Name of column containing master
#' proteins. (Only used when verbose=TRUE)
#' @return Returns a `SummarisedExperiment` with the standardised assignments
#' @export
update_peptide_assignments <- function(obj,
                                       pep_inf,
                                       verbose=FALSE,
                                       master_protein_col = "Master.Protein.Accessions"){

  # check arguments
  check_se(obj)

  if(verbose){
    # print input summary
    message_parse(rowData(obj), master_protein_col, "Input")
  }

  pep_infdf <- read.delim(pep_inf)

  assignment_cols <- c('Protein.Accessions', 'Master.Protein.Accessions',
                       'Positions.in.Master.Proteins', 'Master.Protein.Descriptions')
  assignment_cols_present <- intersect(assignment_cols, colnames(pep_infdf))

  pep2assignments <- pep_infdf %>% select(Sequence, all_of(assignment_cols_present)) %>%
    distinct()

  if('Protein.Accessions' %in% assignment_cols_present){
    # Protein.Accessions column is not always consistent, eg. Sequest may include contaminant IDs, where
    # Comet does not. We will use whichever value is longest
    pep2assignments <-  pep2assignments %>%
      mutate(prot_acc_len=nchar(Protein.Accessions)) %>%
      group_by(across(all_of(setdiff(assignment_cols_present, 'Protein.Accessions')))) %>%
      slice_max(prot_acc_len, n=1) %>%
      select(-prot_acc_len) %>% ungroup()
  }

  if(sum(duplicated(pep2assignments$Sequence))>0){
    duplicated_pep <- pep2assignments$Sequence[duplicated(pep2assignments$Sequence)]
    print(filter(pep2assignments, Sequence==duplicated_pep[[1]]))
    stop(sprintf(
      'Peptide to protein/master protein assignments in
      %s
      are not unique. See example', pep_inf))

  }

  assignment_cols_to_remove <- intersect(assignment_cols_present, colnames(rowData(obj)))

  updated_rowData <- rowData(obj) %>%
    data.frame() %>%
    tibble::rownames_to_column('row_value') %>%
    select(-all_of(assignment_cols_to_remove)) %>%
    merge(pep2assignments, by='Sequence') %>%
    tibble::column_to_rownames('row_value')


  new_obj <- obj[rownames(obj) %in% rownames(updated_rowData),]
  rowData(new_obj) <- updated_rowData[rownames(new_obj),]

  if(verbose) message_parse(rowData(new_obj),
                            master_protein_col,
                            "Updating peptide assignments")

  return(new_obj)

}

#' Remove redundant Peptide Spectrum Matches
#'
#' @description
#' When multiple search engines are used in Proteome Discoverer for TMT proteomics,
#' the PSM-level quantification contains duplicate rows for each search engine match
#' to the same spectrum. We can use the scan number columns
#' and quantitative data to identify the duplicate rows and exclude them.
#'
#' The removal of duplicate quantification should be performed following initial
#' filtering to 1) remove rejected PSMs (see filter_TMT_PSMs) and 2)
#' retain only rank 1 matches (Filter on Rank and Search.Engine.Rank columns)
#'
#' @param obj `SummarisedExperiment` containing PSM-level output from Proteome Discoverer.
#' @param verbose `boolean`. Don't output tallies of PSMs/proteins; Default=FALSE
#' @param master_protein_col `string`. Name of column containing master
#' @param validate `logical`. Check that sequence to protein assignments are consistent. Default=TRUE
#' proteins. (Only used when verbose=TRUE)
#' @return Returns a `SummarisedExperiment` with the redundant quantification removed
#' @export
remove_redundant_psm_quant <- function(obj,
                                       verbose=FALSE,
                                       master_protein_col = "Master.Protein.Accessions",
                                       validate=TRUE){

  # check arguments
  check_se(obj)

  spectra_id_cols <- c("First.Scan", "Last.Scan", "Master.Scans", "File.ID")

  if(length(setdiff(spectra_id_cols, colnames(rowData(obj))))!=0){
    stop(sprintf(
      'Following columns must be present in rowData(obj): %s',
      paste(spectra_id_cols, collapse=', ')))
  }

  if(validate){

    validation_columns <- c('Sequence', 'Protein.Accessions')

    if(length(setdiff(validation_columns, colnames(rowData(obj))))!=0){
      stop(sprintf(
        'Following columns must be present in rowData(obj): %s',
        paste(validation_columns, collapse=', ')))
    }

    seq2acc <- rowData(obj)[,validation_columns] %>%
      data.frame() %>%
      distinct()

    if(sum(duplicated(seq2acc$Sequence))>0){
      stop('Multiple Protein.Accessions values observed for the same Sequence.
           This will occur if you have multiple search engines and can be rectified
           with update_peptide_assignments()')
    }
  }

  if(verbose){
    # print input summary
    message_parse(rowData(obj), master_protein_col, "Input")
  }



  quant_data <- cbind(assay(obj), rowData(obj)[,spectra_id_cols]) %>%
    data.frame()

  duplicated_quant <- duplicated(quant_data)

  new_obj <- obj[!duplicated_quant,]

  if(verbose) message_parse(rowData(new_obj),
                            master_protein_col,
                            "Excluding duplicated quantification")

  return(new_obj)
}
