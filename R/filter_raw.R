#' Filter Proteome Discoverer DDA output
#'
#' @description This function filters the output .txt files (peptide groups or PSMs) from
#' Proteome Discoverer for DDA, based on various criteria:
#'
#' 1. Remove features without a master protein
#' 2. (Optional) Remove features without a unique master protein  (i.e.
#'    Number.of.Protein.Groups == 1)
#' 3. (Optional) Remove features matching a cRAP protein
#' 4. (Optional) Remove features matching any protein associated with
#'   a cRAP protein (see below)
#' 5. Remove features without quantification values (only if TMT or SILAC
#'   are `TRUE` and `level = "peptide"`.)
#'
#' @details **Associated cRAP proteins** are proteins which have at least
#' one feature shared with a cRAP protein. It has been observed that the cRAP
#' fasta files often do not contain all possible cRAP proteins e.g. some features
#' can be assigned to a keratin which is not in the provided cRAP database.
#'
#' In the example below, using `filter_associated_crap = TRUE` will filter out f2 and f3 in
#' addition to f1, regardless of the value in the Master.Protein.Accession column.
#'
#' ```
#' feature  Protein.Accessions         Master.Protein.Accessions
#' f1       protein1, protein2, cRAP,  protein1,
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
#' @param silac `logical`. Is the experiment a SILAC experiment?
#' @param TMT `logical`. Is the experiment a TMT experiment?
#' @param level `string`. Type of input file, must be one of either
#' `"peptide"` or `"PSM"`.
#' @param filter_crap `logical`. Filter out features which match a cRAP
#' protein.
#' @param crap_proteins `character vector`. The protein IDs form the cRAP proteins,
#' @param filter_associated_crap `logical`. Filter out features which
#' match a cRAP associated protein.
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
#' # extract the UniProt accessions from the cRAP FASTA headers
#' crap_accessions <- get_crap_fasta_accessions(crap_fasta_inf)
#'
#' # filter thg peptides
#' psm2 <- filter_features_pd_dda(
#'   obj = tmt_qf[['psms_raw']],
#'   master_protein_col = "Master.Protein.Accessions",
#'   protein_col = "Protein.Accessions",
#'   unique_master = TRUE,
#'   TMT = TRUE,
#'   level = "PSM",
#'   filter_crap = TRUE,
#'   crap_proteins = crap_accessions,
#'   filter_associated_crap = TRUE
#' )
#'
#'
#' }
#' @export
filter_features_pd_dda <- function(obj,
                                   master_protein_col = "Master.Protein.Accessions",
                                   protein_col = "Protein.Accessions",
                                   unique_master = TRUE,
                                   silac = FALSE,
                                   TMT = FALSE,
                                   level = "peptide",
                                   filter_crap = TRUE,
                                   crap_proteins = NULL,
                                   filter_associated_crap = TRUE,
                                   remove_no_quant = TRUE) {


  # check arguments
  check_se(obj)

  stopifnot(level %in% c("PSM", "peptide"))
  if (filter_crap) {
    if (is.null(crap_proteins)) {
      stop("must supply the crap_proteins argument to filter cRAP proteins")
    }
  }

  # print input summary
  message("Filtering data...")
  message_parse(rowData(obj), master_protein_col, "Input")

  # remove crap proteins
  if (filter_crap) {
    message(sprintf(
      "%s cRAP proteins supplied",
      length(crap_proteins)
    ))

    # Identify features for crap_proteins
    crap_features <- rowData(obj)[[master_protein_col]] %in% crap_proteins |
      grepl("cRAP", rowData(obj)[[protein_col]])


    # identify associated crap proteins first if necessary
    if (filter_associated_crap) {
      associated_crap <- rowData(obj[crap_features,])[[protein_col]] %>%
        strsplit("; ") %>%
        unlist()

      associated_crap <- associated_crap[!grepl("cRAP", associated_crap)]

      message(sprintf(
        "%s proteins identified as 'cRAP associated'",
        length(associated_crap)
      ))
    }

    # then remove normal crap proteins
    obj <- obj[!crap_features,]
    message_parse(rowData(obj), master_protein_col, "cRAP features removed")

    # then remove associated crap proteins if necessary
    if (filter_associated_crap) {
      if (length(associated_crap) > 0) {
        # remove isoforms
        associated_crap_no_isoform <- unique(sapply(strsplit(associated_crap, "-"), "[[", 1))
        associated_crap_regex <- paste(associated_crap_no_isoform, collapse = "|")
        obj <- obj[!grepl(associated_crap_regex, rowData(obj)[[protein_col]]), ]

        message_parse(rowData(obj), master_protein_col, "associated cRAP features removed")
      }
    }
    if('Contaminant' %in% colnames(rowData(obj))){
      obj <- obj[rowData(obj)$Contaminant=='False',]
      message_parse(rowData(obj), master_protein_col, "PD-labelled 'Contaminants' removed")
    }
  }

  # remove features without a master protein
  missing_master_protein <- (is.na(rowData(obj)[[master_protein_col]]) |
                               rowData(obj)[[master_protein_col]] == '')

  obj <- obj[!missing_master_protein,]
  message_parse(rowData(obj), master_protein_col, "features without a master protein removed")

  # remove features with non-unique master proteins
  if (unique_master) {
    # Add option(s) here to handle cases where column does not exist
    obj <- obj[rowData(obj)[["Number.of.Protein.Groups"]] == 1, ]
    message_parse(rowData(obj), master_protein_col, "features with non-unique master proteins removed")
  }

  # remove features with quantification warnings if necessary
  if (remove_no_quant) {
    obj <- obj[rowSums(is.finite(assay(obj)))>0,]
     message_parse(rowData(obj), master_protein_col, "features without quantification removed")
  }

  return(obj)

}


