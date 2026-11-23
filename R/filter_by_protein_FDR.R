#' Filter to remove peptides from proteins failing the protein FDR threshold
#'
#' @description This function filters the content of a summarised experiment
#' based on a separate file detailing the proteins passing the FDR confidence threshold
#'
#' This function is used to filter peptides to retain only those from proteins
#' passing the FDR threshold
#'
#' Note that peptides whose master protein has no row in
#' `protein_fdr_filename` are also removed, since no FDR confidence can be
#' established for them. With the default `multi_accession = "remove"`, this
#' includes every peptide whose master protein is a group of several accessions,
#' which cannot match a single-accession row in the protein-level output.
#' Applying this function after filtering to unique master proteins (see
#' `filter_features_pd_dda(unique_master = TRUE)`) keeps the two effects
#' separate. Where groups have been kept (`unique_master = FALSE`),
#' `multi_accession = "all"` or `"any"` assesses them on the confidence of their
#' member accessions. Accessions listed in `retain_proteins` are exempt.
#'
#' Groups retained by `multi_accession = "all"` or `"any"` keep their compound
#' accession string; they are not split or reassigned. On aggregation to protein
#' level they become separate compound rows, which can overlap the rows for
#' their member accessions.
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
#' @param multi_accession `string`. How to treat features whose master protein is
#' a group of several accessions, separated by ";". `"remove"` (default) matches
#' the whole string, so a group is removed unless the protein-level output
#' contains that same string. `"all"` retains the feature if every accession has
#' 'High' confidence, with an accession absent from the protein-level output
#' counting as not 'High'. `"any"` retains the feature if at least one accession
#' has 'High' confidence. Features with a single accession are unaffected.
#' @return Returns a `SummarisedExperiment` with the filtered Proteome Discoverer output.
#' @export
filter_by_protein_fdr <- function(obj,
                                  protein_fdr_filename,
                                  protein_col_peptide = "Master.Protein.Accessions",
                                  protein_col_protein = "Accession",
                                  protein_FDR_col = 'Protein.FDR.Confidence.Combined',
                                  retain_proteins = NULL,
                                  multi_accession = c("remove", "all", "any")){

  check_se(obj)

  multi_accession <- match.arg(multi_accession)

  proteinData <- read.delim(file = protein_fdr_filename)

  ## Extract protein accessions and corresponding confidence from search output file
  protein_search_output <- proteinData %>%
    select(!!sym(protein_col_protein), !!sym(protein_FDR_col)) %>%
    distinct()

  ## A blank accession is not an identifier, so may be repeated
  protein_accessions <- as.character(protein_search_output[[protein_col_protein]])
  conflicting <- protein_accessions[duplicated(protein_accessions)] %>%
    unique() %>%
    setdiff(c(NA, ""))

  if (length(conflicting) > 0) {
    stop(sprintf("Accessions with more than one value for %s in protein_fdr_filename: %s",
                 protein_FDR_col, paste(conflicting, collapse = ", ")))
  }

  ## Look up the confidence for each master protein
  confidence <- setNames(protein_search_output[[protein_FDR_col]], protein_accessions)
  accessions <- as.character(rowData(obj)[[protein_col_peptide]])

  rowData(obj)$Protein.Confidence <- unname(confidence[accessions])

  ## Assess groups on the confidence of their member accessions
  if (multi_accession != "remove") {
    is_group <- grepl(";", accessions, fixed = TRUE)
    agg <- if (multi_accession == "all") all else any
    pass <- vapply(strsplit(accessions[is_group], ";\\s*"), function(x)
      agg(confidence[x] %in% "High"), logical(1))
    rowData(obj)$Protein.Confidence[is_group][pass] <- "High"
  }

  obj <- obj[(!is.na(rowData(obj)$Protein.Confidence) | rowData(obj)[[protein_col_peptide]] %in% retain_proteins), ]

  obj <- obj[(rowData(obj)$Protein.Confidence == 'High' | rowData(obj)[[protein_col_peptide]] %in% retain_proteins),]

  message_parse(rowData(obj), protein_col_peptide,
                'Removing peptides from proteins not passing the FDR threshold')

  return(obj)
}
