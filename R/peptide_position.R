#' Add a column describing the position of the peptide sequence with respect
#' to the protein
#'
#' @description Identify the position of the peptide sequence in the protein.
#'
#' The peptide position is undefined (NA) if:
#' 1. The peptide sequence has multiple master proteins
#' 2. The protein does not exist in the proteome fasta
#' 3. The peptide is repeated in the protein sequence.
#'
#' The peptide_position_info column details why the peptide positions columns are NA
#'
#' @param obj `SummarizedExperiment` with PD output at PSM/peptide level
#' @param proteome_fasta `character` Filepath for proteome fasta
#' @param master_protein_col `character` Column name for master protein
#' @param protein_col_split `character` Delimiter for multiple proteins in master_protein_col
#' @param sequence_col `character` Column name for peptide sequence
#'
#' @return `SummarizedExperiment` with extra row columns detailing the peptide position
#' @export
add_peptide_positions <- function(obj,
                                  proteome_fasta,
                                  master_protein_col = "Master.Protein.Accessions",
                                  protein_col_split='; ',
                                  sequence_col = 'Sequence') {
  warning('This function will provide the positions of peptides within proteins,
          irrespective of whether the peptide could have originated from those locations
          given the digestion enzyme employed. If you want to use in-silico digestion
          to filter the peptide positions, use the add_peptide_positions_from_cleavage
          function instead.')
  check_se(obj)

  proteome <- readAAStringSet(proteome_fasta)
  names(proteome) <- sapply(base::strsplit(names(proteome), split = "\\|"), "[[", 2)

  combine_protein_peptide_positions <- function(proteome, protein, sequence) {

    # Given a master protein(s) and AA sequence,
    # return the AA position with respect to protein sequence
    # If more than one position, possible, return NA

    if (length(strsplit(protein, 'protein_col_split')[[1]]) > 1) {
      return(c(NA, NA, NA, 'Multiple proteins'))
    }

    if (!protein %in% names(proteome)) {
      return(c(NA, NA, NA, 'Protein not in fasta'))
    }

    m <- gregexpr(sprintf('(?=(%s))', sequence), proteome[[protein]], perl=T)[[1]]
    peptide_start <- attr(m,"capture.start")
    peptide_end <- attr(m,"capture.start") + attr(m,"capture.length") - 1

    if (length(peptide_start) != 1) {
      return(c(NA, NA, NA, 'Multiple peptide locations'))
    }

    else {
      return(c(peptide_start, peptide_end, length(proteome[[protein]]), ''))
    }
  }

  rowData(obj)[, c('peptide_start', 'peptide_end', 'protein_length', 'peptide_position_info')] <- t(apply(
    rowData(obj),
    MARGIN = 1, function(x) combine_protein_peptide_positions(
      proteome, x[[master_protein_col]], x[[sequence_col]]
    )
  ))

  rowData(obj)$peptide_start <- as.numeric(rowData(obj)$peptide_start)
  rowData(obj)$peptide_end <- as.numeric(rowData(obj)$peptide_end)
  rowData(obj)$protein_length <- as.numeric(rowData(obj)$protein_length)

  return(obj)
}
