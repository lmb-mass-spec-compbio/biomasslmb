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


#' Extract Peptide Position Start and End Columns from a QFeatures Object
#'
#' Parses the \code{Positions.in.Master.Proteins} column from the \code{rowData}
#' of a specified \code{SummarizedExperiment} within a \code{QFeatures} object,
#' and extracts the start and end positions of each peptide into two new
#' \code{rowData} columns. When multiple entries are present
#' (semicolon-delimited), the function returns semicolon-delimited lists of
#' start and end values.
#'
#' @param qf A \code{QFeatures} object.
#' @param i A single integer or character string specifying the index or name of
#'   the \code{SummarizedExperiment} within \code{qf} to operate on.
#' @param start_col A single character string giving the name of the new column
#'   to store peptide start positions. Defaults to \code{"start"}.
#' @param end_col A single character string giving the name of the new column
#'   to store peptide end positions. Defaults to \code{"end"}.
#'
#' @return The input \code{QFeatures} object with two additional columns in the
#'   \code{rowData} of the specified experiment:
#'   \describe{
#'     \item{start_col}{Character. Semicolon-delimited start position(s) for the peptide.}
#'     \item{end_col}{Character. Semicolon-delimited end position(s) for the peptide.}
#'   }
#'
#' @examples
#' \dontrun{
#' library(QFeatures)
#' # Using default column names
#' qf <- extract_peptide_positions(qf, i = "peptides")
#'
#' # Using custom column names
#' qf <- extract_peptide_positions(qf, i = 1L,
#'                                 start_col = "pep_start",
#'                                 end_col   = "pep_end")
#' }
#'
#' @importFrom QFeatures QFeatures
#' @importFrom SummarizedExperiment rowData rowData<-
#' @importFrom S4Vectors DataFrame
#' @export
extract_peptide_positions <- function(qf, i, start_col = "start", end_col = "end") {

  check_q(qf)

  if (missing(i)) {
    stop("'i' must be specified as an integer index or experiment name.")
  }

  if (length(i) != 1L) {
    stop("'i' must be a single integer or character value.")
  }

  check_se_exists(qf, i)

  if (!is.character(start_col) || length(start_col) != 1L || !nzchar(start_col)) {
    stop("'start_col' must be a single non-empty character string.")
  }

  if (!is.character(end_col) || length(end_col) != 1L || !nzchar(end_col)) {
    stop("'end_col' must be a single non-empty character string.")
  }

  if (start_col == end_col) {
    stop("'start_col' and 'end_col' must be different.")
  }

  # Extract the SummarizedExperiment at index i
  se <- qf[[i]]
  rd <- rowData(se)

  if (!"Positions.in.Master.Proteins" %in% colnames(rd)) {
    stop(
      "'Positions.in.Master.Proteins' column not found in rowData of experiment '",
      i, "'."
    )
  }

  raw <- as.character(rd[["Positions.in.Master.Proteins"]])

  # Parse a single entry string, returning a list with start and end strings
  parse_positions <- function(entry) {
    if (is.na(entry) || !nzchar(trimws(entry))) {
      return(list(start = NA_character_, end = NA_character_))
    }

    matches        <- gregexpr("\\[(\\d+)-(\\d+)\\]", entry, perl = TRUE)
    matched_strings <- regmatches(entry, matches)[[1]]

    if (length(matched_strings) == 0) {
      return(list(start = NA_character_, end = NA_character_))
    }

    parsed <- regmatches(
      matched_strings,
      regexpr("(\\d+)-(\\d+)", matched_strings, perl = TRUE)
    )
    parts  <- strsplit(parsed, "-")
    starts <- sapply(parts, `[`, 1)
    ends   <- sapply(parts, `[`, 2)

    list(
      start = paste(starts, collapse = ";"),
      end   = paste(ends,   collapse = ";")
    )
  }

  results           <- lapply(raw, parse_positions)
  rd[[start_col]]   <- sapply(results, `[[`, "start")
  rd[[end_col]]     <- sapply(results, `[[`, "end")

  rowData(se) <- rd
  qf[[i]]     <- se

  qf
}

#' Evaluate peptide termini against theoretical cleavage positions
#'
#' @param obj SummarizedExperiment
#' @param proteome_fasta character. Path to FASTA file
#' @param digest_enzyme character. Enzyme for cleaver::cleavageRanges
#' @param missed_cleavages numeric vector. Allowed missed cleavages
#' @param master_protein_col character. Column in rowData containing protein IDs
#' @param start_col character. Column containing peptide start positions
#' @param end_col character. Column containing peptide end positions
#'
#' @return SummarizedExperiment with additional rowData columns:
#' N_tryp and C_tryp
#' @export
add_tryptic_termini_flags <- function(obj,
                                      proteome_fasta,
                                      digest_enzyme = "trypsin-simple",
                                      missed_cleavages = c(0, 1, 2),
                                      master_protein_col = "Leading.razor.protein",
                                      start_col = "start",
                                      end_col = "end") {

  check_se(obj)

  proteins <- Biostrings::readAAStringSet(proteome_fasta)

  # Standardise FASTA names to UniProt IDs
  names(proteins) <- gsub('(sp|tr)\\|(\\S*)\\|.*', '\\2', names(proteins))

  pep_ranges <- cleaver::cleavageRanges(
    proteins,
    enzym = digest_enzyme,
    missedCleavages = missed_cleavages
  )

  # Build lookup of valid cleavage coordinates
  cleavage_lookup <- lapply(names(pep_ranges), function(p) {

    df <- as.data.frame(pep_ranges[[p]])

    data.frame(
      protein = p,
      start = df$start,
      end = df$end,
      stringsAsFactors = FALSE
    )
  }) |> dplyr::bind_rows()

  valid_positions <- cleavage_lookup |>
    dplyr::group_by(protein) |>
    dplyr::summarise(
      valid_starts = list(unique(start)),
      valid_ends   = list(unique(end)),
      .groups = "drop"
    )

  rd <- as.data.frame(SummarizedExperiment::rowData(obj))

  merged <- rd |>
    dplyr::left_join(
      valid_positions,
      by = setNames("protein", master_protein_col)
    )

  N_tryp <- mapply(function(s, valid_s) {
    if (is.null(valid_s) || is.na(s)) return(NA)
    s %in% valid_s
  }, merged[[start_col]], merged$valid_starts)

  C_tryp <- mapply(function(e, valid_e) {
    if (is.null(valid_e) || is.na(e)) return(NA)
    e %in% valid_e
  }, merged[[end_col]], merged$valid_ends)

  SummarizedExperiment::rowData(obj)$N_tryp <- as.logical(N_tryp)
  SummarizedExperiment::rowData(obj)$C_tryp <- as.logical(C_tryp)

  return(obj)
}
