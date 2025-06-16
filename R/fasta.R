#' Extract the contaminants protein accessions from a cRAP fasta file
#'
#' @param cont_fasta_inf `character` Filepath to the contaminants fasta
#' @return Returns a `list` of contaminant IDs
#' @export
get_crap_fasta_accessions <- function(cont_fasta_inf){
  # Load the cRAP FASTA used for the PD search
  crap.fasta <- Biostrings::fasta.index(cont_fasta_inf, seqtype = "AA")

  # Define a base R version of stringr::str_extract_all()
  # base R str_extract
  str_extract_all <- function(pattern, string) {
    gregexpr(pattern, string, perl = TRUE) %>%
      regmatches(string, .) %>%
      unlist()
  }

  # Extract the non cRAP UniProt accessions associated with each cRAP protein
  crap.accessions <- crap.fasta %>%
    pull(desc) %>%
    str_extract_all("(?<=\\|).*?(?=\\|)", .) %>%
    unlist()

  return(crap.accessions)
}

#' Extract the contaminants protein accessions from a MaxQuant contaminants fasta file
#'
#' @param cont_fasta_inf `character` Filepath to the contaminants fasta. Defaults to NULL
#' in which case, contaminants.fasta.gz from this package will be used
#' @return Returns a `list` of contaminant IDs
#' @export
get_maxquant_cont_accessions <- function(cont_fasta_inf=NULL){
  if(is.null(cont_fasta_inf)){

    message('No contaminants fasta supplied, defaulting to contaminants.fasta.gz from biomasslmb')

    cont_fasta_inf <- system.file(
      "extdata", "contaminants.fasta.gz",
      package = "biomasslmb"
    )
  }
  cont.fasta <- Biostrings::fasta.index(cont_fasta_inf, seqtype = "AA")

  cont.accessions <- cont.fasta %>%
    pull(desc) %>%
    sapply(function(x) gsub(' .*', '', x)) %>%
    unname()

  cont.accessions <- paste0('CON__', cont.accessions)

  return(cont.accessions)
}
