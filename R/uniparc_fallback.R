get_uniparc_fallback_one <- function(accession){

  entry <- httr::content(
    httr::GET(sprintf("https://rest.uniprot.org/uniprotkb/%s.json", accession)),
    as = "parsed"
  )

  if(!identical(entry$entryType, "Inactive")) return(NULL)

  upi <- entry$extraAttributes$uniParcId
  if(is.null(upi)) return(NULL)

  xrefs <- httr::content(
    httr::GET(sprintf("https://rest.uniprot.org/uniparc/%s.json", upi)),
    as = "parsed"
  )$uniParcCrossReferences

  named <- Filter(function(x) !is.null(x$geneName), xrefs)
  if(length(named)==0) return(NULL)

  best <- named[[which.max(as.Date(sapply(named, `[[`, "lastUpdated")))]]

  data.frame(
    UniprotID = accession,
    Gene.Names = best$geneName,
    Gene.Names.First = best$geneName,
    Protein.names = best$proteinName,
    Annotation.Source = "UniParc (retired UniProtKB entry)",
    stringsAsFactors = FALSE
  )
}

#' Recover annotations for retired UniProtKB accessions via UniParc
#'
#' @description MaxQuant searches are run against a FASTA snapshot of UniProt
#' taken at search time. By the time the annotation step runs, some
#' accessions have since been deleted from UniProtKB (most commonly because
#' the underlying Ensembl gene model was withdrawn or superseded), so a
#' normal `uniprotREST::uniprot_map()` call (as used by
#' \code{\link{get_go_terms}}) returns blank `Gene.Names`/`Protein.names` for
#' them even though the protein may be a real, well-supported hit.
#'
#' For each accession, this function queries the live UniProtKB entry. If it
#' is inactive, it follows the entry's `uniParcId` to UniParc and returns the
#' most recently-updated UniParc cross-reference that has a gene name. This
#' two-step lookup (rather than bulk ID mapping to UniParc) is used because
#' accessions are versioned internally (e.g. `A0A5G2QPJ4.1` vs `.2`), and bulk
#' mapping to UniParc can return multiple UniParc entries for one accession
#' with no clean way to disambiguate; resolving via the entry's own
#' `uniParcId` is unambiguous.
#'
#' Intended to be called only with the subset of accessions that came back
#' empty/`"deleted"` from a prior `uniprot_map()` call, not a full accession
#' list.
#'
#' @param accessions `character vector` UniProt accessions to resolve via UniParc.
#' @return `data.frame` with columns `UniprotID`, `Gene.Names`,
#' `Gene.Names.First`, `Protein.names`, `Annotation.Source`. Accessions that
#' are not inactive, have no `uniParcId`, or have no UniParc cross-reference
#' with a gene name are omitted from the result.
#' @export
#' @examples
#' \dontrun{
#' get_uniparc_fallback("A0A5G2QPJ4")
#' }
get_uniparc_fallback <- function(accessions){

  results <- Filter(Negate(is.null), lapply(accessions, get_uniparc_fallback_one))

  if(length(results)==0){
    return(data.frame(
      UniprotID = character(), Gene.Names = character(),
      Gene.Names.First = character(), Protein.names = character(),
      Annotation.Source = character(), stringsAsFactors = FALSE
    ))
  }

  do.call(rbind, results)
}

#' Get UniProt annotation details, with UniParc fallback for retired accessions
#'
#' @description Queries `uniprotREST::uniprot_map()` for gene/protein name
#' annotations for a vector of UniProt accessions, following the same
#' `uniprot_map()` usage as \code{\link{get_go_terms}}. Because MaxQuant
#' searches are run against a UniProt snapshot that may be weeks to years out
#' of date by the time annotation runs, some accessions may since have been
#' deleted from UniProtKB; these come back from `uniprot_map()` with blank
#' `Gene.Names` and `Protein.names == "deleted"`. For any such accession,
#' \code{\link{get_uniparc_fallback}} is used to recover the last-known
#' annotation from UniParc.
#'
#' @param accessions `character vector` UniProt accessions.
#' @param verbosity `integer` Verbosity level for uniprotREST::uniprot_map.
#' @return `data.frame` with one row per accession and columns `UniprotID`,
#' `Entry.Name`, `Reviewed`, `Protein.names`, `Gene.Names`, `Organism`,
#' `Length`, `Gene.Names.First`, `Annotation.Source`. `Annotation.Source` is
#' `"UniProtKB (live)"` for accessions resolved directly from UniProtKB, or
#' `"UniParc (retired UniProtKB entry)"` for accessions recovered via
#' \code{\link{get_uniparc_fallback}}. Accessions that are deleted and have no
#' recoverable UniParc annotation keep their original (blank/`"deleted"`)
#' values from `uniprot_map()`.
#' @export
#' @examples
#' \dontrun{
#' get_uniprot_details(c('O76024', 'Q03135', 'A0A5G2QPJ4'))
#' }
get_uniprot_details <- function(accessions, verbosity=0){

  uniprot2details <- uniprotREST::uniprot_map(
    accessions,
    from='UniProtKB_AC-ID',
    method='stream',
    verbosity=verbosity) %>%
    dplyr::rename('UniprotID'=From) %>%
    select(-Entry) %>%
    mutate(Gene.Names.First=gsub(' .*','', Gene.Names),
           Annotation.Source='UniProtKB (live)')

  needs_fallback <- uniprot2details$UniprotID[
    uniprot2details$Gene.Names=='' | uniprot2details$Protein.names=='deleted'
  ]

  if(length(needs_fallback)==0) return(uniprot2details)

  fallback_result <- get_uniparc_fallback(needs_fallback)

  if(nrow(fallback_result)==0) return(uniprot2details)

  dplyr::rows_update(uniprot2details, fallback_result, by='UniprotID')
}
