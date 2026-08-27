# Select the most recently-updated UniParc cross-reference that has both a
# gene name and a protein name (some cross-references, e.g. Ensembl, carry a
# gene name but no protein name; Protein.names is a required output column).
pick_best_uniparc_xref <- function(xrefs){
  named <- Filter(function(x) !is.null(x$geneName) && !is.null(x$proteinName), xrefs)
  if(length(named)==0) return(NULL)
  named[[which.max(as.Date(sapply(named, `[[`, "lastUpdated")))]]
}

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

  best <- pick_best_uniparc_xref(xrefs)
  if(is.null(best)) return(NULL)

  data.frame(
    UniprotID = accession,
    Gene.Names = best$geneName,
    Gene.Names.First = best$geneName,
    Protein.names = best$proteinName,
    Annotation.Source = "UniParc (retired UniProtKB entry)",
    stringsAsFactors = FALSE
  )
}

# GET a batch of UniProt JSON entries in parallel via httr2. UniProt doesn't
# publish a rate limit for these single-entry read endpoints; 10 req/s is
# comfortably above what get_uniparc_fallback_one()'s one-at-a-time httr::GET
# calls could ever sustain, while still leaving room below any undocumented
# server-side limit. on_error = "continue" means a single failed/erroring
# fetch (e.g. a 404) drops out as NULL rather than aborting the whole batch,
# matching get_uniparc_fallback_one()'s existing "unresolvable accessions are
# just omitted" behaviour.
fetch_uniprot_json_parallel <- function(urls){
  if(length(urls)==0) return(list())

  reqs <- lapply(urls, uniprotREST::uniprot_request, rate = 10)
  resps <- httr2::req_perform_parallel(reqs, on_error = "continue")

  lapply(resps, function(resp){
    if(!inherits(resp, "httr2_response")) return(NULL)
    httr2::resp_body_json(resp)
  })
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
#' most recently-updated UniParc cross-reference that has both a gene name
#' and a protein name (some cross-references, e.g. Ensembl, carry a gene name
#' but no protein name). This two-step lookup (rather than bulk ID mapping to
#' UniParc) is used because
#' accessions are versioned internally (e.g. `A0A5G2QPJ4.1` vs `.2`), and bulk
#' mapping to UniParc can return multiple UniParc entries for one accession
#' with no clean way to disambiguate; resolving via the entry's own
#' `uniParcId` is unambiguous.
#'
#' Intended to be called only with the subset of accessions that came back
#' empty/`"deleted"` from a prior `uniprot_map()` call, not a full accession
#' list.
#'
#' The two lookups this requires (UniProtKB entry, then UniParc entry) are
#' each done for all `accessions` in parallel via `httr2::req_perform_parallel()`,
#' rather than one accession at a time, since this is the dominant runtime
#' cost of \code{\link{get_uniprot_details}} for inputs with many retired
#' accessions.
#'
#' @param accessions `character vector` UniProt accessions to resolve via UniParc.
#' @return `data.frame` with columns `UniprotID`, `Gene.Names`,
#' `Gene.Names.First`, `Protein.names`, `Annotation.Source`. Accessions that
#' are not inactive, have no `uniParcId`, or have no UniParc cross-reference
#' with both a gene name and a protein name are omitted from the result.
#' @export
#' @examples
#' \dontrun{
#' get_uniparc_fallback("A0A5G2QPJ4")
#' }
get_uniparc_fallback <- function(accessions){

  empty_result <- data.frame(
    UniprotID = character(), Gene.Names = character(),
    Gene.Names.First = character(), Protein.names = character(),
    Annotation.Source = character(), stringsAsFactors = FALSE
  )

  if(length(accessions)==0) return(empty_result)

  kb_entries <- fetch_uniprot_json_parallel(
    sprintf("https://rest.uniprot.org/uniprotkb/%s.json", accessions)
  )

  upis <- vapply(kb_entries, function(entry){
    if(is.null(entry) || !identical(entry$entryType, "Inactive")) return(NA_character_)
    upi <- entry$extraAttributes$uniParcId
    if(is.null(upi)) NA_character_ else upi
  }, character(1))

  needs_uniparc <- !is.na(upis)
  if(!any(needs_uniparc)) return(empty_result)

  sel_accessions <- accessions[needs_uniparc]
  uc_entries <- fetch_uniprot_json_parallel(
    sprintf("https://rest.uniprot.org/uniparc/%s.json", upis[needs_uniparc])
  )

  results <- Map(function(accession, entry){
    if(is.null(entry)) return(NULL)
    best <- pick_best_uniparc_xref(entry$uniParcCrossReferences)
    if(is.null(best)) return(NULL)
    data.frame(
      UniprotID = accession,
      Gene.Names = best$geneName,
      Gene.Names.First = best$geneName,
      Protein.names = best$proteinName,
      Annotation.Source = "UniParc (retired UniProtKB entry)",
      stringsAsFactors = FALSE
    )
  }, sel_accessions, uc_entries)

  results <- Filter(Negate(is.null), results)
  if(length(results)==0) return(empty_result)

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
#' Some accessions come back from `uniprot_map()` as multiple rows (they have
#' subsequently been 'demerged'); these are collapsed to one row per
#' accession, with `;`-separated values where the demerged rows differed.
#' @param check_complete `string`, one of `"error"` (default), `"warn"`, or
#' `"none"`. Passed to `uniprotREST::uniprot_map()`, which compares the
#' number of submitted accessions against the number mapped plus the number
#' UniProt reports as failed, to catch ID mapping jobs that silently return a
#' truncated result. See `?uniprotREST::uniprot_map` for details.
#'
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
get_uniprot_details <- function(accessions, verbosity=0, check_complete='error'){

  uniprot2details_raw <- uniprotREST::uniprot_map(
    accessions,
    from='UniProtKB_AC-ID',
    method='stream',
    verbosity=verbosity,
    check_complete=check_complete) %>%
    dplyr::rename('UniprotID'=From) %>%
    select(-Entry) %>%
    mutate(Gene.Names.First=gsub(' .*','', Gene.Names))

  # detect accessions needing fallback on the raw, one-row-per-mapping data,
  # since collapsing (below) joins blank/'deleted' values with ';' and would
  # break this exact-match check for accessions with multiple rows
  needs_fallback <- unique(uniprot2details_raw$UniprotID[
    uniprot2details_raw$Gene.Names=='' | uniprot2details_raw$Protein.names=='deleted'
  ])

  uniprot2details <- uniprot2details_raw %>%
    # some UniprotIDs generate multiple rows, likely because they have
    # subsequently been 'demerged' -- merge them
    group_by(UniprotID) %>%
    summarise_all(.funs=function(x) paste(x, collapse=';')) %>%
    mutate(Annotation.Source='UniProtKB (live)')

  if(length(needs_fallback)==0) return(uniprot2details)

  fallback_result <- get_uniparc_fallback(needs_fallback)

  if(nrow(fallback_result)==0) return(uniprot2details)

  dplyr::rows_update(uniprot2details, fallback_result, by='UniprotID')
}

#' Map per-accession UniProt details onto protein IDs with multiple accessions
#'
#' @description Takes the `data.frame` returned by
#' \code{\link{get_uniprot_details}} (one row per single accession) and a
#' vector of protein IDs where an element may itself be several accessions
#' joined by `sep` -- the `Master.Protein.Accessions` column from Proteome
#' Discoverer `PeptideGroups.txt` output is the common source of this. Each
#' `protein_ids` element is split into its constituent accessions, the
#' matching rows of `uniprot2details` are looked up, then re-collapsed back
#' to one row per original (possibly multi-accession) `protein_ids` element,
#' with `sep`-separated values where the constituent accessions' details
#' differ. Does not call \code{\link{get_uniprot_details}} itself, so the
#' same `uniprot2details` result can also be used on its own, without
#' fetching annotations twice.
#'
#' @param uniprot2details `data.frame`, the output of
#' \code{\link{get_uniprot_details}}, covering (at least) every accession
#' that appears in `protein_ids` once split by `sep`.
#' @param protein_ids `character vector` protein IDs, each optionally
#' containing several UniProt accessions joined by `sep`.
#' @param sep `string` delimiter separating accessions within a
#' `protein_ids` element, and also used to join differing values in the
#' collapsed output. Default `"; "`, matching Proteome Discoverer's
#' `Master.Protein.Accessions` column.
#'
#' @return `data.frame` with one row per unique `protein_ids` element and
#' the same columns as `uniprot2details`. Where the constituent accessions
#' agree on a value it appears once; where they differ the distinct values
#' are joined by `sep`.
#' @export
#' @examples
#' \dontrun{
#' uniprot2details <- get_uniprot_details(c('O76024', 'Q03135', 'A0A5G2QPJ4'))
#' collapse_uniprot_details_multi_accession(
#'   uniprot2details, c('O76024', 'Q03135; A0A5G2QPJ4'))
#' }
collapse_uniprot_details_multi_accession <- function(uniprot2details, protein_ids, sep='; '){

  protein_ids <- unique(protein_ids)

  data.frame(
    UniprotID = protein_ids,
    UniprotID_single = protein_ids,
    stringsAsFactors = FALSE) %>%
    separate_rows(UniprotID_single, sep=sep) %>%
    merge(uniprot2details, by.x='UniprotID_single', by.y='UniprotID') %>%
    select(-UniprotID_single) %>%
    group_by(UniprotID) %>%
    summarise_all(.funs=function(x) paste(unique(x), collapse=sep))
}
