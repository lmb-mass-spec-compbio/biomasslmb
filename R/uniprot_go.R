#' Obtain GO term annotations for proteins
#'
#' @description Given a set of UniProt IDs, this function queries UniProt to obtain
#' the annotated GO terms. Optionally, these GO terms can be expanded to include all ancestors too, which
#' can be helpful when using GO over-representation/enrichment tools that do not
#' consider the GO term heirarchy. Note that expansion will significantly increase the run-time
#'
#' @param uniprotIDs `character vector` Uniprot IDs
#' @param expand_terms `logical` Should GO terms be expanded to include all ancestors
#' @param verbosity `integer` Verbosity level for uniprotREST::uniprot_map
#' @return `data.frame` object.
#' @export
#' @importFrom uniprotREST uniprot_map
#' @examples
#' uniprotIDs <- c('O76024', 'Q03135', 'Q96T23')
#' get_go_terms(uniprotIDs, expand_terms=TRUE)
get_go_terms <- function(UniprotID, expand_terms=FALSE, verbosity=0){

  go_res <- uniprot_map(
    ids = UniprotID,
    from = "UniProtKB_AC-ID",
    verbosity=verbosity,
    method='stream',
    to = "UniProtKB",
    fields = "go",
  ) %>%
    dplyr::rename(UNIPROTKB=From)

  # The Uniprot API returns a single row for each protein, with all GO terms pasted together
  # Below, we separate each GO term into a separate row
  go_res <- go_res %>%
    filter(Gene.Ontology..GO.!='') %>% separate_rows(Gene.Ontology..GO., sep = '; ') %>%
    tidyr::separate_wider_delim(Gene.Ontology..GO., ' [GO:', names=c('GO_desc', 'GO.ID')) %>%
    mutate(GO.ID=paste0('GO:', gsub(']$', '', GO.ID)))

  if(expand_terms){
    go_res <- get_ancestor_go(go_res)
  }

  return(go_res)

}


