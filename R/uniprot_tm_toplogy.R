#' Query UniProt to determine protein transmembrane-domains and topology
#'
#' @description Given a set of UniProt IDs, this function returns the information regarding transmembrane domains and topology.
#'
#' @param uniprotIDs `character vector` Uniprot IDs
#' @param verbosity `integer` Verbosity level for uniprotREST::uniprot_map
#' @return `data.frame` object.
#' @export
#' @examples
#' uniprotIDs <- c('O76024', 'Q03135', 'Q96T23')
#' query_protein_tm_topology(uniprotIDs)
query_protein_tm_topology <- function(uniprotIDs, verbosity=0){

  uniprotREST::uniprot_map(
    ids = uniprotIDs,
    from = "UniProtKB_AC-ID", method='stream',
    to = "UniProtKB", verbosity = verbosity,
    fields = c("accession", "length", "ft_intramem", "ft_transmem", "ft_topo_dom")) %>%
    dplyr::rename(UniprotID=From) %>%
    dplyr::select(-Entry)

}


get_tm_info <- function(trans_info_string){

  if(length(trans_info_string)==0){
    return(rep(NA, 5))
  }


  start_stop <- sapply(strsplit(trans_info_string, ';'), '[[', 1)
  start_stop <- strsplit(start_stop, '\\.\\.')

  start <- as.numeric(sapply(start_stop, '[[', 1))
  end <- as.numeric(sapply(start_stop, '[[', 2))

  tm_lengths <- end-start


  tm_types <- gsub('^ /note=', '', sapply(strsplit(trans_info_string, ';', 3), '[[', 2))
  tm_types_set <- unique(tm_types)

  start <- paste(start, collapse=';')
  end <- paste(end, collapse=';')
  tm_lengths <- paste(tm_lengths, collapse=';')
  tm_types <- paste(tm_types, collapse=';')
  tm_types_set <- paste(sort(tm_types_set), collapse=';')

  if('Beta stranded' %in% tm_types_set){
    return(rep(NA, 5))
  }

  return(c(start, end, tm_lengths, tm_types, tm_types_set))
}

#' Add transmembrane domain details
#'
#' @description Parse the Transmembrane column from the output of query_protein_tm_topology to add columns decscribing the transmembrane domains.
#'
#' @param query_resul `data.frame` output of query_protein_tm_topology
#' @return `data.frame` object.
#' @export
#' @examples
#' uniprotIDs <- c('O76024', 'Q03135', 'Q96T23')
#' query_result <- query_protein_tm_topology(uniprotIDs)
#' query_result_parsed <- add_tm_info(query_result)
add_tm_info <- function(query_result){

  trans_info_split <- strsplit(gsub('^TRANSMEM ', '', query_result$Transmembrane), split='; TRANSMEM ')
  n_tms <- sapply(trans_info_split, length)

  query_result$n_tms <- n_tms

  query_result[,c('tm_start', 'tm_end', 'tm_length','tm_types',
                  'tm_types_set')] <- t(sapply(trans_info_split, get_tm_info))
  return(query_result)
}


get_topology <- function(topology_string){

  topology <- gsub('/note=', '', sapply(strsplit(topology_string, '; '), '[[', 2))

  if(length(topology_string)==0){
    return(rep(NA, 6))
  }

  n_term <- topology[1]
  c_term <- topology[length(topology)]

  topology <- paste(topology, collapse=';')

  start_stop <- sapply(strsplit(topology_string, ';'), '[[', 1)
  start_stop <- strsplit(start_stop, '\\.\\.')

  start <- as.numeric(sapply(start_stop, '[[', 1))
  end <- as.numeric(sapply(
    start_stop, FUN = function(x) ifelse(length(x)>1, x[2], as.numeric(x[1])+1)))

  top_lengths <- end-start

  start <- paste(start, collapse=';')
  end <- paste(end, collapse=';')
  top_lengths <- paste(top_lengths, collapse=';')

  return(c(topology, n_term, c_term, start, end, top_lengths))
}

#' Add topology details
#'
#' @description Parse the Topological.domain column from the output of query_protein_tm_topology to add columns decscribing the topology
#'
#' @param query_resul `data.frame` output of query_protein_tm_topology
#' @return `data.frame` object.
#' @export
#' @examples
#' uniprotIDs <- c('O76024', 'Q03135', 'Q96T23')
#' query_result <- query_protein_tm_topology(uniprotIDs)
#' query_result_parsed <- add_topology_info(query_result)
add_topology_info <- function(result){
  topology_info_split <- strsplit(gsub('^TOPO_DOM ', '', result$Topological.domain), split='; TOPO_DOM ')

  result[,c('topology', 'n_term', 'c_term',
            'topology_start', 'topology_end', 'topology_length')] <- t(
              sapply(topology_info_split, get_topology))

  return(result)

}

#' Obtain protein transmembrane-domains and topology information
#'
#' @description Given a set of UniProt IDs, this function queries UniProt to obtain
#' the information regarding transmembrane domains (TMDs) and topology and then parses
#' this information to define additional columns describining the TMDs and topology
#'
#' @param uniprotIDs `character vector` Uniprot IDs
#' @param verbosity `integer` Verbosity level for uniprotREST::uniprot_map
#' @return `data.frame` object.
#' @export
#' @examples
#' uniprotIDs <- c('O76024', 'Q03135', 'Q96T23')
#' get_protein_tm_topology(uniprotIDs)
get_protein_tm_topology <- function(uniprotIDs, verbosity=0){
  query_result <- query_protein_tm_topology(uniprotIDs, verbosity=verbosity)
  add_topology_info(add_tm_info(query_result))
}
