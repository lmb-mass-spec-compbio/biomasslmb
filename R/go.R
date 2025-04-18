#' Determine GO ancestor object
#'
#' @description For a single GO term, get its ontology and then return the
#' correct ancestor mapping object for that ontology.
#'
#' @param term `string`. A single GO.ID.
#' @param ontology `string`. A single ontology, one of: BP, MF, CC.
#'
#' @return Returns an `AnnDbBimap` object which maps GO terms (BP, MF, or
#' CC) to all ancestor terms. Each GO term is mapped to a vector of ancestor
#' GO terms.
#' @keywords internal
#' @export
determine_ancestor_function <- function(term, ontology) {
  ontology <- ontology
  if (is.na(ontology)) {
    return(NA)
  }
  if (ontology == "MF") {
    return(GO.db::GOMFANCESTOR)
  }
  else if (ontology == "CC") {
    return(GO.db::GOCCANCESTOR)
  }
  else if (ontology == "BP") {
    return(GO.db::GOBPANCESTOR)
  }
  else {
    return(NA)
  }
}

#' Determine GO offspring object
#'
#' @description For a single GO term, get its ontology and then return the
#' correct offpsring mapping object for that ontology.
#'
#' @param term `string`. A single GO.ID.
#' @param ontology `string`. A single ontology, one of: BP, MF, CC.
#'
#' @return Returns an `AnnDbBimap` object which maps GO terms (BP, MF, or
#' CC) to all offspring terms. Each GO term is mapped to a vector of offsrping
#' GO terms.
#' @keywords internal
#' @export
determine_offspring_function <- function(term, ontology) {
  if (is.na(ontology)) {
    return(NA)
  }
  if (ontology == "MF") {
    return(GO.db::GOMFOFFSPRING)
  }
  else if (ontology == "CC") {
    return(GO.db::GOCCOFFSPRING)
  }
  else if (ontology == "BP") {
    return(GO.db::GOBPOFFSPRING)
  }
  else {
    return(NA)
  }
}

#' Get all mappings for GO terms
#'
#' For a vector of GO terms, obtain all of the ancestor or offspring terms.
#'
#' @param go_ids `character vector`. GO terms to use in the format `GO:1234567`.
#' @param ontologies `named character vector`. Names = `go_ids` and values =
#' ontologies e.g. `BP`, `CC`, or `MF`.
#' @param verbose `logical`.
#' @param direction `string` Either `"ancestor"` or `"offspring"`.
#'
#' @return Returns a `named list` of `character vectors`. Names == GO terms and
#' values == vectors containing all ancestor GO terms for the particular input
#' GO term.
#'
#' @importFrom purrr map2
#' @keywords internal
#' @export
get_all_mappings <- function(go_ids, ontologies, verbose = TRUE, direction = "ancestor") {
  go2relations <- list()

  if (direction == "ancestor") {
    determineFunction <- determine_ancestor_function
  }

  else if (direction == "offspring") {
    determineFunction <- determine_offspring_function
  }

  else {
    stop("direction must be `ancestor` or `offspring`")
  }

  if (verbose) {
    print(sprintf("Getting all %s GO terms for %i observed terms. This may take a while!", direction, length(go_ids)))
  }

  result <- map2(
    ontologies,
    go_ids,
    function(y, x) {

      go_relations <- AnnotationDbi::get(x, determineFunction(x, y))

      if (inherits(go_relations, "character")) {
        go2relations[[x]] <- go_relations
      } else {
        if (verbose) {
          cat(sprintf("Could not get %s GO terms for GO.ID=%s\n", direction, x))
        }
        go2relations[[x]] <- c(NA)
      }
    }
  )
  return(result)
}

#' Expand data.frame GO terms
#'
#' Input is a data.frame with all the GO terms annotated to a single protein in
#' a single column. Output is a new data.frame with all of the GO terms for
#' that protein (annotated and ancestor).
#'
#' @param go_df `data.frame` for a single protein with a single column where
#' each row is a single GO term.
#' @param go_col `variable`. Name of column from the data.frame that
#' contains the GO terms.
#' @param go2Ancestor `named list`. Returned by \code{\link{get_all_mappings}}.
#' Names == GO terms and values == vectors containing all ancestor GO terms
#' for the particular input GO term.
#'
#' @return Returns a `data.frame` containing all ancestor GO terms for a
#' single protein.
#'
#' @keywords internal
#' @export
expand_terms <- function(go_df, go_col, go2Ancestor) {
  observed_go_ids <- unique(go_df[[go_col]]) # unique() shouldn't be req. but does not harm

  unprocessed_ids <- observed_go_ids
  all_ancestors <- observed_go_ids

  while (length(unprocessed_ids) > 0) {
    query_ancestors <- unique(go2Ancestor[unprocessed_ids[1]])
    all_ancestors <- c(all_ancestors, unlist(query_ancestors))

    # removed processed go id
    # and any which were found in ancestors (we don't want to re-query these)
    unprocessed_ids <- unprocessed_ids[-1]
    unprocessed_ids <- base::setdiff(unprocessed_ids, query_ancestors)
  }
  all_ancestors <- setdiff(unique(all_ancestors), "all")

  return(data.frame("GOID" = all_ancestors))
}

#' Get all ancestor GO terms
#'
#' Given a data.frame with with a column containing GO terms, this function
#' will output a data.frame with with the ontologies of those GO terms, their
#' description, and all their ancestor terms.
#'
#' @param go_df `data.frame`. Contains all initial GO terms for proteins
#' of interest with `columns == [feature_col, go_col]`.
#' @param feature_col `string`. Name of the column with the features, e.g.
#' `"UNIPROTKB"`.
#' @param go_col `string`. Name of the column with the GO ids, e.g. `"GO.ID"`.
#' @param return_early_debug `logical`. Stop function early and return a
#' list of GO terms and ancestor GO terms for debugging purposes.
#' @param verbose `logical`.
#'
#' @return Returns a `data.frame` containing all ancestor GO terms for all
#' proteins plus GO term descriptions and ontologies.
#' @importFrom stats setNames
#' @importFrom utils head
#' @export
get_ancestor_go <- function(go_df, feature_col = "UNIPROTKB", go_col = "GO.ID",
                            return_early_debug = FALSE, verbose = TRUE) {

  all_observed_go <- unique(go_df[[go_col]])
  all_observed_go <- all_observed_go[!is.na(all_observed_go)]

  ontologies <- AnnotationDbi::select(
    GO.db::GO.db,
    all_observed_go,
    columns = c("ONTOLOGY"),
    keytype = "GOID"
  )

  ontologies <- setNames(ontologies$ONTOLOGY, ontologies$GOID)
  if (any(is.na(ontologies))) {
    # get names of GO terms that failed to be matched
    bad.terms <- names(ontologies[is.na(ontologies) == TRUE])

    warning(
      "Warning: the following GO terms failed to match to a primary ontology ",
      "in GO.db and were removed: ",
      paste(bad.terms, collapse = ", ")
    )
    # remove NA ontologies
    ontologies <- ontologies[!is.na(ontologies)]

    # remove corresponding GO ids
    all_observed_go <- all_observed_go[!all_observed_go %in% bad.terms]
  }

  go2Ancestor <- get_all_mappings(
    all_observed_go,
    ontologies,
    direction = "ancestor",
    verbose = verbose
  )

  if (return_early_debug) {
    return(list("go_df" = go_df, "go2Ancestor" = go2Ancestor))
  }

  if (verbose) message("Expanding GO terms to include all ancestors for all entries")

  full_go_df <- go_df %>%
    filter(!is.na(all_of(go_col))) %>%
    group_by(.data[[feature_col]]) %>%
    summarise(expand_terms(across(), go_col, go2Ancestor), .groups = "drop_last")

  full_go_details <- AnnotationDbi::select(
    GO.db::GO.db,
    as.character(unique(full_go_df$GOID)),
    columns = c("TERM", "ONTOLOGY"),
    keytype = "GOID"
  )

  full_go_df <- merge(full_go_df, full_go_details, by = "GOID", all.x = TRUE)

  full_go_df <- full_go_df[, c(feature_col, "GOID", "TERM", "ONTOLOGY")]
  colnames(full_go_df)[colnames(full_go_df) == "GOID"] <- go_col

  return(full_go_df)
}

#' GO term enrichment using goseq
#'
#' @description A wrapper function around \code{\link[goseq]{goseq}} to perform
#' GO term enrichment analysis. See the \code{\link[goseq]{goseq}} documentation
#' for details. `pwf` can be made using \code{\link[goseq]{nullp}}.
#'
#' Over/underrepresented p-values are automatically
#' adjusted using `method = "BH"`. If `gene2cat` is not provided then this
#' function will default to using the Homo sapiens genome `hg19` and will
#' expect Ensembl gene IDs to have been used to construct the `pwf` input.
#'
#' @param pwf `data.frame` with 3 columns (`DEgenes` = logical,
#' `bias.data` = numeric/integer, `pwf` = numeric) and row names (usually UniProt
#' accessions, Ensembl gene IDs or similar).
#' Typically constructed using \code{\link[goseq]{nullp}}.
#' @param gene2cat `data.frame` with 2 columns containing the mapping between
#' genes (usually UniProt accessions, Ensembl gene IDs or similar) and GO terms.
#' Alternatively, a `named list` where the names are genes and each entry is
#' a `character vector` of GO terms.
#' @param ... Other arguments to be passed to \code{\link[goseq]{goseq}}.
#' @param shorten_term `logical`. Should an extra column with a substring of
#' the output GO terms be added to the output data.frame? Default is `TRUE`.
#' @param shorten_lims `integer vector` of length 2. The start and stop
#' coordinates of the substring.
#'
#' @return Returns a `data.frame` of over/underrepresented GO terms.
#' @export
get_enriched_go <- function(pwf, gene2cat = NULL, ...,
                            shorten_term = TRUE, shorten_lims = c(1L, 30L)) {
  # perform GO term enrichment with or without gene2cat
  if(!is.null(gene2cat)) {
    message(sprintf("Number of DE genes input: %i", sum(pwf$DEgenes)))
    out <- goseq::goseq(pwf = pwf, gene2cat = gene2cat, ...)
  } else {
    message(sprintf("Number of DE genes input: %i", sum(pwf$DEgenes)))
    message('gene2cat not provided. Defaulting to genome = "hg19" and id = "ensGene"')
    out <- goseq::goseq(pwf, genome = "hg19", id = "ensGene", ...)
  }

  # adjust p-values
  out$over_represented_adj_pval <- stats::p.adjust(out$over_represented_pvalue, method = "BH")
  out$under_represented_adj_pval <- stats::p.adjust(out$under_represented_pvalue, method = "BH")

  # add column with shortened terms if necessary
  if (shorten_term) {
    out$term_short <- substr(out$term, start = shorten_lims[1], stop = shorten_lims[2])
  }

  # remove GO terms without any genes assigned to them
  filter(out, .data$numDEInCat > 0)
}

#' Estimate effect size of over-representation
#'
#' @description This is a crude function to estimate the effect size of
#' over-representation i.e. we know a feature is over-represented, but we want to
#' estimate the effect size/_how_ over-represented it is. This function is typically
#' run after \code{\link{get_enriched_go}} but will work from the output of `goseq`,
#' regardless whether the functional enrichment tested was for GO terms.
#'
#' @param obj `data.frame` containing `goseq` results as generated by
#' \code{\link{get_enriched_go}} or \code{\link[goseq]{goseq}}.
#' @param pwf `data.frame` as used in \code{\link{get_enriched_go}} or
#' \code{\link[goseq]{goseq}}.
#' @param gene2cat `data.frame` as used in \code{\link{get_enriched_go}} or
#' \code{\link[goseq]{goseq}}.
#'
#' @return Returns `obj` with an extra column added called `adj_overrep`. This
#' column is calculated for each GO term by:
#'
#' `numDEInCat / numInCat / (avgTermWeight / avgNonTermWeight) / (totalDEFeatures / totalFeatures)`
#'
#' where:
#' - `numDEInCat` is the number of differentially expressed genes (aka. proteins)
#' assigned to that GO term.
#' - `numInCat` is the total number of genes (aka. proteins) annotated to that
#' GO term.
#' - `avgTermWeight` is the average `pwf$pwf` value for all the differentially
#' expressed genes that were assigned to that GO term.
#' - `avgNonTermWeight` is the average `pwf$pwf` for all the other genes supplied
#' in `pwf`.
#' - `totalDEFeatures` is the total number of differentially expressed genes
#' indicated in `pwf`.
#' - `totalFeatures` is the total number of genes indicated in `pwf`, i.e. the
#' number of rows.
#'
#' @export
estimate_overrep <- function(obj, pwf, gene2cat) {

  # if gene2cat is a list, convert to data.frame of correct format
  if (is.list(gene2cat) & !is.data.frame(gene2cat)) {
    gene2cat <- gene2cat %>%
      tibble::enframe(name = "annotation", value = "IDs") %>%
      tidyr::unnest(cols = c("IDs", "annotation")) %>%
      as.data.frame() %>%
      select(IDs, annotation)
  }

  n_de_genes <- sum(pwf$DEgenes)
  n_genes <- length(pwf$DEgenes)

  # gene2cat can have variable column names so we rely on column positions
  # 1 = gene ids e.g. UniProt accessions, 2 = GO terms
  gene2cat_subset <- gene2cat[gene2cat[, 2] %in% obj$category, 1:2]

  # filter gene2cat for GO terms present in obj, then output a named list of
  # vectors where the name is a GO term and the element is a vector of gene ids
  gene2cat_long <- with(
    gene2cat_subset,
    split(gene2cat_subset[, 1], gene2cat_subset[, 2])
  )

  # sort obj in order of GO term i.e. GO:0000002 first
  obj_sorted <- obj[match(names(gene2cat_long), obj$category), ]

  # Add unadjusted over-representation
  obj_sorted$overrep <- (obj_sorted$numDEInCat / obj_sorted$numInCat) / (n_de_genes / n_genes)

  # calculate overrepresentation score
  out <- vector(mode = "numeric", length = nrow(obj))
  for (i in seq_len(nrow(obj_sorted))) {
    out[i] <- as.numeric(obj_sorted[i, "numDEInCat"]) /
      as.numeric(obj_sorted[i, "numInCat"]) /
      (mean(pwf[rownames(pwf) %in% gene2cat_long[[i]], "pwf"]) / # term weight
         mean(pwf[!rownames(pwf) %in% gene2cat_long[[i]], "pwf"])) / # non-term weight
      (n_de_genes / n_genes)
  }

  # add overrepresentation score to input and sort in increasing p-value
  obj_sorted$adj_overrep <- out

  dplyr::arrange(obj_sorted, .data$over_represented_pvalue)
}

#' Remove redundant GO terms
#'
#' @description Given the output of \code{\link{get_enriched_go}} or
#' \code{\link[goseq]{goseq}}, remove redundant GO terms.
#'
#' @param obj `data.frame` containing `goseq` results as generated by
#' @param `go_category_col` `character` column name for GO terms
#' @param `p_value_col` `character` column name for p-values
#' \code{\link{get_enriched_go}} or \code{\link[goseq]{goseq}}.
#'
#' @return Returns `obj` with redundant GO terms filtered out.
#' @export
remove_redundant_go <- function(obj, go_category_col='category', p_value_col='over_represented_pvalue') {

  obj <- obj %>% arrange(!!sym(p_value_col))

  all_observed_go <- unique(obj[[go_category_col]]) # identify all GO terms
  all_observed_go <- all_observed_go[!is.na(all_observed_go)] # remove any NA

  # get ontologies for all GO terms (i.e. BP, MF, CC)
  ontologies <- AnnotationDbi::select(
    GO.db::GO.db,
    all_observed_go,
    keytype = c("GOID"),
    columns = c("ONTOLOGY")
  )
  ontologies <- setNames(ontologies$ONTOLOGY, ontologies$GOID)

  if (any(is.na(ontologies))) {
    # get names of GO terms that failed to be matched
    bad.terms <- names(ontologies[is.na(ontologies) == TRUE])

    warning(
      "Warning: the following GO terms failed to match to a primary ontology ",
      "in GO.db and were removed: ",
      paste(bad.terms, collapse = ", ")
    )
    # remove NA ontologies
    ontologies <- ontologies[!is.na(ontologies)]

    # remove corresponding GO ids
    all_observed_go <- all_observed_go[!all_observed_go %in% bad.terms]
  }

  # get ancestor and offspring GO terms of parent GO terms
  offspring <- get_all_mappings(all_observed_go, ontologies, verbose = FALSE, direction = "offspring")
  ancestors <- get_all_mappings(all_observed_go, ontologies, verbose = FALSE, direction = "ancestor")

  # start by retaining all input GO terms (except NA)
  retained <- all_observed_go

  # keep track of GO terms processed
  processed <- NULL

  # if any GO term has no offspring or ancestor terms, mark them as processed
  # which means they are always retained
  no_off_anc <- dplyr::setdiff(all_observed_go, dplyr::union(names(offspring), names(ancestors)))
  if (length(no_off_anc) > 0) {
    message(sprintf("No offspring or ancestors could be found for these terms: %s", no_off_anc))
    processed <- no_off_anc
  }

  # stop 'while' loop when all GO terms have been processed
  while(length(dplyr::setdiff(all_observed_go, processed)) != 0) {
    # take a GO term
    go_id <- dplyr::setdiff(all_observed_go, processed)[1]

    # find all offspring and ancestor terms = go_tree
    # (only include those also observed as over-represented GO)
    go_tree <- dplyr::union(ancestors[[go_id]], offspring[[go_id]]) %>%
      dplyr::intersect(all_observed_go) %>%
      c(go_id)

    top_go <- obj %>%
      dplyr::filter(.data[[go_category_col]] %in% go_tree) %>% # subset to terms in go_tree
      arrange(.data[[p_value_col]]) %>% # order by ascending p-value
      pull(.data[[go_category_col]]) %>% # pull out category
      head(1) # keep the top GO term

    # remove all offspring and ancestor terms within go_tree for the top GO term
    terms_to_remove <- dplyr::union(offspring[[top_go]], ancestors[[top_go]]) %>%
      dplyr::intersect(go_tree)

    processed <- dplyr::union(processed, c(top_go, terms_to_remove)) # all terms in go_tree now 'processed'
    retained <- dplyr::setdiff(retained, terms_to_remove) # remove 'processed' terms from retained
  }

  out <- obj %>%
    filter(.data[[go_category_col]] %in% retained) # keep only 'retained' terms

  out
}

#' Plot selected GO terms
#'
#' @description This function plots a set of GO terms of interest.
#' be run after \code{\link{get_enriched_go}} and \code{\link{estimate_go_overrep}}.
#' To avoid plotting too many terms, you may wish to use
#' \code{\link{remove_redundant_go}} first too. A warning is shown if you try
#' and plot more than 100 GO terms.
#'
#' @param obj `data.frame` containing `goseq` results as generated by
#' \code{\link{get_enriched_go}} then \code{\link{estimate_go_overrep}}.
#' See below for an example.
#' @param term_col `character` column name for GO term description.
#' @param ontology_col `character` column name for GO term ontology.
#' @param annotate_n `logical` Include the number of features per term
#' @param make_terms_unique `logical` Make GO terms unique. Useful if truncated GO terms are being plotted
#' (i.e. numDEInCat column) in the plot? (Default is `TRUE`)
#'
#' @return Returns a `ggplot` object.
#'
#' @examples
#' # Make a data.frame of the correct format
#' df <- data.frame(
#'   "term_short" = c("A GO term", "Another GO term", "One more GO term"),
#'   "ontology" = c("BP", "MF", "CC"),
#'   "over_represented_adj_pval" = c(0.0001, 1, 0.01),
#'   "adj_overrep" = c(15, 3, 1),
#'   "numDEInCat" = c(304, 22, 78)
#' )
#'
#' # Plot the GO terms
#' plot_go(df)
#'
#' @export
plot_go <- function(obj,
                    term_col = "term_short",
                    ontology_col = "ontology",
                    annotate_n = TRUE,
                    make_terms_unique = TRUE) {
  # throw an error if peptide input does not contain the required columns
  req_cols <- c(term_col, ontology_col, "over_represented_adj_pval",
                "adj_overrep", "numDEInCat")

  if (!all(req_cols %in% colnames(obj))) {
    stop(
      paste("The input data.frame is missing the following required columns:",
            req_cols[!req_cols %in% colnames(obj)])
    )
  }

  # throw a warning if trying to plot more than 100 GO terms
  if (nrow(obj) > 100) warning("Trying to plot more than 100 GO terms is not advisable.")

  # print message if term_col contains non-unique values
  if (any(duplicated(obj[, term_col]))) {
    message(sprintf("The column '%s' contains duplicates.", term_col))

    if(make_terms_unique){
      obj[, term_col] <- make.unique(obj[, term_col])
    }
  }

  # construct y axis text
  to_plot <- obj %>%
    arrange(desc(.data$over_represented_adj_pval)) %>%
    mutate(term_ontology = paste0(!!sym(term_col), ' (', !!sym(ontology_col), ')')) %>%
    mutate(term_ontology = factor(.data$term_ontology, levels = unique(.data$term_ontology)))

  # construct plot
  p <- to_plot %>%
    ggplot(aes(x = .data$adj_overrep, y = .data$term_ontology,
               fill = -log10(.data$over_represented_adj_pval))) +
    geom_bar(stat = 'identity') +
    theme_biomasslmb(base_size = 15, border = FALSE, aspect_square = FALSE) +
    scale_fill_continuous(limits = c(0, NA),
                          high = get_cat_palette(1),
                          low = 'grey',
                          na.value = get_cat_palette(1)) +
    theme(legend.title = element_text(size = 10),
          legend.text = element_text(size = 10),
          legend.key.size = unit(0.75, "lines")) +
    labs(x = "Over-representation",
         y = "",
         fill = "-log10(FDR)")

  if (annotate_n) {
    p <- p + geom_text(aes(label = .data$numDEInCat),
                       x = max(to_plot$adj_overrep)*1.2, hjust = 1) +
      xlim(NA, max(to_plot$adj_overrep)*1.2)
  }

  return(p)
}

#' Plot selected GO terms
#'
#' @description This function plots a set of GO terms of interest.
#' be run after \code{\link{get_enriched_go}} and \code{\link{estimate_go_overrep}}.
#' To avoid plotting too many terms, you may wish to use
#' \code{\link{remove_redundant_go}} first too. A warning is shown if you try
#' and plot more than 100 GO terms.
#'
#' @param goi `character` GO terms of interest
#' @param foi `character` features (genes/proteins etc) of interest
#' @param gene2cat `data.frame` with 2 columns containing the mapping between
#' genes (usually UniProt accessions, Ensembl gene IDs or similar) and GO terms.
#' @param term_col `character` column name for GO term description.
#' @param gene_col `character` column name for genes (or protein etc).
#'
#' @return Returns a `ggplot` object.
#'
#' @export
plot_go_terms_upset <- function(goi,
                                foi,
                                gene2cat,
                                term_col='TERM',
                                gene_col='UNIPROTKB'){

  to_plot_upset <- gene2cat %>% filter(!!sym(term_col) %in% goi) %>%
    filter(!!sym(gene_col) %in% foi) %>%
    mutate(present=1) %>%
    select(!!sym(gene_col), !!sym(term_col), present) %>%
    pivot_wider(names_from=!!sym(term_col), values_from=present) %>%
    select(-!!sym(gene_col))

  to_plot_upset[is.na(to_plot_upset)] <- 0

  p <- UpSetR::upset(data.frame(to_plot_upset),
                nsets = length(goi),
                nintersects = 50,
                mb.ratio = c(0.5, 0.5),
                order.by = "freq",
                decreasing = TRUE)
  return(p)
}




#' Perform independent filtering to threshold the over-representation testing
#' on the number of features in each category to limit the multiple testing burden
#'
#' @description When testing for enrichment of feature annotations, we frequently test
#' thousands of annotations simultaneously, which increases the multiple testing burden considerably.
#' Given many annotations have very few annotated features, there will be too little power to detect
#' over-representation. Inspired by DESeq2, which was in turn inspired by \href{https://www.pnas.org/doi/10.1073/pnas.0914005107}{Bourgon, Gentleman, and Huber (2010)},
#' this function applies the principle of 'independent filtering' to threshold the minimum
#' number of features per annotation (in the background) to optimally limit the multiple testing burden. This is
#' statistically valid because the filtering criteria (number of features per annotation) is
#' statistically independent from the test statistic under the null hypothesis.
#' The independent filtering works because the filtering criteria is correlated
#' with the test statistic under the alternative hypothesis. For further justification and discussion of independent filtering,
#' see \href{https://www.pnas.org/doi/10.1073/pnas.0914005107}{Bourgon, Gentleman, and Huber (2010)} and
#' the \href{https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#indfilttheory}{DESeq2 vignette}.
#'
#' The code function is largely lifted from `DEseq:::pvalueAdjustment`
#'
#' Please cite \href{https://www.pnas.org/doi/10.1073/pnas.0914005107}{Bourgon, Gentleman, and Huber (2010)} when using this function.
#'
#' @param obj `data.frame` containing `goseq` results as generated by
#' \code{\link{get_enriched_go}}
#' See below for an example.
#' @param alpha `numeric` alpha value to use for rejection of null hypothesis.
#' Note this is only used for optimising the threshold and any alternative value for alpha can be used downstream
#' @param theta `numeric` vector of thresholds (fractions of the data to remove)
#' @param p_value_col `character` Column with p-values
#' @param filter_col `character` Column with filtering criteria values
#' @param plot `logical` Plot the relationship between the
#'
#' @return Returns the `data.frame` with an added column denoting the optimised adjusted p-value ('padj_if')
#' @importFrom genefilter filtered_p
#'
#' @export
add_independent_filtering_padj <- function(obj,
                                           alpha=0.1,
                                           theta=seq(0, 1, .01),
                                           p_value_col='over_represented_pvalue',
                                           filter_col='numInCat',
                                           plot=TRUE){

  p <- obj[[p_value_col]]
  filter <- obj[[filter_col]]

  filtPadj <- filtered_p(filter=filter, test=p, theta=theta, method='BH')

  numRej  <- MatrixGenerics::colSums(filtPadj < alpha, na.rm = TRUE)
  # prevent over-aggressive filtering when all genes are null,
  # by requiring the max number of rejections is above a fitted curve.
  # If the max number of rejection is not greater than 10, then don't
  # perform independent filtering at all.
  lo.fit <- lowess(numRej ~ theta, f=1/5)
  if (max(numRej) <= 5) { # More relaxed than `DEseq:::pvalueAdjustment` value of 10
    j <- 1
  } else {

    residual <- if (all(numRej==0)) {
      0
    } else {
      numRej[numRej > 0] - lo.fit$y[numRej > 0]
    }

    # this usually works: find the threshold at which num rejections
    # surpasses the root mean squared error around the fitted curve.
    # it may not work if there is a sharp uptick in the curve at
    # the end of the grid, and there is very little variation.
    maxFit <- max(lo.fit$y)
    rmse <- sqrt(mean(residual^2))
    thresh <- maxFit - rmse

    j <- if (any(numRej > thresh)) {
      # backup case: if low variation and uptick at end,
      # pick the first point at which num rejections reaches
      # 90% of the fitted curve, or 80% of the fitted curve
      which(numRej > thresh)[1]
    } else if (any(numRej > .9 * maxFit)) {
      which(numRej > .9 * maxFit)[1]
    } else if (any(numRej > .8 * maxFit)) {
      which(numRej > .8 * maxFit)[1]
    } else {
      1
    }
  }

  padj <- filtPadj[, j, drop=TRUE]

  if(plot){
    cutoffs <- quantile(filter, theta)
    filterThreshold <- cutoffs[j]
    filterNumRej <- data.frame(theta=theta, numRej=numRej)
    filterTheta <- theta[j]

    plot(filterNumRej,
         type="b", ylab="number of rejections",
         xlab="quantiles of filter")
    lines(lo.fit, col="red")
    abline(v=filterTheta)

    plot(log2(cutoffs), filterNumRej$numRej,
         type="b", ylab="number of rejections",
         xlab="Threshold value (log2)")
    abline(v=log2(filterThreshold))
    mtext(filterThreshold, side=3, at=log2(filterThreshold))
    }

  obj$padj_if <- padj

  return(obj)
}
