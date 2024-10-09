#' Identify how many features (PSMs/Peptides) are quantified for each protein
#'
#' @description For summarisation of PSM or peptide to protein, we need a
#' minimum number of finite values per protein per sample. This function simply tallies how many we have.
#'
#'
#'
#' @param obj `SummarizedExperiment` with PSM or peptide-level quantification
#' @param master_protein_col `character` Column name for master protein
#'
#' @return `data.frame` detailing how many features are present for each protein in each sample
get_n_feature_per_prot <- function(obj,
                                   master_protein_col = "Master.Protein.Accessions") {
  feature2protein <- rowData(obj) %>%
    data.frame() %>%
    dplyr::select(!!sym(master_protein_col)) %>%
    tibble::rownames_to_column('feature_ID')

  n_feature_per_prot <- assay(obj) %>%
    data.frame() %>%
    tibble::rownames_to_column('feature_ID') %>%
    gather(key = 'sample', value = 'value', -"feature_ID") %>%
    merge(feature2protein, by = "feature_ID") %>%
    filter(is.finite(.data$value)) %>%
    group_by(!!sym(master_protein_col), sample) %>%
    tally() %>%
    ungroup() %>%
    data.frame()

  return(n_feature_per_prot)
}

#' Remove features which are assigned to a protein with too few supporting
#' features in total
#'
#' @description For summarisation of PSM or peptide to protein, we need a
#' minimum number of finite values per protein per sample.
#'
#'
#' @param obj `SummarizedExperiment` with PSM or peptide-level quantification
#' @param min_features `numeric` Threshold for minimum features per protein
#' @param master_protein_col `character` Column name for master protein
#' @param plot Set TRUE to plot histogram of features per protein per sample
#'
#' @return `SummarizedExperiment`
#' @export
filter_features_per_protein <- function(obj,
                                        min_features,
                                        master_protein_col = "Master.Protein.Accessions",
                                        plot = FALSE) {

  check_se_psm(obj)

  n_feature_per_prot <- get_n_feature_per_prot(obj)

  if (plot) {
    p <- ggplot(n_feature_per_prot, aes(log2(n))) +
      geom_histogram() +
      theme_biomasslmb() +
      xlab('# features (log2)')
    print(p)
  }

  retain_proteins <- n_feature_per_prot %>%
    filter(n >= min_features) %>%
    pull(!!sym(master_protein_col))

  return_obj <- obj[rowData(obj)[[master_protein_col]] %in% retain_proteins, ]

  invisible(return_obj)
}

#' Identify proteins which have too few features to quantify protein abundance in each sample
#'
#' @description
#'
#' Protein level abundances are more accurately quantified where there are too more
#' features (PSMs/peptides) to summarise from.
#'
#' Usually, we are performing the summarisation
#' from a matrix (columns=samples, rows=features) with an associated feature to protein ID mapping.
#' Within the matrix, some values may be missing (NA). In order to correctly identify which proteins can
#' be quantified, we need to start from the feature level object and create a mask
#' which we can use to replace protein-level quantification values with NA.
#' This is what this function does. This can then be combined with the \code{\link{mask_protein_quant}} function
#' to replace protein level quantification values with NA where they were derived from too few quantification values
#'
#' @param obj `SummarizedExperiment` with PSM or peptide-level quantification
#' @param min_features `numeric` Threshold for minimum features per protein
#' @param master_protein_col `character` Column name for master protein
#' @param plot Set TRUE to plot how many proteins are quantified in each sample.
#'  Horizontal line represents total number of proteins quantified across all samples
#'
#' @return `Matrix` defining whether the protein was quantified from sufficient features
#' @export
get_protein_no_quant_mask <- function(obj,
                                      min_features,
                                      master_protein_col = "Master.Protein.Accessions",
                                      plot = FALSE) {

  check_se_psm(obj)

  n_feature_per_prot <- get_n_feature_per_prot(obj)

  protein_quant_mask <- n_feature_per_prot %>%
    pivot_wider(names_from=sample, values_from=n, values_fill=0) %>%
    tibble::column_to_rownames(master_protein_col)

  protein_quant_mask <- protein_quant_mask>=min_features
  protein_quant_mask[!protein_quant_mask] <- NA

  if(plot){
    p <- protein_quant_mask %>%
      colSums(na.rm=TRUE) %>%
      data.frame() %>%
      tibble::rownames_to_column('Sample') %>%
      ggplot(aes(Sample, .)) + geom_bar(stat='identity') +
      theme_biomasslmb(border=FALSE) +
      theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1)) +
      ylab('Proteins quantified') +
      geom_hline(yintercept=sum(apply(protein_quant_mask, MARGIN=1, any)), linetype=2, colour='grey60')

    print(p)
  }

  return(protein_quant_mask)

}

#' Replace protein-level quantifications with NA if they derive from too few lower feature level quantifications
#'
#' @description
#'
#' Protein level abundances are more accurately quantified where there are more
#' features (PSMs/peptides) to summarise from.
#'
#' Usually, we are performing the summarisation
#' from a matrix (columns=samples, rows=features) with an associated feature to protein ID mapping.
#' Within the matrix, some values may be missing (NA). In order to correctly identify which proteins can
#' be quantified, we need to start from the feature level object and create a mask
#' which we can use to replace protein-level quantification values with NA.
#' This is can be obtained with \code{\link{get_protein_no_quant_mask }}. This can then be combined with this function
#' to replace protein level quantification values with NA where they were derived from too few quantification values
#'
#' @param obj `SummarizedExperiment` with PSM or peptide-level quantification
#' @param retain_mask `matrix` detailing whether a protein-level quantification had sufficient lower level
#'  quantification values for each sample. Can be obtained with \code{\link{get_protein_no_quant_mask }}
#'
#' @return `SummarizedExperiment` with quantification values replaced by NA where they derive from too few lower feature level quantifications
#' @export
mask_protein_level_quant <- function(obj, retain_mask){

  check_se_protein(obj)

  colnames(retain_mask) <- remove_x(colnames(retain_mask))

  retain_mask <- retain_mask[rownames(obj), colnames(obj)]
  masked_exprs <- assay(obj) * retain_mask

  assay(obj) <- as.matrix(masked_exprs)

  return(obj)
  }
