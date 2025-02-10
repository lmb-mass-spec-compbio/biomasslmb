#' Plot quantification at multiple levels
#'
#' @description For the exploration of specific proteins, it's often useful to visualise the quantification at multiple levels,
#' e.g PSM, peptide, filtered peptides and protein.
#' This function plots the quantitative data at multiple levels (assays) for a given list of proteins of interest.
#'
#' @param obj `QFeatures` object
#' @param poi `list`. Proteins of interest
#' @param experiments_to_plot `list`. experiments (assays) to plot. Defaults to all assays
#' @param protein_id_col `string`. Column with the protein ids to search for `poi` values
#' @param label_col `string`. Column with labels to use for proteins in plot
#' @param rename_labels `named list`. Mapping from labels to renamed labels
#' @param log2transform_cols `list`. Assays which need to be log2-transforms (all values should ultimately be transformed)
#' @param norm_quant `logical`. Should the quantifications be normalised to the fold-change vs mean abundance
#' @param add_mean_summary `logical`. Add a line summarising the mean value over all features in the assay (Not recommended if norm_quant=FALSE)
#' @param colour_assays `logical`. Column each assay
#' @param alpha_assays `logical`. Set sensible alpha value for each assay
#' @return `ggplot` object.
#' @export

plot_protein_assays <- function(obj,
                                poi,
                                experiments_to_plot=NULL,
                                protein_id_col='Master.Protein.Accessions',
                                label_col='Master.Protein.Accessions',
                                rename_labels=NULL,
                                log2transform_cols='',
                                norm_quant=FALSE,
                                add_mean_summary=FALSE,
                                colour_assays=FALSE,
                                alpha_assays=TRUE){

  check_q(obj)

  if(is.null(experiments_to_plot)){
    experiments_to_plot <- names(obj)
  }

  to_plot <- obj[, , experiments_to_plot] %>%
    longFormat(rowvars=c(protein_id_col, label_col))  %>%
    as_tibble() %>%
    filter(!!sym(protein_id_col) %in% poi) %>%
    mutate(label=!!sym(label_col)) %>%
    mutate(value=ifelse(assay %in% log2transform_cols, log2(value), value)) %>%
    mutate(assay_order = factor(assay,
                                levels = experiments_to_plot))

  if(!is.null(rename_labels)){
    to_plot <- to_plot %>%
      mutate(label=sapply(label, function(x) rename_labels[[x]]))
  }

  if(norm_quant){
    to_plot <- to_plot %>%
      group_by(label, assay_order, rowname) %>%
      mutate(value=value-mean(value, na.rm=TRUE)) %>%
      ungroup()
    ylab <- 'Abundance relative to mean (log2FC)'
  } else{
    ylab <- 'Abundance (log2)'
  }

  n_features <- to_plot %>%
    select(assay, rowname, label) %>%
    distinct() %>%
    group_by(assay, label) %>%
    tally() %>%
    # n used to define alpha in plot. This prevents alpha getting too low!
    mutate(n=ifelse(n>100, 100, n))

  p <- to_plot %>%
    merge(n_features, by=c('assay', 'label')) %>%
    ggplot(aes(x = colname, y = value)) +
    geom_point() +
    geom_line(aes(group = rowname)) +
    facet_grid(assay_order~label, scales='free_y')  +
    theme_biomasslmb(border=FALSE, base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1, size=8),
          strip.text.y = element_text(size=10, angle=0),
          strip.background=element_blank(),
          panel.spacing.x=unit(10, 'mm')) +
    labs(x='', y=ylab)

  if(alpha_assays){
    p <- p + aes(alpha=1/n) + scale_alpha_identity()
  }
  if(colour_assays){
    p <- p + aes(colour=assay)
    if(length(experiments_to_plot)<=12){
      p <- p +
        scale_colour_manual(values=get_cat_palette(length(experiments_to_plot)),
                            guide=FALSE)
    }
  }

  if(add_mean_summary){
    p <- p +
      stat_summary(geom='line', fun='mean', aes(group=assay), colour='grey30', alpha=1)
  }

  return(p)
}
