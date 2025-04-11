#' Extract the number of samples each feature was detected in for each experiment in a Qfeatures object
#'
#' @description When exploring the proteome coverage, it is instructive to consider how many samples a protein
#' was present in for each level of filtering. This function takes a Qfeatures object and tallies how many samples a variable in the `rowData`
#' was present in. When the row variables are the protein names, the tally is thus the number of samples the protein was present in.
#'
#' @param obj `QFeatures` object
#' @param rowVars `character vector` row variables to group by, e.g 'Master.Protein.Accessions'
#' @param rename_cols `named character vector` optional named list to rename the experiments. List values should be current experiment names and list names should be updated experiment names
#' @return `data.frame` object.
#' @export
get_samples_present <- function(obj, rowVars, rename_cols=NULL){

  check_q(obj)

  samples_present <- longFormat(obj, rowvars=rowVars) %>%
    data.frame() %>%
    filter(is.finite(value)) %>%
    group_by_at(.vars=c('assay', 'colname', rowVars)) %>%
    tally()  %>%
    group_by_at(.vars=c('assay', rowVars)) %>%
    tally() %>%
    pivot_wider(names_from=assay, values_from = n)

  if(!is.null(rename_cols)){

    samples_present <- samples_present %>%
      select(all_of(rowVars), all_of(rename_cols))

  }

  else{

    samples_present <- samples_present %>%
      select(all_of(rowVars), names(obj))

  }

  ungroup(samples_present)
}

#' Plots the number of samples each feature was detected in for each experiment in a Qfeatures object
#'
#' @description When exploring the proteome coverage, it is instructive to consider how many samples a protein
#' was present in for each level of filtering. This function plots the output of `get_samples_present`
#'
#' @param samples_present `data.frame` output from `get_samples_present`
#' @param rowVars `character vector` row variables on which the tallying was grouped by
#' @return `ggplot` object.
#' @export
plot_samples_present <- function(samples_present, rowvars, breaks=NULL){

  n_samples <- select(samples_present, -all_of(rowvars)) %>% max(na.rm=TRUE)

  p <- samples_present %>%
    pivot_longer(cols=-all_of(rowvars)) %>%
    group_by(name, value) %>%
    tally() %>%
    mutate(value=ifelse(value==0, NA, value)) %>%
    mutate(name=factor(name, levels=colnames(samples_present))) %>%
    arrange(!is.na(value)) %>%
    ggplot(aes(name, n, fill=value)) +
    geom_bar(stat='identity', position='stack', colour='grey30') +
    theme_biomasslmb() +
    theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1)) +
    labs(x='Analysis stage', y='Proteins') +
    scale_fill_gradient2(low='grey99', mid=get_cat_palette(2)[2], high=get_cat_palette(5)[5], midpoint=floor(n_samples/2),
                         na.value='grey40', name='Samples', limits=c(0, n_samples))

  if(!is.null(breaks)){
    p <- p + scale_fill_gradient2(low='grey99', mid=get_cat_palette(2)[2], high=get_cat_palette(5)[5], midpoint=floor(n_samples/2),
                              na.value='grey40', name='Samples', limits=c(0, n_samples), breaks=breaks)
  }

  return(p)
}
