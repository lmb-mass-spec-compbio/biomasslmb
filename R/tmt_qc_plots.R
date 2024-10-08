
#' Plot the missing values vs signal:noise
#'
#' @description Missing values are more frequent with low signal:noise (S:N).
#' This function visualises this relationship to aid selection of thresholds for
#' minimal S:N filtering
#'
#' @param obj `MSnSet` containing PSM level TMT intensities
#' @param sn_column `character` column name for signal:noise values
#' @param bins `numeric` Number of bins to plot
#'
#' @return `ggplot` stacked bar plot to show S:N vs # missing values
#' @export
#' @import dplyr
#' @import tidyr
#' @importFrom grDevices colorRampPalette
plot_missing_SN <- function(obj,
                          sn_column="Average.Reporter.SN",
                          bins=20){


  check_se_psm(obj)

  pal <- colorRampPalette(get_cat_palette(2))

  n_missing <- obj |> SummarizedExperiment::assay() |> is.na() |> rowSums()

  p <- data.frame('n_missing'=n_missing,
                  'sn'=rowData(obj)[[sn_column]]) %>%
    mutate(binned_sn=cutr::cutf2(sn, g=bins, digits=1)) %>%
    filter(is.finite(sn)) %>%
    group_by(n_missing, binned_sn) %>%
    tally() %>%
    ggplot(aes(binned_sn, n, fill=factor(n_missing))) +
    geom_bar(stat='identity', position='fill', colour='grey20', linewidth=0.1) +
    scale_fill_manual(values=c('grey70', pal(max(n_missing))), name='Missing values') +
    guides(colour=guide_legend(override.aes = list(size = 1.5))) +
    xlab('Signal:Noise') +
    scale_y_continuous(name='Proportion', expand = c(0, 0)) +
    theme_biomasslmb(base_size=12, border=FALSE) +
    theme(axis.text.x=element_text(size=10, angle=45, vjust=1, hjust=1))

  return(p)

}


#' Plot the missing values vs signal:noise for each sample
#'
#' @description Missing values are more frequent with low signal:noise (S:N).
#' This function visualises this relationship for each sample to aid selection
#' of thresholds for minimal S:N filtering.
#'
#' @param obj `MSnSet` containing PSM level TMT intensities
#' @param sn_column `character` column name for Signal:noise values
#' @param bins `numeric` Number of bins to plot
#'
#' @return `ggplot` tile plot to show S:N vs # missing values for each sample
#' @importFrom magrittr %>%
#' @export
plot_missing_SN_per_sample <- function(obj,
                                       sn_column="Average.Reporter.SN",
                                       bins=20){

  check_se_psm(obj)

  sn_per_sample <- obj %>%
    SummarizedExperiment::assay() %>%
    data.frame() %>%
    tibble::rownames_to_column('PSM_id') %>%
    gather(key='sample', value="value", -"PSM_id") %>%
    merge(rowData(obj)[,sn_column, drop=FALSE], by.x='PSM_id', by.y='row.names') %>%
    filter(is.finite(Average.Reporter.SN))

  sn_per_sample$binned_sn <- cutr::cutf2(sn_per_sample[[sn_column]], g=bins, digits=1)

  p <- sn_per_sample %>%
    group_by(sample, binned_sn,
             missing=ifelse(is.na(value), 'missing', 'present')) %>%
    tally() %>%
    spread(key=missing, value=n, fill=0) %>%
    mutate(percentage_missing=(100*missing)/(missing+present),
                  sample=remove_x(sample)) %>%
    ggplot(aes(binned_sn, sample, fill=percentage_missing)) +
    geom_tile(colour='grey20', linewidth=0.1) +
    scale_fill_gradient(low='grey97', high=get_cat_palette(1), name='Missing (%)',
                                 limits=c(0,100)) +
    theme_biomasslmb(base_size=10) +
    theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1)) +
    labs(x='Signal:Noise', y='Sample')

  return(p)
}





