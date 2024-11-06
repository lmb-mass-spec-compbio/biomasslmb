#' Plot the correlation between sample
#'
#' @description The correlation between feature quantifications in each sample can
#' be a useful QC step to assess whether the experiment has worked. This function plots
#' the correlation values in a heatmap. Samples are not clustered and are ordered
#' on the axes in the order they are in the SummarizedExperiment
#'
#' @param obj `QFeatures`. Proteomics dataset
#' @param i `string`. Index for the SummarizedExperiment you wish to plot
#' @param ... addiional arguments passed onto corrplot::corrplot
#' @return Returns a _ggplot_ object.
#' @export
#' @importFrom corrplot corrplot
#' @examples
#'
#' tmt_qf <- QFeatures::readQFeatures(assayData = psm_tmt_total,
#'   quantCols = 36:45,
#'   name = "psms_raw")
#'
#' cor_sample(tmt_qf, 'psms_raw')
plot_cor_samples <- function(obj, i, order='original', ...){

  check_q(obj)

  check_se_exists(obj, i)

  corrplot(cor(assay(obj[[i]]),use='complete', method='spearman'),
           method = 'circle',
           type = 'upper',
           order=order,
           tl.col='grey20',
           tl.cex=1,
           cl.cex=1,
           ...)

}

