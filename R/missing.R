#' Plot the most common missing value patterns
#'
#' @description The patterns in missing values can be informative with respect to
#' whether the experiment has worked, or if particular samples are outliers. This
#' function uses an 'upset' plot to show the top 50 most common missing value patterns
#' across the samples
#'
#' @param obj `QFeatures`. Proteomics dataset
#' @param i `string`. Index for the SummarisedExperiment you wish to plot
#'
#' @return Returns a _ggplot_ object.
#' @export
#' @importFrom naniar gg_miss_upset
#' @examples
#' set.seed(11)
#' library(ggplot2)
#'
#' df <- diamonds[sample(nrow(diamonds), 1000), ]
#'
#' tmt_qf <- QFeatures::readQFeatures(assayData = psm_tmt_total,
#'   quantCols = 36:45,
#'   name = "psms_raw")
#'
#'
#'
#'
#'
#'
plot_missing_upset <- function(obj, i){

  check_q(obj)

  check_se_exists(obj, i)

  missing_data <- obj[[i]] %>%
    assay() %>%
    data.frame()

  p <- gg_miss_upset(missing_data,
                     sets = paste0(colnames(missing_data), '_NA'),
                     keep.order = TRUE,
                     nintersects = 50)
  return(p)
}
