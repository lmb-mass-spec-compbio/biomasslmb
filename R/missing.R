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
#' tmt_qf <- QFeatures::readQFeatures(assayData = psm_tmt_total,
#'   quantCols = 36:45,
#'   name = "psms_raw")
#'
#' # which combinations of samples are missing together?
#' plot_missing_upset(tmt_qf, "psms_raw")
#'
plot_missing_upset <- function(obj, i){

  check_q(obj)

  check_se_exists(obj, i)

  missing_data <- obj[[i]] %>%
    assay() %>%
    data.frame()

  missing_data <- missing_data[,colSums(is.na(missing_data))>0]
  missing_data <- missing_data[,sort(colnames(missing_data))]

  p <- gg_miss_upset(missing_data,
                     sets = paste0(colnames(missing_data), '_NA'),
                     keep.order = TRUE,
                     nintersects = 50)
  return(p)
}
