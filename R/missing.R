#' Plot the most common missing value patterns
#'
#' @description The patterns in missing values can be informative with respect to
#' whether the experiment has worked, or if particular samples are outliers. This
#' function uses an 'upset' plot to show the most common missing value patterns
#' across the samples, the top 50 by default
#'
#' @param obj `QFeatures`. Proteomics dataset
#' @param i `string`. Index for the SummarisedExperiment you wish to plot
#' @param ... additional arguments passed onto `naniar::gg_miss_upset`, and from
#' there onto `UpSetR::upset`. Arguments given here override the defaults set by
#' this function (`sets`, `keep.order` and `nintersects`).
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
#' # just the ten most common patterns
#' plot_missing_upset(tmt_qf, "psms_raw", nintersects = 10)
#'
plot_missing_upset <- function(obj, i, ...){

  check_q(obj)

  check_se_exists(obj, i)

  missing_data <- obj[[i]] %>%
    assay() %>%
    data.frame()

  missing_data <- missing_data[,colSums(is.na(missing_data))>0]
  missing_data <- missing_data[,sort(colnames(missing_data))]

  upset_args <- list(data = missing_data,
                     sets = paste0(colnames(missing_data), '_NA'),
                     keep.order = TRUE,
                     nintersects = 50)

  dots <- list(...)
  upset_args[names(dots)] <- dots

  p <- do.call(gg_miss_upset, upset_args)
  return(p)
}
