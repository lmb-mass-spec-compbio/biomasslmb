#' @noRd
check_se_psm <- function(obj){
  if(class(obj)!="SummarizedExperiment"){
    stop("`obj` must be a SummarizedExperiment object. Typically it would be the PSM-level quantification, e.g obj['psm']")
  }
}
