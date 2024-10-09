#' @noRd
check_se <- function(obj){
  if(class(obj)!="SummarizedExperiment"){
    stop("`obj` must be a SummarizedExperiment object. E.g a single experiment from a QFeatures object")
  }
}

#' @noRd
check_se_psm <- function(obj){
  if(class(obj)!="SummarizedExperiment"){
    stop("`obj` must be a SummarizedExperiment object. It should contain the PSM-level quantification, e.g obj[['psm']].")
  }
}

#' @noRd
check_se_protein <- function(obj){
  if(class(obj)!="SummarizedExperiment"){
    stop("`obj` must be a SummarizedExperiment object. It should contain the Protein-level quantification, e.g obj[['psm']].")
  }
}

