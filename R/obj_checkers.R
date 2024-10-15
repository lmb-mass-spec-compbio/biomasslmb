#' @noRd
check_null <- function(obj){
  if(is.null(obj)){
    stop("`obj` is a null object. Did you mistype the QFeatures experiment name?")
  }
}

#' @noRd
check_se <- function(obj){

  check_null(obj)

  if(class(obj)!="SummarizedExperiment"){
    stop("`obj` must be a SummarizedExperiment object. E.g a single experiment from a QFeatures object")
  }
}

#' @noRd
check_se_psm <- function(obj){

  check_null(obj)

  if(class(obj)!="SummarizedExperiment"){
    stop("`obj` must be a SummarizedExperiment object. It should contain the PSM-level quantification, e.g obj[['psm']].")
  }
}

#' @noRd
check_se_peptide <- function(obj){

  check_null(obj)

  if(class(obj)!="SummarizedExperiment"){
    stop("`obj` must be a SummarizedExperiment object. It should contain the Peptide-level quantification, e.g obj[['psm']].")
  }
}

#' @noRd
check_se_protein <- function(obj){

  check_null(obj)

  if(class(obj)!="SummarizedExperiment"){
    stop("`obj` must be a SummarizedExperiment object. It should contain the Protein-level quantification, e.g obj[['psm']].")
  }
}

