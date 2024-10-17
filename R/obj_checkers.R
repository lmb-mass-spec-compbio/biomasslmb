#' @noRd
check_null <- function(obj, object_name='obj'){
  if(is.null(obj)){
    stop(sprintf("Summarised experiment `%s` is a null object. Did you mistype the QFeatures experiment name?", object_name))
  }
}

#' @noRd
check_q <- function(obj){

  check_null(obj)

  if(class(obj)!="QFeatures"){
    stop("`obj` must be a QFeatures object")
  }
}



#' @noRd
check_se_exists <- function(obj, i){

  check_null(obj[[i]], object_name=i)

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


#' @noRd
check_colData_col <- function(obj, col){
  if(! col %in% colnames(colData(obj))){
    stop(sprintf("column %s is not in the colData. Available columns are: %s", col, paste(colnames(colData(obj)), collapse=',')))
  }
}

