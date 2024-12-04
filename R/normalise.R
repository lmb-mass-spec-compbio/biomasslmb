#' Extract the assay column medians from an MSnSet
#'
#' @param obj `SummarisedExperiment`
#' @return `vector` of assay column medians in colnames order
#' @importFrom robustbase colMedians
#' @export
get_medians <- function(obj){
  check_se(obj)
  medians <- robustbase::colMedians(assay(obj), na.rm=TRUE)
  return(medians)
}


#' Center-median normalise the expression matrix in an MSnSet using medians
#' from a reference dataset
#'
#' @description Center-median normalisation is a simple normalisation method
#' that is appropriate for relative abundance proteomics such as isobaric tagging.
#' This can be achieved with `QFeatures::normalize(method='diff.median')`. However,
#' for some experimental designs, the normalisation should be against the medians
#' in another dataset. For example, for PTM studies, one may wish to isobaric
#' tag samples, pool, and then PTM-enriched, with the enriched sample quantified
#' in a separate run to the non-enriched (total) sample. In this case, it may make
#' more sense to center-median normalise the PTM-enriched samples using the median
#' from the total samaples.
#'
#' @param obj `SummarisedExperiment`
#' @param medians `vector, numeric`. Sample medians from reference dataset
#' @param center_to_zero `logical`. Centre the data range on zero.
#' If FALSE, normalisation retains original data range.
#' @param on_log_scale `logical`. Input data is log-transformed
#' @return Returns a `SummarisedExperiment` with the assay data column center-median
#' normalised
#'
#' @export
center_normalise_to_ref <- function(obj, medians,
                                    center_to_zero=FALSE,
                                    on_log_scale=FALSE){

  check_se(obj)

  if(!center_to_zero){
    if(on_log_scale){
      medians <- medians - mean(medians)
      } else{
        medians <- medians/mean(medians)
      }
  }

  if(on_log_scale){
    assay(obj) <- t(t(assay(obj)) - medians)
  } else{
    assay(obj) <- t(t(assay(obj)) / medians)
  }

  return(obj)
}
