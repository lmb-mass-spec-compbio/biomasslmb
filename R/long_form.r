#' Compatibility wrapper for QFeatures long-format extraction
#'
#' Dispatches to the appropriate QFeatures long-format function
#' depending on the Bioconductor version and object class.
#'
#' Specifically, this function:
#' * Uses \code{longForm()} (BiocGenerics generic) if a method exists
#'   for the input object,
#' * Falls back to \code{QFeatures::longFormat()}
#'   otherwise.
#'
#' This allows package or script code to remain compatible across
#' QFeatures releases without explicitly checking Bioconductor versions.
#'
#' @param object A QFeatures object (or object class supported by longForm/longFormat)
#' @param ... Additional arguments passed to the underlying long-format function,
#'   typically including \code{colvars} and \code{rowvars}.
#'
#' @return A data structure identical to the output of the underlying
#'   long-format function (usually a \code{LongTable} or \code{DataFrame}).
#'
#' @details
#' The function first checks whether \code{longForm()} is a generic and
#' whether there is a registered method for the input object class.
#' If so, it calls \code{longForm()}. Otherwise, it falls back to the legacy
#' \code{longFormat()} function. An error is thrown if neither is available.
#'
#' @seealso
#' \code{\link[QFeatures]{longFormat}},
#' \code{\link[QFeatures]{longForm}},
#' \code{\link[methods]{isGeneric}},
#' \code{\link[methods]{findMethods}}
#'
#' @keywords internal
qfeatures_long <- function(object, ...) {

  if (methods::isGeneric("longForm") &&
      length(methods::findMethods("longForm", classes = class(object))) > 0) {
    return(longForm(object, ...))
  }

  if (exists("longFormat", where = asNamespace("QFeatures"), inherits = FALSE)) {
    return(QFeatures::longFormat(object, ...))
  }

  stop(
    paste0(
      "Neither a longForm() method nor longFormat() is available for object of class: ",
      paste(class(object), collapse = ", ")
    )
  )
}
