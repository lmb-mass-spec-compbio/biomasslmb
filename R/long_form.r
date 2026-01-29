#' Compatibility wrapper for QFeatures long-format extraction
#'
#' Dispatches to \code{longForm()} or
#' \code{longFormat()}, depending on which
#' function is available in the QFeatures namespace.
#'
#' This helper exists to maintain compatibility across QFeatures
#' releases without relying on Bioconductor version checks.
#' The function is evaluated at runtime and works regardless of
#' whether \pkg{QFeatures} is attached.
#'
#' @param ... Arguments passed directly to
#'   \code{QFeatures::longForm()} or \code{QFeatures::longFormat()}.
#'
#' @return
#' A data structure identical to the return value of the underlying
#' QFeatures long-format function (typically a \code{LongTable}
#' or \code{DataFrame}).
#'
#' @details
#' If both functions were to exist, \code{longForm()} is preferred.
#' An error is raised if neither function is found.
#'
#' @seealso
#' \code{\link[QFeatures]{longForm}},
#' \code{\link[QFeatures]{longFormat}}
#'
#' @keywords internal
qfeatures_long <- function(...) {
  ns <- asNamespace("QFeatures")

  if (exists("longForm", envir = ns, inherits = FALSE)) {
    return(ns$longForm(...))
  }

  if (exists("longFormat", envir = ns, inherits = FALSE)) {
    return(ns$longFormat(...))
  }

  stop(
    "Neither 'longForm' nor 'longFormat' is available in the QFeatures namespace."
  )
}
