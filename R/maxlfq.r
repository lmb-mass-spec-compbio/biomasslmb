#' Wrapper for iq::maxLFQ to use with QFeatures::aggregateFeatures
#'
#' This function wraps `iq::maxLFQ()` so it can be passed to
#' `QFeatures::aggregateFeatures()` as the `fun` argument. It expects
#' a numeric matrix with peptides as rows and samples as columns, and
#' returns a numeric vector of protein-level intensities.
#'
#' @param mat A numeric matrix of peptide intensities (rows = peptides,
#'   columns = samples). Missing values should be `NA`.
#' @param ... Additional arguments passed to `iq::maxLFQ()`.
#'
#' @return A numeric vector of length equal to the number of columns in
#'   `mat`, containing protein-level abundance estimates.
#'
#' @details
#' The MaxLFQ algorithm computes protein abundances from peptide
#' intensities by using peptide ratios across samples and solving a
#' least-squares system. This wrapper extracts the `estimate` component
#' returned by `iq::maxLFQ()` and coerces it to a numeric vector,
#' suitable for `aggregateFeatures()`.
#'
#' @examples
#' \dontrun{
#' library(QFeatures)
#' library(iq)
#'
#' # Example peptide-level matrix (log2 intensities)
#' mat <- matrix(rnorm(20), nrow = 5, ncol = 4)
#' colnames(mat) <- paste0("Sample", 1:4)
#' rownames(mat) <- paste0("Pep", 1:5)
#'
#' protein_estimates <- maxlfq_wrapper(mat)
#' }
#'
#' @export
maxlfq_wrapper <- function(mat, ...) {
  if (!requireNamespace("iq", quietly = TRUE)) {
    stop("Package 'iq' is required for maxLFQ. Please install it with install.packages('iq').")
  }

  # Ensure input is a numeric matrix
  mat <- as.matrix(mat)
  if (!is.numeric(mat)) stop("'mat' must be a numeric matrix")

  # Call iq::maxLFQ
  res <- iq::maxLFQ(mat, ...)

  # Return the numeric estimate vector
  as.numeric(res$estimate)
}
