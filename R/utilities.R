#' Remove duplicated full stops
#'
#' @description A convenience function to remove any duplicated full stops
#' (aka periods) from the elements of a vector, e.g. will change `..` or `...`
#' or `....` etc. to just `.` Mainly used to fix column names.
#'
#' @param x `character` or `string`. Contains duplicate full stops to be removed.
#'
#' @return Returns `character` or `string` with duplicate full stops removed.
#' @examples
#'
#' df <- data.frame(
#'   column...name = c(1, 2, 3)
#' )
#'
#' colnames(df) <- remove_dots(colnames(df))
#'
#' @export
remove_dots <- function(x) {
  gsub("(?<=\\.)\\.+", "", x, perl = TRUE)
}

#' Remove leading X
#'
#' @description A convenience function to remove a leading capital X.
#' Is case sensitive.
#'
#' @param x `character` or `string`.
#' @return Returns `character` or `string` with leading X removed.
#' @examples
#'
#' df <- data.frame('X1'=c(1,2))
#'
#' remove_x(colnames(df))
#'
#' @export
remove_x <- function(x) {
  gsub("^X", "", x)
}

#' Report how many features and master proteins remain
#'
#' @description Prints the number of rows and the number of distinct master
#' proteins in a feature-level annotation table, followed by a short note
#' saying which step the count describes. The `filter_features_*` functions call this
#' after each filter they apply, so the same message format can be used to
#' report a count at a point where no filter function ran, such as after
#' `filterNA()` or a manual subset.
#'
#' @param x `data.frame` or `DataFrame`. Feature-level annotations, normally
#'   the output of `rowData()` on an assay.
#' @param column `string`. Name of the column in `x` holding the master protein
#'   accession, e.g. `"Master.Protein.Accessions"`.
#' @param note `string`. Short description of the step being reported.
#'
#' @return Invisibly `NULL`; called for the message it prints.
#' @examples
#' tmt_qf <- QFeatures::readQFeatures(assayData = psm_tmt_total,
#'   quantCols = 36:45,
#'   name = "psms_raw")
#'
#' message_parse(SummarizedExperiment::rowData(tmt_qf[["psms_raw"]]),
#'               "Master.Protein.Accessions",
#'               "Input")
#'
#' @export
message_parse <- function(x, column, note) {
  message(sprintf("%s features found from %s master proteins => %s",
                  nrow(x), length(unique(x[[column]])), note))
}
