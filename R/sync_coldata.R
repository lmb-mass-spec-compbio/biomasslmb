#' Copy the experimental design onto one or more assays
#'
#' @description `QFeatures::readQFeatures()` attaches `colData` at the
#' `QFeatures` level only, so the individual `SummarizedExperiment` assays
#' inside it start with an empty `colData`. Functions which plot or model a
#' single assay — `plot_quant()`, `plot_pca()`, `condition_miss_score()` and
#' anything using `limma` — need the design on the assay itself.
#'
#' `sync_coldata()` copies the object-level `colData` onto the named assays,
#' matching on column name so that an assay holding a subset of the samples
#' gets the rows belonging to it.
#'
#' Assays created by `QFeatures::aggregateFeatures()` and
#' `QFeatures::joinAssays()` also start without `colData`, so this is worth
#' calling after those as well as after reading the data in.
#'
#' @param obj `QFeatures`. Proteomics dataset
#' @param i `character` or `numeric`. Assays to copy the `colData` onto.
#'   Defaults to every assay in `obj`.
#'
#' @return `QFeatures` with `colData` attached to the named assays
#'
#' @examples
#' tmt_qf <- QFeatures::readQFeatures(assayData = psm_tmt_clock,
#'   colData = tmt_clock_design,
#'   quantCols = rownames(tmt_clock_design),
#'   name = "psms_raw")
#'
#' # the assay starts without the design attached
#' dim(SummarizedExperiment::colData(tmt_qf[["psms_raw"]]))
#'
#' tmt_qf <- sync_coldata(tmt_qf)
#'
#' dim(SummarizedExperiment::colData(tmt_qf[["psms_raw"]]))
#' @export
sync_coldata <- function(obj, i = names(obj)) {

  check_q(obj)

  object_coldata <- colData(obj)

  for (assay_name in i) {

    check_se_exists(obj, assay_name)

    assay_samples <- colnames(obj[[assay_name]])

    missing_samples <- setdiff(assay_samples, rownames(object_coldata))

    if (length(missing_samples) > 0) {
      stop(sprintf(
        paste("Assay `%s` has samples with no row in the object-level colData: %s.",
              "Check the sample names match the rownames of the colData."),
        assay_name, paste(missing_samples, collapse = ', ')))
    }

    colData(obj[[assay_name]]) <- object_coldata[assay_samples, , drop = FALSE]

  }

  obj
}
