#' Restrict imputed values to specific conditions
#'
#' @description This function takes a `Qfeatures` object with unimputed data and
#' imputed data in separate assays and creates a new assay where the imputed values are
#' only used in specified circumstances. Note that this restriction occurs post imputation,
#' which may not be suitable for imputation methods where the imputed values are not independently derived
#'
#' @details
#' The `data.frame` provided to `use_imputed_df` should specify the conditions
#' under which imputation should be performed, with respect to the experimental
#' conditions and the number of missing values. `data.frame` must contain a column
#' called `n_finite` and all other columns must be columns in `colData(obj)`
#'
#' For example, in an enrichment vs
#' control experiment, to only input in control samples where zero or one replicate
#' are quantified, this should be specified thusly
#'
#' condition  n_finite
#' control  0
#' control  1
#'
#' @param obj `QFeatures`. Proteomics dataset
#' @param i_unimputed `string`. Index for the SummarizedExperiment with data without imputation
#' @param i_imputed `string`. Index for the SummarizedExperiment with data with imputation
#' @param i_restricted_imputed `string`. Index for the output assay with restricted imputation

#' @param use_imputed_df `data.frame` with conditions where imputation should be used. see @details for more information
#' @param verbose `logical` Describe how many missing values in input and output assays
#' @return Returns a `QFeatures` with restricted imputation with specified assay name
#' @export
restrict_imputation <- function(obj, i_unimputed, i_imputed,
                                i_restricted_imputed,
                                use_imputed_df,
                                verbose=TRUE){

  check_q(obj)

  non_imputed_data <- assay(obj[[i_unimputed]])
  imputed_data <- assay(obj[[i_imputed]])

  group_cols <- colnames(use_imputed_df) %>%
    setdiff(c('n_finite', 'use_imputed'))

  missing_variables <- setdiff(group_cols, colnames(colData(obj)))
  if(!length(missing_variables)==0){
    stop(sprintf(
      'Not all grouping columns are in the colData. Missing: %s',
      paste(missing_variables, collapse = ', ')))
  }

  use_imputed_df <- use_imputed_df %>%
    mutate(use_imputed=TRUE)

  use_imputed <- non_imputed_data %>% data.frame() %>%
    tibble::rownames_to_column('rowname') %>%
    pivot_longer(cols=-rowname) %>%
    mutate(name=remove_x(name)) %>%
    merge(colData(obj), by.x='name', by.y='row.names') %>%
    data.frame() %>%
    group_by(rowname, across(all_of(group_cols))) %>%
    mutate(n_finite=sum(is.finite(value))) %>%
    merge(use_imputed_df, by=c(group_cols, 'n_finite'), all.x=TRUE) %>%
    mutate(use_imputed=replace_na(use_imputed, FALSE)) %>%
    ungroup() %>%
    select(rowname, name, use_imputed) %>%
    pivot_wider(names_from=name, values_from=use_imputed) %>%
    tibble::column_to_rownames('rowname')


  # Ensure all objects have the rownames and colnames in the same order
  use_imputed <- use_imputed[rownames(non_imputed_data), colnames(non_imputed_data)]
  imputed_data <- imputed_data[rownames(non_imputed_data), colnames(non_imputed_data)]

  # create a new quantification dataset with a restricted set of imputation used
  imputed_restricted_data <- non_imputed_data
  imputed_restricted_data[as.matrix(use_imputed)] <- imputed_data[as.matrix(use_imputed)]

  if(verbose){
    # Checking how many missing values are in each dataset
    message(sprintf('Fraction NA in imputed data = %s', round(mean(is.na(imputed_data)), 4)))
    message(sprintf('Fraction NA in non-imputed data = %s', round(mean(is.na(non_imputed_data)), 4)))
    message(sprintf('Fraction NA in restricted-imputed data = %s', round(mean(is.na(imputed_restricted_data)), 4)))
  }

  # Create a new experiment in the QFeatures object with the restricted imputation
  obj[[i_restricted_imputed]] <- obj[[i_unimputed]]

  assay(obj[[i_restricted_imputed]]) <- imputed_restricted_data

  return(obj)
}

#' Create long format data with column defining if quantification value is imputed
#'
#' @description This function takes a `Qfeatures` object with unimputed data and
#' imputed data in separate assays and creates a long format data.frame with a column
#' called `is.imputed` which describes whether the quantification value was imputed

#' @param obj `QFeatures`. Proteomics dataset
#' @param i_unimputed `string`. Index for the SummarizedExperiment with data without imputation
#' @param i_imputed `string`. Index for the SummarizedExperiment with data with imputation
#' @param column_variables `character vector`. Variables in `colData(obj)` to include
#' @param row_variables `character vector`. Variables in `rowData(obj)` to include
#' @return Returns a `data.frame` with quantification data in long form and including
#' a `is.imputed` column
#' @export
create_long_form_imputed_data <- function(obj,
                                          i_unimputed,
                                          i_imputed,
                                          column_variables,
                                          row_variables){

  check_q(obj)

  long_format_protein <- longFormat(
    obj[,,c(i_imputed, i_unimputed)],
    colvars=column_variables,
    rowvars=row_variables) %>%
    data.frame()

  long_format_protein_imputed <- long_format_protein %>%
    filter(assay==i_imputed) %>%
    select(-assay)

  long_format_protein_not_imputed <- long_format_protein %>%
    filter(assay==i_unimputed) %>%
    select(-assay)

  long_format_protein_imputed_annot <- merge(long_format_protein_imputed,
                                             long_format_protein_not_imputed,
                                             by=setdiff(colnames(long_format_protein_imputed), c('value', 'assay'))) %>%
    mutate(is.imputed=is.na(value.y) & !is.na(value.x)) %>%
    mutate(value=value.x) %>%
    select(-c(value.x, value.y))

  return(long_format_protein_imputed_annot)

}
