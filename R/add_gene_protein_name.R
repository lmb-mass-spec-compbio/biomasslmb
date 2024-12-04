#' Extract the gene name and long-form protein name from the master protein descriptions column
#'
#' @description This function extracts the gene name and long-form protein name
#' from the master protein descriptions column in output from
#' Proteome Discoverer for DDA.
#'
#' It assumes the master protein description column is in the format
#' (.*) .* GN=(.*) PE.*, where the first group is the long format of the protein
#' name and the second group is the gene name. Where the gene name is not included,
#' an empty string is returned
#' An example of the expected input in the column is:
#' Serum albumin OS=Bos taurus GN=ALB PE=1 SV=4
#'
#'
#' ```
#' @param obj `SummarisedExperiment` containing output from Proteome Discoverer
#' @param master_protein_desc_col `string`. Name of column containing master
#' proteins descriptions.
#' @param gene_name_out_col `string`. Name of output column containing the gene names
#' @param long_protein_name_out_col `string`. Name of output column containing the long format
#' protein names

#' @return Returns a `QFeatures` with the filtered Proteome Discoverer output.
#' @export
add_gene_long_protein_name_pd <- function(obj,
                                          master_protein_desc_col='Master.Protein.Descriptions',
                                          gene_name_out_col='Master.Protein.gene.name',
                                          long_protein_name_out_col='Master.Protein.long.name'
         ){

  check_se(obj)

  rowData(obj)$Master.Protein.gene.name <-
    sapply(strsplit(rowData(obj)[[master_protein_desc_col]], ';'),
           function(x){
             if(any(grepl("GN=", x))){
               paste(gsub(".* GN=(.*) PE.*", "\\1", x[grepl("GN=", x)]), collapse = ";")
             } else ''
           })

  rowData(obj)$Master.Protein.long.name <-
    sapply(strsplit(rowData(obj)[[master_protein_desc_col]], ';'),
           function(x) paste(
             gsub(' OS=.*', '', x),
             collapse=';'))

  return(obj)
}
