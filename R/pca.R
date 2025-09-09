#' Create a Principle Component plot from the feature quantification
#'
#' @description A PCA visualisation of feature quantifications in each sample can
#' allow one to see how the experimental conditions relate to the sources of variance (principle components).
#' This function plots a PCA, with the option to colour and/or shape the points by experimental conditions.
#' The percentage values indicated on the axes are the variance explained by the PCs.
#'
#' @param obj `QFeatures`. Proteomics dataset
#' @param i `string`. Index for the SummarizedExperiment you wish to plot
#' @param allowing_missing `logical`. If TRUE, will use pcaMethods::pca to allow for missing values. If FALSE (default), will use stats::prcomp and remove any features with missing values
#' @param colour_by `string`. ColData column to colour points by
#' @param shape_by `string`. ColData column to shape points by
#' @param x `numeric`. Principle component to plot on x-axis
#' @param y `numeric`. Principle component to plot on x-axis
#' @return Returns a _ggplot_ object.
#' @export
plot_pca <- function(obj,
                     i,
                     allowing_missing=FALSE,
                     colour_by=NULL,
                     shape_by=NULL,
                     x = 1,
                     y = 2,
                     ...){


  check_q(obj)

  check_se_exists(obj, i)
  # Switch between prcomp and use pca from pcaMethods to handle cases with missing values

  if(allowing_missing){
    if(!requireNamespace("pcaMethods", quietly = TRUE)) {
      stop("pcaMethods package needed for this function to work. Please install it.",
           call. = FALSE)
    }
    pca <- pcaMethods::pca(t(assay(obj[[i]])), nPcs = ncol(assay(obj[[i]])), ...)
    proj <- pca@scores
    var_explained <- 100*pca@R2
  } else{
    pca <- prcomp(t(assay(filterNA(obj[[i]], pNA=0))), ...)
    var_explained <- 100*pca$sdev^2/sum(pca$sdev^2)
    proj <- pca$x
  }


  pca_x = sprintf('PC%s', x)
  pca_y = sprintf('PC%s', y)

  p <- proj %>%
    merge(data.frame(colData(obj)), by='row.names') %>%
    ggplot(aes(!!sym(pca_x), !!sym(pca_y))) +

    scale_shape_manual(values=21:25) +
    xlab(sprintf('%s (%s %%)', pca_x, round(var_explained[x], 2))) +
    ylab(sprintf('%s (%s %%)', pca_y, round(var_explained[y], 2))) +
    theme_biomasslmb()

  if(!is.null(colour_by)){

    check_colData_col(obj, colour_by)

    p <- p + aes(fill=!!sym(colour_by)) +
      guides(fill = guide_legend(override.aes = list(shape = 21)))
  }

  if(!is.null(shape_by)){

    check_colData_col(obj, shape_by)

    p <- p + geom_point(size=5, colour='grey20', stroke=0.5) + aes(shape=!!sym(shape_by))
  } else{
    p <- p + geom_point(size=5, colour='grey20', stroke=0.5, pch=21)
  }


  return(p)
}



