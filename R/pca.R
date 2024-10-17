#' Create a Principle Component plot from the feature quantification
#'
#' @description A PCA visualisation of feature quantifications in each sample can
#' allow one to see how the experimental conditions relate to the sources of variance (principle components).
#' This function plots a PCA, with the option to colour and/or shape the points by experimental conditions.
#' The percentage values indicated on the axes are the variance explained by the PCs.
#'
#' @param obj `QFeatures`. Proteomics dataset
#' @param i `string`. Index for the SummarizedExperiment you wish to plot
#' @param colour_by `string`. ColData column to colour points by
#' @param shape_by `string`. ColData column to shape points by
#' @param x `numeric`. Principle component to plot on x-axis
#' @param y `numeric`. Principle component to plot on x-axis
#' @return Returns a _ggplot_ object.
#' @export
plot_pca <- function(obj,
                     i,
                     colour_by=NULL,
                     shape_by=NULL,
                     x = 1,
                     y = 2){

  check_q(obj)

  check_se_exists(obj, i)

  pca <- prcomp(t(assay(filterNA(obj[[i]]))))

  var_explained <- 100*pca$sdev^2/sum(pca$sdev^2)

  pca_x = sprintf('PC%s', x)
  pca_y = sprintf('PC%s', y)

  p <- pca$x %>%
    merge(data.frame(colData(obj)), by='row.names') %>%
    ggplot(aes(!!sym(pca_x), !!sym(pca_y))) +
    geom_point(size=5, colour='grey20', stroke=0.5) +
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

    p <- p + aes(shape=!!sym(shape_by))
  }

  return(p)
}



