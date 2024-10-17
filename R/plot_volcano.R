#' Volcano plot for differential abundance testing
#'
#' @description Statistical testing results for 'omics' differential abundance
#' are frequently displayed in a 'volcano' plot, where x = log fold-change
#' and y = -log10(p-value). This function generates a volcano plot from a `data.frame`
#' containing the statistical testing results
#'
#' @param obj `data.frame`. Statistical testing results
#' @param lf_col `string`. Column for log fold-change
#' @param p_col `string`. Column for p-value
#' @param sig_col `string`. Column for significance annotation (optional)
#' @param xlab `string`. Label for x axis
#' @param title `string`. Plot title
#'
#' @return Returns a _ggplot_ object.
#' @export
#' @importFrom naniar gg_miss_upset
plot_volcano <- function(obj,
                         lf_col='logFC',
                         p_col='adj.P.Val',
                         sig_col=NULL,
                         title=NULL,
                         xlab='log2 Fold Change'){

  p <- ggplot(obj, aes(x=!!sym(lf_col), y=-log10(!!sym(p_col)))) +
    geom_point(size=2, pch=21, colour='grey', stroke=0.2) +
    theme_biomasslmb() +
    labs(x=xlab, y='-log10(P-value)')

  if(!is.null(sig_col)){
    p <- p + aes(fill=!!sym(sig_col))
  }

  if(!is.null(title)){
    p <- p + ggtitle(title) +
      theme(plot.title=element_text(hjust=0.5))
  }

  return(p)
}
