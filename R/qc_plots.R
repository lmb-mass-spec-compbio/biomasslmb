#' Plot distributions for feature intensities per sample.
#'
#' @description Given a `SummarizedExperiment`, return a plot
#' of the feature quantifications per sample.
#'
#' @param obj `SummarizedExperiment`.
#' @param method `string`. Plot style. Choice of box, density or histogram plot.
#' @param log2transform `logical`. Should feature quantifications be log-transformed?
#' @param facet_by_sample `logical`. Facet the plot by sample.
#'
#' @return `ggplot` object.
#' @export
plot_quant <- function(obj,
                       method=c('box', 'density', 'histogram'),
                       log2transform=FALSE,
                       facet_by_sample=FALSE){
  # check method argument
  method = match.arg(method)

  check_se(obj)

  e_data <- data.frame(assay(obj))

  e_data[e_data==""] <- NA
  e_data <- e_data %>% gather(key='sample', value='intensity') %>%
    mutate(sample = factor(remove_x(sample), levels = remove_x(colnames(e_data))))

  if(log2transform){
    e_data$intensity <- log2(e_data$intensity)
    intensity_label <- 'Feature intensity (log2)'
  } else{
    intensity_label <- 'Feature intensity'
  }

  p <- ggplot(e_data) + theme_bw()

  if(method=='box'){
    p <- p +
      geom_boxplot(aes(sample, intensity)) +
      theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1)) +
      ylab(intensity_label) +
      xlab("")
  }
  else if(method=='density'){
    p <- p +
      stat_density(aes(intensity, col=sample), geom="line",position="identity") +
      xlab(intensity_label) +
      ylab("Density")
  }
  else if(method=='histogram'){
    p <- p +
      geom_histogram(aes(intensity), bins=100) +
      scale_y_continuous(expand = c(0, 0)) +
      xlab(intensity_label) +
      ylab("Count")
  }

  if(facet_by_sample){
    p <- p + facet_wrap(~sample)
  }
  return(p)
}
