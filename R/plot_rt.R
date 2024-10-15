#' Plot the retention time vs delta precursor mass
#'
#' @description It's useful to assess the relationship between retention time and delta precursor
#' mass as a quality control step. There should be little or no relationship
#'
#' @param obj `SummarisedExperiment` containing peptide-level output from Proteome Discoverer.
#' @param rt_col `string`. Name of column with retention time
#' @param delta_ppm_col `string`. Name of column with the Delta precursor mass
#' @return Returns a `ggplot` with the RT vs delta PPM
#' @export
plot_rt_vs_delta <- function(obj,
                             rt_col = 'RT.in.min.by.Search.Engine.Sequest.HT',
                             delta_ppm_col = 'Delta.M.in.ppm.by.Search.Engine.Sequest.HT'){

  check_se_peptide(obj)

  obj %>%
    rowData() %>%
    as_tibble() %>%
    ggplot(aes(x = !!sym(rt_col),
               y = !!sym(delta_ppm_col))) +
    geom_point(size = 0.5, shape = 4) +
    geom_hline(yintercept = 5, linetype = "dashed", color = "red") +
    geom_hline(yintercept = -5, linetype = "dashed", color = "red") +
    labs(x = "RT (min)", y = "Delta precursor mass (ppm)") +
    theme_biomasslmb(border=FALSE)
}

#' Plot the retention time distribution
#'
#' @description It's useful to assess the retention time distribution as a quality control step.
#' The peptides should be spread across the retention times, without any clear 'gaps' or very clear
#' trend towards one end of the gradient.
#'
#' @param obj `SummarisedExperiment` containing peptide-level output from Proteome Discoverer.
#' @param rt_col `string`. Name of column with retention time
#' @return Returns a `ggplot` with the RT vs delta PPM
#' @export
plot_rt_dist <- function(obj,
                         rt_col = 'RT.in.min.by.Search.Engine.Sequest.HT'){

  check_se_peptide(obj)

  obj %>%
    rowData() %>%
    as_tibble() %>%
    ggplot(aes(x = !!sym(rt_col))) +
    geom_histogram(binwidth = 1) +
    labs(x = "RT (min)", y = "Frequency") +
    theme_biomasslmb(border=FALSE)

}
