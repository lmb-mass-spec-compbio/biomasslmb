% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tmt_qc_plots.R
\name{plot_missing_SN_per_sample}
\alias{plot_missing_SN_per_sample}
\title{Plot the missing values vs signal:noise for each sample}
\usage{
plot_missing_SN_per_sample(obj, sn_column = "Average.Reporter.SN", bins = 20)
}
\arguments{
\item{obj}{\code{SummarizedExperiment} containing PSM level TMT intensities}

\item{sn_column}{\code{character} column name for Signal:noise values}

\item{bins}{\code{numeric} Number of bins to plot}
}
\value{
\code{ggplot} tile plot to show S:N vs # missing values for each sample
}
\description{
Missing values are more frequent with low signal:noise (S:N).
This function visualises this relationship for each sample to aid selection
of thresholds for minimal S:N filtering.
}
