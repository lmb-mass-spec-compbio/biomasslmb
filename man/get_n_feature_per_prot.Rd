% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/restrict_features_per_protein.R
\name{get_n_feature_per_prot}
\alias{get_n_feature_per_prot}
\title{Identify how many features (PSMs/Peptides) are quantified for each protein}
\usage{
get_n_feature_per_prot(obj, master_protein_col = "Master.Protein.Accessions")
}
\arguments{
\item{obj}{\code{SummarizedExperiment} with PSM or peptide-level quantification}

\item{master_protein_col}{\code{character} Column name for master protein}
}
\value{
\code{data.frame} detailing how many features are present for each protein in each sample
}
\description{
For summarisation of PSM or peptide to protein, we need a
minimum number of finite values per protein per sample. This function simply tallies how many we have.
}
