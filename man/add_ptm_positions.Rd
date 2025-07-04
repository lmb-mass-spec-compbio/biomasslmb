% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ptm.R
\name{add_ptm_positions}
\alias{add_ptm_positions}
\title{Add rowData columns with details of PTMs positions}
\usage{
add_ptm_positions(
  obj,
  proteome_fasta,
  digest_enzyme = "trypsin-simple",
  missed_cleavages = c(0, 1, 2),
  master_protein_col = "Leading.razor.protein",
  sequence_col = "Sequence"
)
}
\arguments{
\item{obj}{\code{SummarizedExperiment}. Proteomics dataset}

\item{proteome_fasta}{\code{character} filepath to fasta with protein sequences
digest_enzyme = "trypsin-simple",}

\item{digest_enzyme}{\code{character}. Enzyme used. See \code{?cleaver::cleave}}

\item{missed_cleavages}{\code{numeric}. Vector of allowed number of missed cleavages}

\item{master_protein_col}{\code{character}. Name of column containing master proteins}

\item{sequence_col}{\code{character}. Name of column containing peptide sequences}
}
\value{
Returns a \code{SummarizedExperiment} with an additional column in the RowData describing the position of the PTMs with respect to the protein
}
\description{
This function takes a \code{SummarizedExperiment} object ...
}
