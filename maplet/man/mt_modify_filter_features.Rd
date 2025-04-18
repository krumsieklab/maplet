% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mt_modify_filter.R
\name{mt_modify_filter_features}
\alias{mt_modify_filter_features}
\title{Filter features}
\usage{
mt_modify_filter_features(D, filter)
}
\arguments{
\item{D}{\code{SummarizedExperiment} input.}

\item{filter}{Logical expression that selects which features to keep. Can use columns from rowData.}
}
\value{
assay: Filtered features (rows) removed.

rowData: Filtered features removed.
}
\description{
Filters features according to an expression. Expression can access entries of rowData. Note: the filtering
expression selects which features to KEEP.
}
\examples{
\dontrun{... \%>\% mt_modify_filter_features(filter=SUPER_PATHWAY=="Nucleotide") \%>\% ...}

}
\author{
JK
}
