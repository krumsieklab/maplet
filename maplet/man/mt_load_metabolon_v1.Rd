% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mt_load_metabolon_v1.R
\name{mt_load_metabolon_v1}
\alias{mt_load_metabolon_v1}
\title{Load Metabolon-format data (old format)}
\usage{
mt_load_metabolon_v1(D, file, sheet, copy_nan_sheet = "")
}
\arguments{
\item{D}{\code{SummarizedExperiment} input. Missing if first step in pipeline.}

\item{file}{Name of input Excel file.}

\item{sheet}{Name of sheet.}

\item{copy_nan_sheet}{If given, which sheet to copy the NA pattern from. Default: '' (i.e. None).}
}
\value{
If first step in pipeline, creates \code{SummarizedExperiment} object. Populates empty assay, colData, and rowData data
frames.
}
\description{
For Metabolon-format version used prior to 2020-10-29. Loads data from a Metabolon-format Excel file. Needs to be in the original
"Client Data Table" format that they deliver.
}
\examples{
\dontrun{D <-
  # load data
  mt_load_metabolon_v1(file=system.file("extdata", "example_data/sampledata.xlsx", package = "maplet"),
                       sheet="OrigScale") \%>\%
  ...}

}
\author{
JK
}
