% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mt_load_ucd.R
\name{mt_load_ucd}
\alias{mt_load_ucd}
\title{Load UC Davis-format data}
\usage{
mt_load_ucd(D, file, sheet, zero_to_na = T, gen_valid_varnames = F)
}
\arguments{
\item{D}{\code{SummarizedExperiment} input. Missing if first step in pipeline.}

\item{file}{Name of input Excel file in UCD format.}

\item{sheet}{Name or number of sheet.}

\item{zero_to_na}{Replace zeros by NAs? Default: T.}

\item{gen_valid_varnames}{Enforce valid R variable names? Default: F.}
}
\value{
If first step in pipeline, creates \code{SummarizedExperiment} object. Populates empty assay, colData, and rowData data frames.
}
\description{
Loads data from a UC Davis-format Excel file.
}
\examples{
\dontrun{D <-
  # load data
  mt_load_ucd(file=system.file("extdata", "example_data/sampledata.xlsx", package = "maplet"),
              sheet="OrigScale") \%>\%
  ...}

}
\author{
JK
}
