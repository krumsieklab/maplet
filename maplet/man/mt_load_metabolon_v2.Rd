% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mt_load_metabolon_v2.R
\name{mt_load_metabolon_v2}
\alias{mt_load_metabolon_v2}
\title{Load Metabolon-format data (new format)}
\usage{
mt_load_metabolon_v2(D, file, data_sheet, met_sheet, samp_sheet)
}
\arguments{
\item{D}{\code{SummarizedExperiment} input. Missing if first step in pipeline.}

\item{file}{Name of input Excel file.}

\item{data_sheet}{Name of sheet with raw data.}

\item{met_sheet}{Name of sheet with metabolite info.}

\item{samp_sheet}{Name of sheet with sample info.}
}
\value{
If first step in pipeline, creates \code{SummarizedExperiment} object. Populates empty assay, colData, and rowData data frames.
}
\description{
For Metabolon-format version used after 2020-10-29. Loads data from a NEW Metabolon-format Excel file. Needs to be in the original
"Client Data Table" new format that they deliver.
}
\examples{
\dontrun{D <-
  # load data
  mt_load_metabolon_v2(file=system.file("extdata", "example_data/sampledata.xlsx", package = "maplet"),
                       data_sheet="Raw Data",
                       met_sheet="Chemical Annotation",
                       samp_sheet="Sample Meta Data") \%>\%
  ...}

}
\author{
RB
}
