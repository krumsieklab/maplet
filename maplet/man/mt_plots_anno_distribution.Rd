% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mt_plots_anno_distribution.R
\name{mt_plots_anno_distribution}
\alias{mt_plots_anno_distribution}
\title{Annotation Distribution Plots}
\usage{
mt_plots_anno_distribution(
  D,
  anno_type,
  col_list,
  single = T,
  pairwise = F,
  cat_limit = 25
)
}
\arguments{
\item{D}{\code{SummarizedExperiment} input.}

\item{anno_type}{Must be either "samples" (colData) or "features" (rowData).}

\item{col_list}{List of column names from either colData or rowData. If missing, uses all available columns.}

\item{pairwise}{Create pairwise plots for annotations? Default: F.}

\item{cat_limit}{Limit the number of unique values allowed for categorical variables. Default: 25.}
}
\value{
$results$output: plot list, annotation bar and / or histogram plots
}
\description{
Create distribution plots for annotation columns. Produces bar plots for categorical annotations and histogram plots for numeric
annotations. Works on one annotation type at a time. If list of column names not provided, produces plots for all columns in
selected annotation data frame.
}
\examples{
ADD

}
\author{
KC
}
