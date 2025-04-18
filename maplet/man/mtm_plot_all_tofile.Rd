% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mtm_access.R
\name{mtm_plot_all_tofile}
\alias{mtm_plot_all_tofile}
\title{Plot all plots from a pipeline}
\usage{
mtm_plot_all_tofile(input, dev = pdf, ...)
}
\arguments{
\item{input}{List of plots or SummarizedExperiment}

\item{dev}{Device to plot to (default: PDF)}

\item{...}{Further paramaters to be passed to dev() function.}
}
\value{
None
}
\description{
Opens a device (default: PDF), plots all plots, closes device.
Works either on list of plots, or on SE
}
