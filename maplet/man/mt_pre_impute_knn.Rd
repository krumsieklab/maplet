% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mt_pre_impute_knn.R
\name{mt_pre_impute_knn}
\alias{mt_pre_impute_knn}
\title{Impute using KNN method}
\usage{
mt_pre_impute_knn(
  D,
  method = "knn.obs.euc.sel",
  k = 10,
  verbose = F,
  use_multicore = F,
  n_cores = 5
)
}
\arguments{
\item{D}{\code{SummarizedExperiment} input.}

\item{method}{Specific KNN method to use. Default: "knn.obs.euc.sel" (recommendation: do not touch).}

\item{k}{Number of nearest neighbors to consider. Recommend not changing this parameter. Default: 10 (recommendation: do not touch).}

\item{verbose}{Output intermediate steps? Default: F.}

\item{use_multicore}{Use multicore KNN imputation method? Default: F.}

\item{n_cores}{If use_multicore==T, number of cores to use for imputation. Default: 5.}
}
\value{
assay: Imputed data.
}
\description{
Default settings are from the winning approach of our paper.
https://link.springer.com/article/10.1007%2Fs11306-018-1420-2

Multi-core functionality is available by setting the flag use_multicore=T.
}
\examples{
\dontrun{# in the context of a SE pipeline
... \%>\% mt_pre_impute_knn() \%>\% ...    # standard call
}

}
\author{
JK, PG
}
