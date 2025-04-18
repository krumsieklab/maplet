% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mt_ml_repeat.R
\name{mt_ml_repeat}
\alias{mt_ml_repeat}
\title{Run Cross-Validation / Sampling Methods}
\usage{
mt_ml_repeat(
  D,
  ml_fun,
  ml_name,
  response_col,
  response_type = "binary",
  covar_cols,
  sampling_method = "cv",
  num_folds,
  rand_seed,
  mod_args = list(),
  pred_args = list(),
  pos_class
)
}
\arguments{
\item{D}{\code{SummarizedExperiment} input.}

\item{ml_fun}{The maplet-implemented algorithm to use. See mtml_algorithms for a list of implemented algorithms.}

\item{ml_name}{Name under which the results will be stored.}

\item{response_col}{Name of the colData column containing the data labels.}

\item{response_type}{Indicates the type of learning task. One of "binary" or
"continuous". Default: "binary".}

\item{covar_cols}{An optional vector of covariate columns taken from colData.}

\item{sampling_method}{Type of cross-validation / sampling method using. Must be one of: cv (add more later).
Default: cv.}

\item{num_folds}{For k-fold cross-validation, this represents the number of folds (k).}

\item{rand_seed}{Value to set random seed.}

\item{mod_args}{Named list of parameters passed to maplet-implemented machine learning algorithm. See algorithm
implementations for default parameters used. Default: empty list.}

\item{pred_args}{Named list of parameters passed to predict() for a specific machine learning function. For example,
predict.glmnet accepts the parameter s, which supplies the penalty parameter lambda at which prediction are
required. See details of each function in mtml_algorithms for parameter lists for each algorithm-specific
predict function. See mpalet-implemented ml functions for default parameters used. Default: empty list.}

\item{pos_class}{Required if response_type = "binary". Indicates which value in the response column
represents the positive class.}
}
\value{
result$output: List containing a list of results for each fold.
}
\description{
Performs one of the following cross-validation / sampling methods for a provided machine learning
algorithm: k-fold cross-validation, (TO-DO: extend to more algorithms later).
}
\author{
KC
}
