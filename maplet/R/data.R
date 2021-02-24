#' A toy SE object with duplicated samples
#'
#' A \code{SummarizedExperiment} object with 100 samples and 100 metabolites. 5 out of the 100 samples are duplicated. Created with randomly
#' generated data.
#'
#' @format A \code{SummarizedExperiment} object with the following fields:
#' \describe{
#'   \item{assay}{Dataframe with 100 rows and 100 columns}
#'   \item{colData}{Dataframe with 100 rows and 2 columns}
#'   \item{rowData}{Dataframe with 100 rows and 0 columns}
#'   \item{metadata}{empty list}
#' }
"D_ex1"


#' A toy SE object with duplicated metabolites
#'
#' A \code{SummarizedExperiment} object with 50 samples and 50 metabolites. 5 out of the 100 samples are duplicated. The data for this object
#' is subset from the simulated_data.xlsx file.
#'
#' @format A \code{SummarizedExperiment} object with the following fields:
#' \describe{
#'   \item{assay}{Dataframe with 50 rows and 50 columns}
#'   \item{colData}{Dataframe with 50 rows and 6 columns}
#'   \item{rowData}{Dataframe with 50 rows and 15 columns}
#'   \item{metadata}{empty list}
"D_ex2"
