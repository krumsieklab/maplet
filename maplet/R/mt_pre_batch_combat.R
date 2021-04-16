#' ComBat batch correction
#'
#' Performs batch correction via COMBAT method from \code{sva} package. Function will check that data is logged. If data is not logged,
#' this function will log the data, run ComBat, and exponentiate the data before returning it into the assay data frame.
#'
#' @param D \code{SummarizedExperiment} input.
#' @param batch_col Sample annotation (colData) column name that contains batch assignment.
#'
#' @return assay: Batch-corrected version.
#'
#' @examples
#' \dontrun{... %>% mt_pre_batch_combat(batch_col="BATCH") %>% ...
#' }
#'
#' @author JK
#'
#' @export
mt_pre_batch_combat = function(D, batch_col) {

  # validate and extract arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  X = assay(D) # combat uses data probe s sample

  # crash if there are any missing values
  if (any(is.na(X))) {
    stop("ComBat batch correction cannot work with missing/NA values.")
  }

  # check if data is logged, if not log the data
  is_logged <- maplet:::mti_check_is_logged(D)
  if(!is_logged){
    X <- log(X, base = 2)
  }

  # get variable
  if (!(batch_col %in% colnames(colData(D)))) stop(sprintf("'%s' not found in sample annotations.", batch_col))
  b = colData(D)[[batch_col]]
  # ensure it's a factor or character vector
  if (!(is.character(b) || is.factor(b))) stop(sprintf("'%s' has to be character, factor, or numeric.", batch_col))
  b = as.factor(b)

  # run ComBat
  # invisible (capture.output(     # seems to be the only way to silence ComBat() output to stdout and stderr
  #   invisible( capture.output(
      Y <- sva::ComBat(
        dat=X,
        batch=b,
        mod=NULL,
        par.prior=TRUE,
        prior.plots=FALSE)


  # add status information
  funargs <- mti_funargs()
  D %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("ComBat batch correction, variable: '%s'",batch_col)
    )

  # if data was logged above, exponentiate it before returning
  if(!is_logged){
    Y <- 2^Y
  }

  # return
  assay(D) = Y
  D


}
