#' Log, base 2 by default
#'
#' Transform the entire dataset log_based(x), i.e. log(x)/log(base).
#'
#' @param D \code{SummarizedExperiment} input.
#' @param base Operation: log(x)/log(base) for every data point. Default: 2.
#'
#' @return assay: Logged data.
#'
#' @examples
#' \dontrun{# in the context of a SE pipeline
#' ... %>% mt_pre_trans_log() %>% ...    # standard call, base 2
#' ... %>% mt_pre_trans_log(base=10) %>% ...    # base 10
#' }
#'
#' @author JK
#'
#' @export
mt_pre_trans_log <- function(D, base=2) {

  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  if(base%%1!=0) warning(sprintf("Base not integer. Please double check. Base %f.", base))

  # log
  assay(D) = log(assay(D), base = base)

  # add status information
  funargs <- maplet:::mti_funargs()
  metadata(D)$results %<>%
    maplet:::mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("log%f", base)
    )

  # return
  D

}
