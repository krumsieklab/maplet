#' Exponentiate, base 2 by default
#'
#' Transform the the entire dataset base^x.
#'
#' @param D  \code{SummarizedExperiment} input.
#' @param base Operation: base^x for every data point. Default: 2.
#'
#' @return assay: Exponentiated data.
#'
#' @examples
#' \dontrun{# in the context of a SE pipeline
#' ... %>% mt_pre_trans_exp() %>% ...    # standard call, base 2
#' ... %>% mt_pre_trans_exp(base=10) %>% ...    # base 10
#' }
#'
#' @author JK
#'
#' @export
mt_pre_trans_exp <- function(D, base=2) {

  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  stopifnot(base%%1==0) # integer number

  # exp
  assay(D) = base^(assay(D))

  # add status information
  funargs <- mti_funargs()
  metadata(D)$results %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf('exp, base %d', base),
    )


  # return
  D

}
