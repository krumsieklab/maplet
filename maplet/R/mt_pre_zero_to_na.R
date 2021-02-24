#' Set zeros in dataset to NA
#'
#' Every 0 value in the dataset is replaced by an NA. Used for platforms that represent missing/sub-LOD values as zeros.
#'
#' @param D \code{SummarizedExperiment} input.
#'
#' @return assay: Zeros replaced by NA.
#'
#' @examples
#' \dontrun{# in the context of a SE pipeline
#' ... %>% mt_pre_zero_to_na() %>% ...
#' }
#'
#' @author JK
#'
#' @export
mt_pre_zero_to_na <- function(D) {

  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))

  # replace
  X <- assay(D)
  X[X==0] <- NA
  assay(D) <- X

  # add status information
  funargs <- mti_funargs()
  metadata(D)$results %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = "zeros replaced by NAs"
    )

  # return
  D

}
