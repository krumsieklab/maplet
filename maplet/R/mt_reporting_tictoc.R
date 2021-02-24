#' Timing functionality
#'
#' Call mt_reporting_tic to start timing anywhere in pipeline, then mt_reporting_toc to show the time elapsed in the status log.
#'
#' @param D \code{SummarizedExperiment} input.
#'
#' @return Does not change the \code{SummarizedExperiment} object.
#'
#' @examples
#' \dontrun{%>%
#'   mt_reporting_tic() %>%
#'   ... %>% ... # pipeline steps
#'   ... %>% ... # pipeline steps
#'   ... %>% ... # pipeline steps
#'   mt_reporting_toc()}
#'
#' @author JK
#'
#' @import tictoc
#'
#' @export
mt_reporting_tic <- function(D) {

  # validate argument
  stopifnot("SummarizedExperiment" %in% class(D))

  # tic
  tic()
  # add status information & plot
  funargs <- mti_funargs()
  metadata(D)$results %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("tic")
    )

  # return
  D

}

#' Timing functionality
#'
#' Call mt_reporting_tic to start timing anywhere in pipeline, then mt_reporting_toc to show the time elapsed in the status log.
#'
#' @param D \code{SummarizedExperiment} input.
#'
#' @return Does not change the \code{SummarizedExperiment} object.
#'
#' @examples
#' \dontrun{%>%
#'   mt_reporting_tic() %>%
#'   ... %>% ... # pipeline steps
#'   ... %>% ... # pipeline steps
#'   ... %>% ... # pipeline steps
#'   mt_reporting_toc()}
#'
#' @author JK
#'
#' @import tictoc
#'
#' @export
mt_reporting_toc <- function(D) {

  # validate argument
  stopifnot("SummarizedExperiment" %in% class(D))

  # toc
  t <- toc(quiet=T)
  if (is.null(t)) {
    logtxt <- "toc without prior tic, no timing recorded"
  } else {
    logtxt = sprintf("toc, elapsed: %.2fs", t$toc-t$tic)
  }
  # add status information & plot
  funargs <- mti_funargs()
  metadata(D)$results %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = logtxt
    )

  # return
  D

}
