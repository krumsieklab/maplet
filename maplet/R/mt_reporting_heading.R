#' Store heading for html report
#'
#' Store heading text and level to be used when calling \code{mt_reporting_html}.
#'
#' @param D  \code{SummarizedExperiment} input. Missing if first step in pipeline.
#' @param heading Heading text.
#' @param lvl Heading level. Can be used for nested outline structures. Default: 1.
#'
#' @return If first step in pipeline, creates \code{SummarizedExperiment} object.
#' @return $result$output: Stores info about heading.
#'
#' @examples
#' \dontrun{... %>%
#' # add first and second level of heading
#' mt_reporting_heading("Preprocessing") %>%
#' mt_reporting_heading("Part 1", lvl=2) %>%
#' ...}
#'
#' @author JK
#'
#' @export
mt_reporting_heading <- function(D, heading, lvl=1) {

  # if first step in pipeline, create SE
  if(missing(D)){
    # create an empty SummarizedExperiment object
    D <- SummarizedExperiment()
  }else{
    # validate argument
    stopifnot("SummarizedExperiment" %in% class(D))
  }

  # add status information & heading info
  funargs <- mti_funargs()
  metadata(D)$results %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = glue::glue("reporting heading, level {lvl}: {heading}"),
      output = list(lvl=lvl,title=heading)
    )

  # return
  D

}


