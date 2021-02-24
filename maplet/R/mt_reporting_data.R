#' Output information about dataset to log
#'
#' Leaves a log entry containing the number of samples, features and annotation fields at the current stage of the pipeline.
#'
#' @param D \code{SummarizedExperiment} input.
#'
#' @return Does not change the \code{SummarizedExperiment} object.
#'
#' @examples
#' \dontrun{...  %>%
#'   mt_reporting_data() %>% ...}
#'
#' @author JK
#'
#' @export
mt_reporting_data <- function(D) {

  # validate argument
  stopifnot("SummarizedExperiment" %in% class(D))

  # toc
  logtxt <- sprintf("Dataset info: %d samples, %d features; %d sample annotation fields, %d feature annotation fields",
                    ncol(D), nrow(D), ncol(colData(D)), ncol(rowData(D)))

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
