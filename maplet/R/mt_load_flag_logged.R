#' Flag data as logged
#'
#' Required when loading already-logged data - signal to other functions (e.g. fold change calculation) that this is logged data.
#' Technically, function does not do anything, but leaves its call footprint in the pipeline to be found by others.
#'
#' @param D \code{SummarizedExperiment} input.
#'
#' @return Does not change the \code{SummarizedExperiment} object.
#'
#' @examples
#' \dontrun{# simply insert into pipeline
#' ... %>% mt_load_flag_logged() %>%
#' ...}
#'
#' @author JK
#'
#' @export
mt_load_flag_logged <- function(D){

  # validate argument
  stopifnot("SummarizedExperiment" %in% class(D))

  # add status information & plot
  funargs <- mti_funargs()
  metadata(D)$results %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = "Data flagged as logged"
    )

  # return
  D

}


