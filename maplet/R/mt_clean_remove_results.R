#' Remove metadata results from \code{SummarizedExperiment} object
#'
#' @description
#' Removes either the entire result list or all plots from the metadata of a \code{SummarizedExperiment}.
#' Can be used to obtain a more lightweight object for further processing.
#'
#' @description
#' remove=="all" will only leave a single result in the object, the one from this function.
#'
#' @description
#' remove=="plots" will set all ggplot objects to NULL, still allowing HTML reports to be generated.
#'
#' @param D \code{SummarizedExperiment} input.
#' @param remove What to remove. Has to be "all" or "plots". Default: "all".
#'
#' @return $results: Either all plots removed or all results entries removed.
#'
#' @examples
#' \dontrun{# at the end of a pipeline
#' ... %>% mt_clean_remove_results(remove="plots")
#' ...}
#'
#' @author JK
#'
#' @export
mt_clean_remove_results <- function(D, remove="all") {
  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  stopifnot(remove %in% c("all", "plots"))

  # delete all plots?
  if (remove=="plots") {
    # indices to plot results
    inds <- D %>% metadata() %>% .$results %>% purrr::map("fun") %>% purrr::map(~"plots" %in% .) %>% unlist() %>% which()
    # delete plots
    for (ind in inds) {
      metadata(D)$results[[ind]]['output'] <- list(NULL) # [] and list() needs to be used, otherwise $output is deleted
    }
  }

  # add status information & plot
  funargs <- maplet:::mti_funargs()
  D %<>% 
    maplet:::mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("Removed metadata results: %s", remove)
    )

  # remove everything but the last entry?
  if (remove=="all") {
    l <- D %>% metadata() %>% .$results %>% length()
    metadata(D)$results <-  metadata(D)$results[l]
  }

  # return
  D


}




