#' Output information about statistical results
#'
#' Leaves a log entry containing the number of samples in each group that was used, and the number of hits according to a formula.
#'
#' @param D \code{SummarizedExperiment} input.
#' @param stat_name Name of the statistical comparison.
#' @param stat_filter Filter formula to display number of results, e.g. p.adj<0.2.
#'
#' @return Does not change the \code{SummarizedExperiment} object.
#'
#' @examples
#' \dontrun{...  %>%
#'   mt_reporting_stats(stat_name="comp1", stat_filter=p.adj<0.2) %>% ...}
#'
#' @author JK
#'
#' @export
mt_reporting_stats <- function(D, stat_name, stat_filter) {

  # validate argument
  stopifnot("SummarizedExperiment" %in% class(D))

  # trick: access argument so that a missing argument error is thrown from here instead of from inside mtm_get_stat_by_name
  stat_name
  # find statistical result (not just the table, but the entire $output block)
  allstats <- D %>% maplet::mtm_res_get_entries(c("stats"))
  statind <- allstats %>% purrr::map("output") %>% purrr::map("name") %>% purrr::map(~.==stat_name) %>% unlist() %>% which()
  if (length(statind)==0) stop(sprintf("comparison '%s' not found", stat_name))
  # retrieve actual structures
  output <- allstats %>% purrr::map("output") %>% .[[statind[1]]]

  # initialize logtxt
  logtxt <- sprintf("'%s' info\n\n", stat_name)
  # tabulate used samples with groups to output n's, if outcome variable was not numeric
  if (!is.numeric((D %>% colData() %>% .[[output$outcome]]))) {
    N <- table(D %>% colData() %>% .[[output$outcome]], output$samples.used)[,"TRUE"]
    N <- N[N>0] # filter out empty groups
    # toc
    logtxt <- sprintf("%sSamples:\n%s",  logtxt, paste(names(N),N,sep=": ",collapse=", "))
  }

  # filter?
  if (!missing(stat_filter)) {
    # extract formula
    s <- dplyr::enquo(stat_filter) %>% as.character()
    # remove all the "~" formula characters
    s <- gsub("~","", s)
    # number of hits
    hits <- output$table %>% dplyr::filter(!!dplyr::enquo(stat_filter)) %>% nrow()
    # add to log text
    logtxt <- sprintf("%s /// features: %s, %d of %d", logtxt, s, hits, output$table %>% nrow())
  }

  # add status information & plot
  funargs <- maplet:::mti_funargs()
  metadata(D)$results %<>%
    maplet:::mti_generate_result(
      funargs = funargs,
      logtxt = logtxt
    )

  # return
  D

}
