#' Generate p-value histogram
#'
#' Generates a p-value histogram, either for a given list of statistical results or for all statistical results.
#'
#' @param D \code{SummarizedExperiment} input.
#' @param stat_list Name of the statistical results. If not given, will generate histogram for all. Default: NULL.
#'
#' @return $result$output: plot, p-value histogram
#'
#' @examples
#' \dontrun{... %>% mt_plots_pval_hist() %>% ...  # for all
#' ... %>% mt_plots_pval_hist(stat_list='comp') %>% ...  # for one
#' }
#'
#' @author JK
#'
#' @import ggplot2
#'
#' @export
mt_plots_pval_hist <- function(D, stat_list=NULL) {

  # validate argument
  stopifnot("SummarizedExperiment" %in% class(D))

  # if no stat_list given -> get all
  if (is.null(stat_list)) stat_list <- D %>% mtm_res_get_entries("stats") %>% purrr::map("output") %>% purrr::map("name") %>% unlist()

  # loop over stat_list
  plots <- lapply(stat_list, function(statname){
    # breaks
    st <- D %>% mtm_get_stat_by_name(statname)
    breaks <- pretty(range(st$p.value), n = grDevices::nclass.FD(st$p.value), min.n = 1)
    bwidth <- breaks[2]-breaks[1]
    # histogram
    p <- st %>%
        ggplot(aes(x=p.value)) +
        geom_histogram(binwidth=bwidth,fill="white",colour="black") +
        xlim(0,1) +
        ggtitle(glue::glue("'{statname}' p-values"))
    p
    })




  # add status information & plot
  funargs <- mti_funargs()
  D %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = glue::glue("P-value histograms for {paste0(stat_list, collapse=', ')}"),
      output = plots
    )

  # return
  D

}


