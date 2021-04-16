#' Multiple testing correction
#'
#' Adjust output of statistical tests for multiple testing.
#'
#' @param D \code{SummarizedExperiment} input.
#' @param stat_name Name of the statistical comparison to adjust.
#' @param p_col Name of p-value column to adjust. Default: p.value.
#' @param method Which method to use for multiple testing. See \code{p.adjust} documentation. Default: "bonferroni".
#'
#' @return $results[[stat_name]]$output: p.adj column added to statistical table.
#'
#' @examples
#' \dontrun{# correct the statistical comparison called "Li's" using Benjamini-Hochberg
#' ... %>%
#'  mt_post_multtest(stat_name="Li's", method="BH") %>% ...
#'  }
#'
#' @author JK, JZ
#'
#' @export
mt_post_multtest <- function(D,
                             stat_name,
                             p_col = p.value,
                             method = "bonferroni"){

   p_col <- dplyr::enquo(p_col)

    ## stat
    if(missing(stat_name))
        stop("stat_name must be given")

    ## FIND ENTRY
    stat_id <- metadata(D)$results %>%
                         purrr::map_lgl(~"stats" %in% .x$fun && .x$output$name == stat_name) %>%
                         which()
    if(length(stat_id) == 0)
        stop("stat element with name ", stat_name, " does not exist")
    if(length(stat_id)  > 1)
        stop("there are multiple stat elements with name ", stat_name)

    ## DO CORRECTION
    metadata(D)$results[[stat_id]]$output$table %<>%
                  dplyr::mutate(p.adj = stats::p.adjust(!!p_col, method = method))


    ## add status information & plot
    funargs <- mti_funargs()
    D %<>% 
                  mti_generate_result(
                      funargs = funargs,
                      logtxt = sprintf("Multiple testing correction of '%s' using '%s'", stat_name, method)
                  )
    ## RETURN
    D

}
