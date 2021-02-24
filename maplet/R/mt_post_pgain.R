#' Compute p-gain from feature ratio test
#'
#' @description
#' This is a specialized post-analysis function for ratio-based datasets (see mt_modify_ratios). A "p-gain" is defined
#' as the factor of p-value decrease of a ratio, compared to the better of the p-values of the two single features. It
#' quantifies the amount of signal gained by working with the ratio.
#'
#' @description
#' For example, if feature A's association p-value is p=1e-5, feature B has p=1e-7, and the ratio has p=1e-9,
#' then pgain=min(1e-5, 1e-7)/1e-9 = 1e2 = 100.
#'
#' @description
#' See also: \href{https://pubmed.ncbi.nlm.nih.gov/22672667/}{https://pubmed.ncbi.nlm.nih.gov/22672667/}.
#'
#' @param D \code{SummarizedExperiment} input.
#' @param stat_name Name of the statistical comparison.
#' @param p_col Name of p-value column to compute p-gain from. Default: p.value.
#'
#' @return $results[[stat_name]]$output: pgain column added to statistical table.
#'
#' @examples
#' \dontrun{# add p-gains to the result table of the statistical comparison called "comparison1"
#' ... %>%
#'  mt_post_pgain(stat_name="comparison1") %>% ...
#'  }
#'
#' @author JZ
#'
#' @export
mt_post_pgain <- function(D, stat_name, p_col = p.value){

    p_col <- dplyr::enquo(p_col)

    ## are these ratio results?
    if (length(mtm_res_get_entries(D, c("modify","ratios"))) != 1)
        stop("must supply ratios to calculate p-gains")

    ## stat
    if(missing(stat_name))
        stop("stat_name must be given")

    ## find entry
    stat_id <- metadata(D)$results %>%
                         purrr::map_lgl(~"stats" %in% .x$fun && .x$output$name == stat_name) %>%
                         which()
    if(length(stat_id) == 0)
        stop("stat element with name ", stat_name, " does not exist")
    if(length(stat_id)  > 1)
        stop("there are multiple stat elements with name ", stat_name)

    ## DO CORRECTION
    rd <- rowData(D) %>%
        as.data.frame() %>%
        dplyr::mutate(var = rownames(D)) %>%
        dplyr::select(var, m1, m2)
    res <- metadata(D)$results[[ stat_id ]]$output$table %>%
                     dplyr::left_join(rd, by = "var")
    res <- res %>%
        dplyr::left_join(res %>% dplyr::filter(is.na(m2)) %>% dplyr::transmute(m1 = var, p_single_1 = p.value), by = "m1") %>%
        dplyr::left_join(res %>% dplyr::filter(is.na(m2)) %>% dplyr::transmute(m2 = var, p_single_2 = p.value), by = "m2") %>%
        dplyr::mutate(pgain = pmin(p_single_1, p_single_2) / p.value) %>%
        dplyr::select(-m1, -m2, -p_single_1, -p_single_2)

    ## update results table
    metadata(D)$results[[ stat_id ]]$output$table <- res
    metadata(D)$results[[ stat_id ]]$output$logtxt %<>% stringr::str_c(., ", pgain")

    ## add status information & plot
    funargs <- mti_funargs()
    metadata(D)$results %<>%
                  mti_generate_result(
                      funargs = funargs,
                      logtxt = sprintf("Calculate p-gains for '%s'", stat_name)
                  )
    ## RETURN
    D

}
