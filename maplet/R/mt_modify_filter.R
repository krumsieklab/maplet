#' Filter features
#'
#' Filters features according to an expression. Expression can access entries of rowData.
#'
#' @param D \code{SummarizedExperiment} input.
#' @param filter Logical expression. Can use columns from rowData.
#'
#' @return assay: Filtered features (rows) removed.
#' @return rowData: Filtered features removed.
#'
#' @examples
#' \dontrun{... %>% mt_modify_filter_features(filter=SUPER_PATHWAY=="Nucleotide") %>% ...}
#'
#' @author JK
#'
#' @export
mt_modify_filter_features <- function(D, filter){

    stopifnot("SummarizedExperiment" %in% class(D))
    if(missing(filter)) stop("feature filter can't be empty")

    ## APPLY FILTER TO ROW DATA
    filter_q <- dplyr::enquo(filter)
    rd <- rowData(D) %>%
        as.data.frame() %>%
        dplyr::mutate(rownames = rownames(D)) %>%
        dplyr::filter(!!filter_q)

    ## SUBSET FEATURES
    excluded <- rownames(D)[ !(rownames(D) %in% rd$rownames) ]
    D <- D[rd$rownames, ]

    ## add status information & plot
    funargs <- mti_funargs()
    D %<>% 
                  mti_generate_result(
                      funargs = funargs,
                      logtxt = sprintf("Filter features: %s",  as.character(filter_q)),
                      output = excluded
                  )
    ## return
    D
}

#' Filter samples
#'
#' Filters samples according to an expression. Expression can access entries of colData.
#'
#' @param D \code{SummarizedExperiment} input.
#' @param filter Logical expression. Can use columns from colData.
#'
#' @return assay: Filtered samples (columns) removed.
#' @return colData: Filtered samples removed.
#'
#' @examples
#' # filter to two specific groups of samples
#' \dontrun{... %>% mt_modify_filter_samples(filter = GROUP %in% c("FL","ctrl")) %>% ...}
#'
#' @author JK
#'
#' @export
mt_modify_filter_samples <- function(D, filter){

    stopifnot("SummarizedExperiment" %in% class(D))
    if(missing(filter)) stop("sample filter can't be empty")

    ## APPLY FILTER TO ROW DATA
    filter_q <- dplyr::enquo(filter)
    cd <- colData(D) %>%
        as.data.frame() %>%
        tibble::rownames_to_column("colnames") %>%
        dplyr::filter(!!filter_q)

    ## SUBSET SAMPLES
    cnames <- colData(D) %>% as.data.frame() %>% tibble::rownames_to_column("colnames") %>% .$colnames
    D <- D[, cnames %in% cd$colnames]

    ## add status information & plot
    funargs <- mti_funargs()
    D %<>% 
                  mti_generate_result(
                      funargs = funargs,
                      logtxt = sprintf("Filter samples: %s",  as.character(filter_q)),
                      output = list(kept=cnames %in% cd$colnames)
                  )
    ## return
    D
}
