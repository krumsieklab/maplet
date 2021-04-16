#' Averages duplicate features
#'
#' Averages the values of features (rows) with the same values in specified rowData column.
#'
#' @param D \code{SummarizedExperiment} input.
#' @param group_col Name of rowData column (feature annotation) by which duplicates can be identified.
#'
#' @return assay: Duplicate features (rows) combined.
#'
#' @examples
#' \dontrun{... %>% mt_modify_avg_features(group_col = 'name') %>% ...}
#'
#' @author RB
#'
#' @import dplyr
#'
#' @export
mt_modify_avg_features <- function(D, group_col) {

  stopifnot("SummarizedExperiment" %in% class(D))
  if(missing(group_col))
    stop("group_col can't be empty")

  total_rows <- nrow(D)
  # TODO: check that rowdata for duplicates are the same, throw mti_logwarning if not

  # mean per feature by group_col column in rowData
  X <- dplyr::bind_cols(D %>% assay() %>% data.frame(), D %>% rowData() %>% data.frame() %>%
                          dplyr::select(!!as.name(group_col))) %>%
    dplyr::group_by(!!as.name(group_col)) %>%
    dplyr::summarise_at(dplyr::vars(-group_cols()), mean) %>% # mean per feature
    dplyr::ungroup()

  # identify the unique rows
  unique_rows <- D %>% rowData() %>% data.frame() %>% dplyr::group_by(!!as.name(group_col)) %>%
    group_rows() %>%  lapply(., FUN=function(x)x[1]) %>% unlist()

  # subset SE for unique rows
  D <- D[unique_rows, ]

  # replace assay data with means
  assay(D, withDimnames = F) <- X[match(rowData(D)[[group_col]], X[[group_col]]), ] %>% select(-!!as.name(group_col))

  ## add status information
  funargs <- mti_funargs()
  D %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf('%d duplicate features averaged', (total_rows-length(unique_rows)))
    )

  ##return
  D

}
