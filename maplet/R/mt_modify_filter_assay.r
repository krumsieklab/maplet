#' filter a maplet object using assay() data
#'
#' Filters assay data according to an expression.
#'
#' @param D \code{SummarizedExperiment} input.
#' @param filter Logical expression that selects which features to keep. Can use features from assay() data.
#'
#' @examples
#' \dontrun{... %>% mt_modify_filter_assay(filter = Glucose > 1) %>% ...}
#'
#' @author JCB
#'
#' @export


mt_modify_filter_assay <- function(D, filter) {
  # 1) Check if object is a summarized experiment
  stopifnot("SummarizedExperiment" %in% class(D))
  if (missing(filter)) stop("'filter' expression can't be empty.")

  # 2) Capture the filter expression
  filter_q <- dplyr::enquo(filter)

  # 3) Convert assay data to data frame with samples as rows
  #    (assuming each column of assay(D) is a sample, each row is a feature)
  assay_df <- assay(D) %>%
    t() %>%                           # transpose so that each sample is a row
    as.data.frame() %>%
    tibble::rownames_to_column("colnames")  # store sample (column) names in "colnames" col

  # 4) Apply filter using dplyr
  assay_df_filtered <- assay_df %>% dplyr::filter(!!filter_q)

  # 5) Identify which columns are kept
  included_cols <- assay_df_filtered$colnames

  # 6) Subset the SummarizedExperiment object
  #    We'll do so by subsetting columns (samples)
  D_new <- D[, included_cols]

  # 7) Identify which columns were excluded
  excluded <- setdiff(colnames(D), included_cols)

  # 8) Log the action (following maplet's pattern)
  funargs <- maplet:::mti_funargs()
  D_new %<>% maplet:::mti_generate_result(
    funargs = funargs,
    logtxt = sprintf("Filtered assay with expression: %s", as.character(filter_q)),
    output = excluded
  )

  # 9) Return the updated SummarizedExperiment
  D_new
}