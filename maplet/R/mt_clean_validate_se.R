#' Validate SE object is maplet-compatible
#'
#' Validate an externally created SE object is maplet-compatible. If not, make compatible.
#'
#' @param D \code{SummarizedExperiment} input.
#'
#' @return maplet-compatible SE
#'
#' @examples
#' # example of how to run function
#' \dontrun{... %>%
#'   mt_clean_validate_se() %>%
#' ...}
#'
#' @author KC
#'
#' @export
mt_clean_validate_se <- function(D){

  # check rowData contains a column called "name"
  feat_col_names <- D %>% rowData() %>% colnames()
  if("name" %in% feat_col_names == F) stop("rowData must contain a column called 'name'.")

  # check rowData column names are R Valid
  valid_feat_names <- feat_col_names %>% make.names()
  m <- match(valid_feat_names, feat_col_names) %>% is.na()
  if(any(m)) stop(glue::glue("The following rowData column names are not valid: {feat_col_names[m]}"))

  # add rownames if missing
  #   missing rownames creates issues with functions such as stats_univ_lm
  if(is.null(rownames(D))) rownames(D) <- paste0("x", seq(1, nrow(D)))

  # check colData column names are R valid
  samp_col_names <- D %>% colData() %>% colnames()
  valid_samp_names <- samp_col_names %>% make.names()
  m <- match(valid_samp_names, samp_col_names) %>% is.na()
  if(any(m)) stop(glue::glue("The following colData column names are not valid: {samp_col_names[m]}"))

  # check metadata empty and initialize entry
  if(metadata(D) %>% length() != 0) stop("SE object metadata must be empty.")
  metadata(D) <- list(sessionInfo=utils::sessionInfo())

  # add status information
  funargs <- mti_funargs()
  D %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = "External SE object validated maplet-compatible."
    )

  # return
  D

}
