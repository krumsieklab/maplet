#' Normalize by an external sample annotation
#'
#' Normalize by a column in colData, e.g. sample weight, DNA content, protein content (BRADFORD) etc.
#'
#' @param D \code{SummarizedExperiment} input
#' @param col_name Numeric-value column in colData (samples) to normalize by.
#'
#' @return assay: Normalized version.
#'
#' @examples
#' \dontrun{#' # in the context of a SE pipeline
#' ... %>% mt_pre_norm_external(col_name='DNA') %>% ...    # normalize by values in column DNA
#' }
#'
#' @author JK
#'
#' @export
mt_pre_norm_external = function(D, col_name) {

  # validate and extract arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  X = t(assay(D))
  if (any(unlist(X)[!is.na(unlist(X))]<0)) stop("Matrix contains negative values. Did you input logged data?")

  # get variable to normalize by
  if (!(col_name %in% colnames(colData(D)))) stop(sprintf("'%s' not found in sample annotations.", col_name))
  vc = colData(D)[[col_name]]
  if (!is.numeric(vc)) stop(sprintf("'%s' is not numeric.", col_name))

  # run normalization
  Y = t(sapply(1:dim(X)[1], function(i)unlist(X[i,]/vc[i])))
  rownames(Y) = rownames(X)

  # add status information
  funargs <- mti_funargs()
  metadata(D)$results %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = glue::glue("normalized by '{col_name}'")
    )

  # return
  assay(D) = t(Y)
  D

}
