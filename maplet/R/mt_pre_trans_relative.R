#' Scale each sample relative to the mean of a given group of samples
#'
#' Usually used to make all concentrations relative to a control group, for example.
#' If data are logged, will use a 'minus' operation, if data are non-logged, will divide.
#' Whether data are logged is determine by whether they have been logged inside the MT pipeline.
#'
#' @param D \code{SummarizedExperiment} input.
#' @param ref_samples Expression filtering reference samples from colData.
#' @param is_logged Is the data logged?
#'
#' @return assay: relatively scaled data
#'
#' @examples
#' \dontrun{# normalize to control group, data not logged
#' ... %>% mt_pre_trans_relative(ref_samples = GROUP=="ctrl", is_logged=F) %>% ...
#' }
#'
#' @author JK
#'
#' @export
mt_pre_trans_relative <- function(D, ref_samples, is_logged) {

  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))

  X = t(assay(D))

  # check if data are logged
   if (!is_logged && (any(unlist(X)[!is.na(unlist(X))]<0))) stop("Data is not logged but contains negative values.")

  # find samples to normalize to
  sample_filter_q <- dplyr::enquo(ref_samples)
  cd <- colData(D) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("colnames") %>%
    dplyr::filter(!!sample_filter_q)
  # define samples to be used
  useref <- (colnames(D) %in% cd$colnames)

  # pick operation
  if (!is_logged) op <- `/`
  else op = `-`
  # calculate average vector across reference samples
  avg <- X[useref,] %>% apply(2, mean, na.rm=T)
  X %<>% apply(1, function(s)op(s,avg)) # comes out transposed the right way

  # save assay
  assay(D) = X

  # add status information
  funargs <- mti_funargs()
  metadata(D)$results %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = glue::glue('scaled samples relative to {sum(useref)} reference samples: {dplyr::enquo(ref_samples) %>% as.character()}')
    )

  # return
  D

}
