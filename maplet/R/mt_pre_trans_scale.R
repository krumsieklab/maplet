#' Scale data, mean 0 / sd 1 by default
#'
#' This function scales each feature by subtracting its mean and dividing by its standard deviation.
#' Both steps are optional. This function is the direct equivalent of R's scale() function.
#'
#' @param D \code{SummarizedExperiment} input.
#' @param center_data Mean-center data? Default: T.
#' @param scale_data Scale data to sd 1? Default: T.
#' @param ref_samples Expression for selecting samples to use for center and scale calculation.
#'
#' @return assay: Scaled data.
#'
#' @examples
#' \dontrun{# in the context of a SE pipeline
#' ... %>% mt_pre_trans_scale() %>% ...    # standard call, center and scale
#' ... %>% mt_pre_trans_scale(scale_data=F) %>% ...    # only mean centering
#' }
#'
#' @author JK
#'
#' @export
mt_pre_trans_scale <- function(D, center_data=T, scale_data=T, ref_samples) {

  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  stopifnot(is.logical(center_data))
  stopifnot(is.logical(scale_data))

  # scale
  if(!missing(ref_samples)){

    # get filtered samples
    filter_q <- dplyr::enquo(ref_samples)
    num_samp <- ncol(D)
    Da <- D %>% mti_format_se_samplewise()
    samples.used <- mti_filter_samples(Da, filter_q, num_samp)

    # get and filter
    Da_filtered <- t(assay(D))[samples.used,]

    if(center_data == T){
      # center by mean of selected samples
      da_f_means <- apply(Da_filtered, 2, mean, na.rm=T)
      Da <- sweep(Da, 2, da_f_means, FUN = "-")
    }
    if(scale_data == T){
      # scale by sd of selected samples
      da_f_sd <- apply(Da_filtered, 2, stats::sd, na.rm=T)
      Da <- sweep(Da, 2, da_f_sd, FUN = "/")
    }

    assay(D) = t(Da)

  } else {
    # if no sample filter given, just use scale() function
    assay(D) = t(scale(t(assay(D)),center=center_data,scale=scale_data))
  }

  # add status information
  funargs <- mti_funargs()
  metadata(D)$results %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf('scaled, center=%d, scale=%d', center_data, scale_data)
    )

  # return
  D

}
