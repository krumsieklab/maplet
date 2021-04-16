#' Sets values of outliers in samples to NA
#'
#' Uses univariate approach to identify outliers in the samples and sets the value of these outliers to NA.
#'
#' @param D \code{SummarizedExperiment} input.
#' @param thresh Number of standard deviations or m/n units to use as threshold to define the outlier, if not provided, default is
#'    correction by sample size. This parameter is only used if sample_num_correction is False.
#' @param sample_num_correction Whether number of outliers should depend on number of samples or hard sd cutoff. If true, thresh is ignored
#' @param alpha The percentage of points we will consider outliers based on the assumption of a normal distribution (alpha/2 on either tail). Only used if sample_num_correction is True
#' @param use_quant Optional, if to flag the outliers in the quantile provided in the quantile parameter
#' @param quant_thresh If use_quant is true provide a quantile
#'
#' @return assay: NA values where outliers used to be
#'
#' @examples
#' \dontrun{... %>%
#'   mt_pre_outlier_to_na(thresh=3, sample_num_correction=F) %>%
#' ...}
#'
#' @author AS, RB
#'
#' @export
mt_pre_outlier_to_na <- function(
  D,
  thresh=NA,
  sample_num_correction = T,
  alpha = 0.05,
  use_quant=FALSE, # Whether to correct by quantile method
  quant_thresh=0.025, # Quantile
  ...
) {

  # check arguments, SummarizedExperiment, and exactly one cutoff argument must be non-NA
  stopifnot("SummarizedExperiment" %in% class(D))

  X <- t(assay(D))
  X <- scale(X)

  if(use_quant==FALSE){
    if(is.na(thresh) & sample_num_correction == F){
      stop("Threshold must be provided if not corrected by sample numbers")
    }
    if(is.na(thresh)){
      numsamp=nrow(X)
      tail=alpha/numsamp
      thresh = stats::qnorm( 1 - (tail/2) )
    }
  } else if (use_quant==TRUE){
    if(is.na(quant_thresh)){
      stop("Quantile must be provided")
    }
    # compute threshold based on the given quantile and sample size
    thresh <- abs(stats::qnorm((quant_thresh/2)/nrow(X)))
  }
  # compute univariate outliers
  H <- matrix(F, dim(X)[1], dim(X)[2])
  H[abs(X)>=thresh] <- T
  # change outliers to NA
  assay(D)[t(H)] <- NA
  # add status information
  funargs <- mti_funargs()
  D %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("flagged %d outliers", sum(H))
    )

  # return
  D

}
