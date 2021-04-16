#' Identifies Sample Outliers Using Univariate Approach
#'
#' A feature value is defined as an outlier if it is more than \code{thresh} standard deviations from the mean. A sample
#' is defined an outlier if more than \code{perc} of its features are univariate outliers.
#'
#' @param D \code{SummarizedExperiment} input.
#' @param reduce_dim If TRUE performs PCA-based dimensionality reduction before outlier detection. Can be used to apply multivariate
#'    outlier detection methods to low-rank datasets.
#' @param thresh Number of standard deviations to use as threshold to define an outlier. Default: 4.
#' @param perc Ratio of features that need to be outliers in order to consider the whole sample an outlier. Default: 0.5.
#'
#' @return colData: New columns including a binary vector and a numeric score vector.
#' @return $results$output: Returns the specific parameters used to determine outliers.
#'
#' @examples
#' \dontrun{# identify samples that have more than 50% univariate outliers
#' ... %>%
#'   mt_pre_outlier_detection_univariate(thresh=4, perc=0.5) %>%
#' ...}
#'
#' @author EB, JK
#'
#' @export
mt_pre_outlier_detection_univariate <- function(D, reduce_dim=F, thresh = 4, perc = 0.5) {

  # check arguments, SummarizedExperiment, and exactly one cutoff argument must be non-NA
  stopifnot("SummarizedExperiment" %in% class(D))

  X <- t(assay(D))
  X <- scale(X)

  # perform dimension reduction?
  if (reduce_dim) {
    ## Calculate the number of "independent features"
    ## As in Li and Ji, Heredity, 2005
    cordat <- stats::cor(X)
    eigenvals <- eigen(cordat)$values
    Meff <- sum( as.numeric(eigenvals >= 1) + (eigenvals - floor(eigenvals)) )
    ## reduce
    pca <- stats::prcomp(X)
    X <- pca[]
  }

  # compute univariate outliers
  H <- matrix(0, dim(X)[1], dim(X)[2])
  H[X>=thresh] <- 1
  # compute percentage of univariate outliers per sample
  score <- rowSums(H)/dim(H)[2]
  # define outliers
  out <- rep(FALSE,length(score))
  out[score > perc] <- TRUE

  # adding to colData
  colData(D)$outlier <- out
  colData(D)$score <- score
  colnames(colData(D))[colnames(colData(D))=="outlier"] <- paste0("outlier_", "univariate", sep="")
  colnames(colData(D))[colnames(colData(D))=="score"] <- paste0("score_", "univariate", sep="")

  l <- list(thresh=thresh, perc=perc)
  names(l) <- paste0(names(l), "_", "univariate", sep="")

  # add status information
  funargs <- mti_funargs()
  D %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("flagged %d %s outliers", sum(out, na.rm = TRUE), "univariate"),
      output = l
    )

  # return
  D

}
