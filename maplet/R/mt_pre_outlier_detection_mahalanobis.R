#' Identifies sample outliers using Mahalanobis approach
#'
#' Multivariate approach that uses the Mahalanobis distance to define outliers. A sample is defined an outlier if its Mahalanobis 
#' distance is in the \code{pval} quantile of the chi-square distribution.
#'
#' @param D \code{SummarizedExperiment} input.
#' @param reduce_dim boolean, if TRUE performs PCA-based dimensionality reduction before outlier detection. Can be used to apply 
#'    multivariate outlier detection methods to low-rank datasets.
#' @param pval Value between 0 and 1. Quantile of the chi-squared distribution to use as threshold to define a sample an outlier. 
#'    Default: 0.01.
#'
#' @return colData: New columns including a binary vector and a numeric score vector.
#' @return $results$output: Returns the specific parameters used to determine outliers.
#'
#' @examples
#' \dontrun{# identify multivariate outliers
#' ... %>%
#'   mt_pre_outlier_detection_mahalanobis(pval=0.01) %>%
#' ...}
#'
#' @author EB, JK
#'
#' @export
mt_pre_outlier_detection_mahalanobis <- function(D, reduce_dim=F, pval=0.01) {
  
  # check arguments, SummarizedExperiment, and exactly one cutoff argument must be non-NA
  stopifnot("SummarizedExperiment" %in% class(D))
  
  X <- t(assay(D))
  X <- scale(X)
  
  if(any(is.na(X)))
    stop("Missing values found in the data matrix. Multivariate outlier detection approaches cannot handle missing values.")
  
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
  
  if(!limma::is.fullrank(X))
    stop("The data matrix is not full-rank, Mahalanobis cannot be computed")

  
  # Calculate covariance matrix
  cov_mat <- stats::cov(X)
  # Get mahalanobis distance
  score <- stats::mahalanobis(X, colMeans(X), cov_mat)
  # Define outliers based on threshold
  out <- rep(FALSE,length(score))
  out[score > stats::qchisq(1 - pval/nrow(X), df=ncol(X))] <- TRUE
  # adding to colData
  colData(D)$outlier <- out
  colData(D)$score <- score
  colnames(colData(D))[colnames(colData(D))=="outlier"] <- paste0("outlier_", "mahalanobis", sep="")
  colnames(colData(D))[colnames(colData(D))=="score"] <- paste0("score_", "mahalanobis", sep="")
  
  l <- list(pval=pval)
  names(l) <- paste0(names(l), "_", "mahalanobis", sep="")
  
  # add status information
  funargs <- mti_funargs()
  D %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("flagged %d %s outliers", sum(out, na.rm = TRUE), "mahalanobis"),
      output = l
    )
  
  # return
  D
  
}

