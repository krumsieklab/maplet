#' Identifies sample outliers using leverage approach
#'
#' Multivariate approach that uses leverage to define outliers. A sample is defined an outlier if its leverage is larger than
#' \code{thresh} times m/n (where m is the number of features and n is the number of samples in the dataset).
#'
#' @param D \code{SummarizedExperiment} input.
#' @param reduce_dim If TRUE performs PCA-based dimensionality reduction before outlier detection. Can be used to apply multivariate
#'    outlier detection methods to low-rank datasets.
#' @param thresh Number of m/n units to use as threshold to define an outlier. Default: 4.
#'
#' @return colData: New columns including a binary vector and a numeric score vector.
#' @return $results$output: Returns the specific parameters used to determine outliers.
#'
#' @examples
#' \dontrun{# identify multivariate outliers with a leverage >4m/n
#' ... %>%
#'   mt_pre_outlier_detection_leverage(thresh=4) %>%
#' ...}
#'
#' @author EB, JK
#'
#' @export
mt_pre_outlier_detection_leverage <- function(D, reduce_dim=F, thresh=4) {

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
    stop("The data matrix is not full-rank, leverage cannot be computed")

  # compute hat matrix with the Pivoted Cholesky factorization
  L <- t(suppressWarnings(chol(crossprod(X), pivot = TRUE)))
  r <- attr(L, "rank")
  piv <- attr(L, "pivot")
  Qt <- forwardsolve(L, t(X[, piv]), r)
  H <- crossprod(Qt)

  # extract leverage values
  score <- diag(H)
  # define outliers
  out <- rep(FALSE,length(score))
  out[score > thresh*sum(score)/dim(X)[1]] <- TRUE

  # adding to colData
  colData(D)$outlier <- out
  colData(D)$score <- score
  colnames(colData(D))[colnames(colData(D))=="outlier"] <- paste0("outlier_", "leverage", sep="")
  colnames(colData(D))[colnames(colData(D))=="score"] <- paste0("score_", "leverage", sep="")

  l <- list(threshold=thresh)
  names(l) <- paste0(names(l), "_", "leverage", sep="")

  # add status information
  funargs <- mti_funargs()
  D %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("flagged %d %s outliers", sum(out, na.rm = TRUE), "leverage"),
      output = l
    )

  # return
  D

}

