#' Median batch correction
#'
#' Same approach that Metabolon uses for runday correction. Works for both logged and non-logged data. Subtracts (logged) or
#' divides (non-logged) feature values by the median value per batch.
#'
#' @param D \code{SummarizedExperiment} input.
#' @param batch_col Sample annotation (colData) column name that contains batch assignment.
#' @param ref_samples Expression to filter out reference samples to use from colData.
#'
#' @return assay: Batch-corrected version.
#'
#' @examples
#' \dontrun{... %>% mt_pre_batch_median(batch_col="BATCH") %>% ...
#' }
#'
#' @author JK
#'
#' @export
mt_pre_batch_median = function(D, batch_col, ref_samples) {

  # validate and extract arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  X = t(assay(D))

  # get variable
  if (!(batch_col %in% colnames(colData(D)))) stop(sprintf("'%s' not found in sample annotations.", batch_col))
  b = colData(D)[[batch_col]]
  # ensure it's a factor or character vector
  if (!(is.character(b) || is.factor(b))) stop(sprintf("'%s' has to be character, factor, or numeric.", batch_col))
  b = as.factor(b)

  # no negative values allowed
  if (min(X,na.rm=T)<0) stop("Matrix contains negative values.")
  # select correct operator for logged or non-logged data
  is_logged <- maplet:::mti_check_is_logged(D)
  if(is_logged){
    op <- `-`
    opstr <- "-"
  } else {
    op <- `/`
    opstr <- "/"
  }


  # get samples to use for median calculation
  if (!missing(ref_samples)) {
    # select those samples according to the expression
    ref_samples_q <- dplyr::enquo(ref_samples)
    inds <- colData(D) %>%
      as.data.frame() %>%
      dplyr::mutate(tmporder=1:ncol(D)) %>%
      dplyr::filter(!!ref_samples_q) %>%
      .$tmporder
    use_samples <- rep(F, ncol(D))
    use_samples[inds] <- T
  } else {
    # all samples
    use_samples <- rep(T, ncol(D))
  }

  # median per feature
  for (i in 1:length(levels(b))) {
    batch = levels(b)[i]
    # check that there are any reference samples available
    if (sum(b==batch & use_samples) == 0) stop(sprintf("No reference samples for batch '%s'", batch))
    # build median vector for this batch
    med <- X[b==batch & use_samples,,drop=F] %>% apply(2, stats::median, na.rm=T)
    # transform into matrix of the same size as the batch
    med_matrix <- replicate(sum(b==batch), med) %>% t()
    # median normalize
    X[b==batch,] <-  op(X[b==batch,], med_matrix)
  }

  # replace original assay with batch corrected assay
  assay(D) = t(X)

  # ref samples logging string
  refadd <- if(missing(ref_samples)){""}else{sprintf(": %s", as.character(dplyr::enquo(ref_samples)))}

  # add status information
  funargs <- maplet:::mti_funargs()
  metadata(D)$results %<>%
    maplet:::mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("median-scaling per batch in '%s', based on %d reference samples%s",batch_col,sum(use_samples),refadd)
    )

  # return
  D


}
