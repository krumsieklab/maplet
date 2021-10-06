#' Remove outliers using local outlier factor
#'
#' Identifies outliers using the Local Outlier Factor method implmented in the \code{bigutilsr} package. 
#' See the following publication for details: 
#' \href{https://www.dbs.ifi.lmu.de/Publikationen/Papers/LOF.pdf}{https://www.dbs.ifi.lmu.de/Publikationen/Papers/LOF.pdf}.
#' For additional details about this method, including how to choose the threshold, refer to the following:
#' \href{https://www.r-bloggers.com/2019/08/detecting-outlier-samples-in-pca/}{https://www.r-bloggers.com/2019/08/detecting-outlier-samples-in-pca/}.
#'
#' @param D \code{SummarizedExperiment} input.
#' @param seq_k Sequence of numbers of nearest neighbors to use. Default: c(5, seq(10, nrow(df), by=10)).
#' 
#' @return assay: Outlier samples removed.
#' 
#' @examples
#' \dontrun{
#'   D <- D %>% mt_pre_outlier_lof() %>%
#' ...} 
#' 
#' @export
mt_pre_outlier_lof <- function(D,
                               seq_k = c(5, seq(10, nrow(df)-1, by=10))){
  
  # Validate arguments
  if("SummarizedExperiment" %in% class(D) == F) stop("D must be a SummarizedExperiment object.")
  
  # transpose data so samples in rows, features in columns
  df <- t(assay(D))
  
  # identify local outlier factor using sequence of different neighbours
  lof <- bigutilsr::LOF(df, seq_k = seq_k)
  
  # detect outliers based on departure from histogram
  out_idx <- which(lof>hist_out(lof)$lim[2])
  # remove any detected outliers
  if(length(out_idx)>0) D <- D[, -out_idx]
  
  # add status information
  funargs <- mti_funargs()
  D %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = glue::glue("Removed {length(out_idx)} outliers out of {nrow(df)} samples using LOF method.")
    )
  
  # return
  D
  
}