#'Impute using minimum value
#'
#'Impute using minimum value of that feature across the samples.
#'
#' @param D \code{SummarizedExperiment} input.
#' @param verbose Whether to output intermediate steps. Default: F.
#'
#' @return assay: imputed data.
#'
#' @examples
#' \dontrun{# in the context of a SE pipeline
#' ... %>% mt_pre_impute_min() %>%
#' }
#'
#' @author RB
#'
#' @export
mt_pre_impute_min <- function(D, verbose=F) {

  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))

  my_min <- function(x) {
    if(all(is.na(x))){
      return(NA)
    } else{
      return(min(x, na.rm=T))
    }
  }

  df <- D %>% assay() %>%  data.frame()
  incom.obs <- which(apply(df,2,function(x) any(is.na(x))))
  if(verbose)message(paste0("Number of imcomplete observations: ", length(incom.obs)))
  all_nas <- length(which(apply(df,1,function(x) all(is.na(x)))))
  # impute
  df = apply(df, 1, function(x) {x[is.na(x)] <-  my_min(x); x} ) %>% t()
  assay(D) = df

  # add status information
  funargs <- mti_funargs()
  metadata(D)$results %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf('imputed via minimum value, %d features with all NAs, returned as NAs', all_nas)
    )

  # return
  D

}
