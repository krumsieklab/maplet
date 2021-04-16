#' Apply function to feature or sample annotation column
#'
#' Works in the same way as lapply/sapply, with a function to be executed on each entry of a given feature or sample annotation.
#'
#' @param D \code{SummarizedExperiment} input.
#' @param anno_type Either "samples" (colData) or "features" (rowData).
#' @param col_name Name of column to apply function to.
#' @param fun Function to be applied.
#'
#' @return colData or rowData: Changes the contents of a column.
#'
#' @examples
#'  \dontrun{... %>%
#'  # ensure factor for casecontrol variable
#'  mt_anno_apply(anno_type='samples',
#'                col_name='casecontrol',
#'                fun=as.factor) %>%
#'  ...}
#'
#' @author JK
#'
#' @export
mt_anno_apply <- function(D, anno_type, col_name, fun) {

  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  if(!(any(c("samples","features")%in%anno_type)))stop("'anno_type' must be 'samples' or 'features'")

  # get data frame
  df = if(anno_type=="samples"){colData(D)}else{rowData(D)}

  # get variable
  if (!(col_name %in% colnames(df))) stop(sprintf("'%s' not found in %s annotations.", col_name, ifelse(anno_type=="samples","sample","feature")))
  p = colData(D)[[col_name]]

  # apply function
  pnew <- sapply(p, fun)

  # write back
  if (anno_type=="samples") {
    colData(D)[[col_name]] <- pnew
  } else {
    rowData(D)[[col_name]] <- pnew
  }


  ## add status information & plot
  funargs <- mti_funargs()
  D %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("Transformed column '%s' of %s annotations", col_name, ifelse(anno_type=="samples","sample","feature"))
    )
  ## return
  D

}






