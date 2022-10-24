#' Rename selected columns
#'
#' Replaces existing column names with new column names using dplyr's rename() function. Both old
#' and new column names should be provided as a vector. Works on either colData or rowData
#' annotation data frames.
#'
#' @param D \code{SummarizedExperiment} input.
#' @param anno_type Either "samples" (colData) or "features" (rowData).
#' @param old_col_names Vector of names of existing annotation columns to be changed.
#' @param new_col_names Vector of new names to replace existing names.
#'
#' @return colData or rowData: Column names replaced.
#'
#' @examples
#' \dontrun{# Replace column names in colData data frame
#' ...  %>%
#'  mt_anno_rename(anno_type="samples", old_col_names=c(), new_col_names = c()) %>% ...}
#'
#' @author KC
#'
#' @export
mt_anno_rename <- function(D, anno_type, old_col_names, new_col_names){

  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  if (!(anno_type %in% c("samples","features"))) stop("anno_type must be either 'samples' or 'features'")
  if(length(old_col_names) != length(new_col_names)) stop("Arguments old_col_names and new_col_names must have the same length!")

  if (anno_type=="samples") {
    cn <- colnames(D) # ensure colnames not destroyed
    cd <- colData(D) %>% as.data.frame()
    cd %<>% dplyr::rename_at(dplyr::vars(old_col_names), function(x) new_col_names)
    colData(D) <- DataFrame(cd)
    colnames(D) <- cn

  } else if (anno_type=="features") {
    cn <- colnames(D) # ensure colnames not destroyed
    rd <- rowData(D) %>% as.data.frame()
    rd %<>% dplyr::rename_at(dplyr::vars(old_col_names), function(x) new_col_names)
    rowData(D) <- DataFrame(rd)
    colnames(D) <- cn

  } else {
    stop('bug')
  }


  ## add status information & plot
  funargs <- mti_funargs()
  D %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("Renamed %s annotation columns from: %s, to: %s", anno_type, old_col_names, new_col_names)
    )

  ## return
  D

}
