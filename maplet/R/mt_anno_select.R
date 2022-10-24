#' Select columns from annotation data frame
#'
#' Selects a subset of columns from either the colData or rowData annotation data frame. The user
#' can either choose to provide a list of columns to keep or to provide a list of columns to
#' remove. Uses the select() function from the dplyr package.
#'
#' @param D \code{SummarizedExperiment} input.
#' @param anno_type Either "samples" (colData) or "features" (rowData).
#' @param term Selection expression to forward to \code{dplyr::select()}. Can include any columns
#'    from colData or rowData.
#'
#' @return colData or rowData: Column(s) selected or removed.
#'
#' @examples
#' \dontrun{# Keep columns from colData data frame
#' ...  %>%
#'  mt_anno_select(anno_type="samples", term= ) %>%
#'  # remove columns from rowData data frame
#'  mt_anno_select(anno_type="features", term = )...}
#'
#' @author JK, KC
#'
#' @export
mt_anno_select <- function(D, anno_type, term){

  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  if (!(anno_type %in% c("samples","features"))) stop("anno_type must be either 'samples' or 'features'")

  quo_term = enquo(term)

  if (anno_type=="samples") {
    cn <- colnames(D) # ensure colnames not destroyed
    cd <- colData(D) %>% as.data.frame()
    cd %<>% dplyr::select(!!quo_term)
    colData(D) <- DataFrame(cd)
    colnames(D) <- cn

  } else if (anno_type=="features") {
    cn <- colnames(D) # ensure colnames not destroyed
    rd <- rowData(D) %>% as.data.frame()
    rd %<>% dplyr::select(!!quo_term)
    rowData(D) <- DataFrame(rd)
    colnames(D) <- cn

  } else {
    stop('bug')
  }

  # KC: To-Do - replace as.character() in logtxt
  ## add status information & plot
  funargs <- mti_funargs()
  D %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("Subset %s annotation columns using expression: %s", anno_type, as.character(quo_term))
    )

  ## return
  D

}
