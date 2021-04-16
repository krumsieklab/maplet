#' Create new variable by expression
#'
#' Creates a new variable using dplyr's mutate() mechanism. Formula can use any column in the respective colData/rowData. Can create
#' a new annotation column or replace an existing one (if col_name and term arguments refer to same column).
#'
#' @param D \code{SummarizedExperiment} input.
#' @param anno_type Either "samples" (colData) or "features" (rowData).
#' @param col_name Name of new (or existing) annotation column to store result of expression.
#' @param term Mutate term to forward to \code{dplyr::mutate()}. Can include any columns from colData or rowData.
#'
#' @return colData or rowData: New (or existing) annotation column added (or replaced).
#'
#' @examples
#' \dontrun{# Convert numeric sample annotation field 'num' to factor
#' ...  %>%
#'  mt_anno_mutate(anno_type="samples", col_name="num", term = as.factor(num)) %>% ...}
#'
#' @author JK
#'
#' @importFrom data.table :=
#'
#' @export
mt_anno_mutate <- function(D, anno_type, col_name, term) {

  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  if (!(anno_type %in% c("samples","features"))) stop("anno_type must be either 'samples' or 'features'")

  # run mutate argument in correct direction
  x <- dplyr::enquo(term)

  if (anno_type=="samples") {
    cn <- colnames(D) # mutate destroys colnames
    cd <- colData(D) %>% as.data.frame()
    cd %<>% dplyr::mutate(!!col_name := !!x)
    colData(D) <- DataFrame(cd)
    colnames(D) <- cn
  } else if (anno_type=="features") {
    cn <- colnames(D) # mutate destroys colnames
    rd <- rowData(D) %>% as.data.frame()
    rd %<>% dplyr::mutate(!!col_name := !!x)
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
      logtxt = sprintf("added variable to %s: %s := %s", anno_type, col_name, as.character(x))
    )

  ## return
  D

}
