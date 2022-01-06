#' Output assay, colData and rowData into an Excel file
#'
#' Exports the current \code{SummarizedExperiment} (not the metadata) to an Excel file.
#'
#' @param D \code{SummarizedExperiment} input.
#' @param file Output file name to write to.
#' @param keep_list_col Whether to keep rowData columns containing list values. If TRUE, convert each list to a
#'    string with values separated by "|". If FALSE, the column is removed. Default: FALSE.
#'
#' @return Does not change the \code{SummarizedExperiment} object.
#'
#' @examples
#' \dontrun{%>% mt_write_se_xls(file = "out.xlsx") %>%}
#'
#' @author JK, KC
#'
#' @export
mt_write_se_xls <- function(D, file, keep_list_col = FALSE) {

  # verify that input is a SummarizedExperiment
  stopifnot("SummarizedExperiment" %in% class(D))
  stopifnot(is.character(file))


  # handle rowData list column(s)
  rd = rowData(D) %>% as.data.frame()
  list_cols <- lapply(seq_along(rd), function(x){if(is.list(rd[,x])){names(rd[x])}}) %>% unlist()
  if(length(list_cols) > 0){
    if(keep_list_col){
      # collapse list values to strings
      rd %<>% mutate_at(list_cols, collapse_list)
    }else{
      # remove list column(s)
      rd %<>% dplyr::select(-all_of(list_cols))
    }
  }

  # write out
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb,"assay")
  openxlsx::writeData(wb, "assay", assay(D), rowNames = T, colNames=T)
  openxlsx::addWorksheet(wb,"rowData")
  openxlsx::writeData(wb, "rowData", rd, rowNames = T)
  openxlsx::addWorksheet(wb,"colData")
  openxlsx::writeData(wb, "colData", colData(D) %>% as.data.frame())
  openxlsx::saveWorkbook (wb, file=file, overwrite=TRUE)

  # add status information
  funargs <- mti_funargs()
  D %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("Data exported to Excel file '%s'", file)
    )


  # pass SummarizedExperiment back, so pipeline can keep running
  D

}


#' Collapse List Value into String
#'
#' Accepts a list value and returns a string with each element of the list separated by '|'.
#'
#' @param x A list.
#'
#' @returns A string of concatenated list values.
#'
#' @noRd
collapse_list <- function(x){
  x <- lapply(x, function(y){if(!is.null(y)){
    paste0(y%>%unlist(),collapse = "|")
  }else{y}})
  x
}
