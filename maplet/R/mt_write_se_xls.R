#' Output assay, rowData and colData into an Excel file
#'
#' Exports the current \code{SummarizedExperiment} (not the metadata) to an Excel file.
#'
#' @param D \code{SummarizedExperiment} input.
#' @param file Output file name to write to.
#' @param keep_list_col Whether to keep rowData columns containing list values. If TRUE, convert each list to a
#'    string with values separated by "|". If FALSE, the column is removed. Default: FALSE.
#' @param sheet_names Vector of names for sheets, in the following order: (1) assay, (2) rowData, and
#'    (3) colData. Default: c("assay", "rowData", "colData").
#'
#' @return Does not change the \code{SummarizedExperiment} object.
#'
#' @examples
#' \dontrun{%>% mt_write_se_xls(file = "out.xlsx", sheet_names = c("data", "metabolites", "samples")) %>%}
#'
#' @author JK, KC
#'
#' @export
mt_write_se_xls <- function(D,
                            file,
                            keep_list_col = FALSE,
                            sheet_names = c("assay", "rowData", "colData")) {

  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  stopifnot(is.character(file))
  if(length(sheet_names) != 3) stop("\'sheet_names\' must be given exactly three values.")

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
  openxlsx::addWorksheet(wb,sheet_names[1])
  openxlsx::writeData(wb, sheet_names[1], assay(D), rowNames = T, colNames=T)
  openxlsx::addWorksheet(wb,sheet_names[2])
  openxlsx::writeData(wb, sheet_names[2], rd, rowNames = T)
  openxlsx::addWorksheet(wb,sheet_names[3])
  openxlsx::writeData(wb, sheet_names[3], colData(D) %>% as.data.frame())
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
