#' Output assay, colData and rowData into an Excel file
#'
#' Exports the current \code{SummarizedExperiment} (not the metadata) to an Excel file.
#'
#' @param D \code{SummarizedExperiment} input.
#' @param file Output filename to write to.
#'
#' @return Does not change the \code{SummarizedExperiment} object.
#'
#' @examples
#' \dontrun{%>% mt_write_se_xls(file = "out.xlsx") %>%}
#'
#' @author JK
#'
#' @export
mt_write_se_xls <- function(D, file) {

  # verify that input is a SummarizedExperiment
  stopifnot("SummarizedExperiment" %in% class(D))
  stopifnot(is.character(file))

  # rowdata lists must be collapsed
  rd = rowData(D) %>% as.data.frame()
  rd[] <- lapply(rd, function(x){if(is.list(x)){paste(x%>%unlist(),collapse="|")}else{x}})

  # write out
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb,"assay")
  openxlsx::writeData(wb, "assay", assay(D), rowNames = T, colNames=T)
  openxlsx::addWorksheet(wb,"rowData")
  openxlsx::writeData(wb, "rowData", rowData(D) %>% as.data.frame(), rowNames = T)
  openxlsx::addWorksheet(wb,"colData")
  openxlsx::writeData(wb, "colData", colData(D) %>% as.data.frame())
  openxlsx::saveWorkbook (wb, file=file, overwrite=TRUE)

  # add status information
  funargs <- mti_funargs()
  metadata(D)$results %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("Data exported to Excel file '%s'", file)
    )


  # pass SummarizedExperiment back, so pipeline can keep running
  D

}
