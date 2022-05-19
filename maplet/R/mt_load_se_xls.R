#' Load data matrix from Excel file - mt_write_se_xls format
#'
#' Load data and annotation matrices from Excel sheet created by the mt_write_se_xls function. This
#' function is intended to be used on data produced by the mt_write_se_xls without manual
#' alteration.
#'
#' @param D \code{SummarizedExperiment} input. Missing if first step in pipeline.
#' @param file Name of input Excel file.
#' @param sheet_names Vector of the sheet names to read. Three sheet names must be provided and must
#'    be ordered as such: (1) assay, (2) rowData, and (3) colData.
#'    Default: c("assay", "rowData", "colData").
#'
#' @return If first step in pipeline, creates \code{SummarizedExperiment} object. Populates empty
#'         assay, colData, and rowData data frames.
#'
#' @examples
#' # Load data, rowData, and colData from Excel file
#' \dontrun{D <-
#'   # load raw data
#'   mt_load_se_xls(file=file, sheet_names=c("data", "metabolites", "samples")) %>%
#'   }
#'
#' @author JK, KC
#'
#' @export
mt_load_se_xls <- function(D,
                           file,
                           sheet_names=c("assay", "rowData", "colData")){

  result=list()
  sheet_text = paste0(sheet_names, collapse = ", ")

  ### validate arguments ------
  if(missing(file)) stop("file must be provided")
  if(length(sheet_names) != 3) stop("\'sheet_names\' must be given exactly three values.")

  ### check for SE ------
  # get metadata from D if present
  if(!missing(D)){
    if ("SummarizedExperiment" %in% class(D) == F) stop("D is not of class SummarizedExperiment")
    if (length(assays(D))!=0) stop("Passed SummarizedExperiment assay must be empty!")
    result$meta <- metadata(D)
  }

  ### get assay, rowData, and colData matrices ------
  df <- as.data.frame(readxl::read_excel(path=file,sheet=sheet_names[1], col_names=T))
  df %<>% tibble::column_to_rownames(colnames(df)[1])
  assay <- as.matrix(df)
  rd <- as.data.frame(readxl::read_excel(path=file,sheet=sheet_names[2],col_names=T))
  rd %<>% tibble::column_to_rownames(colnames(rd)[1])
  cd <- as.data.frame(readxl::read_excel(path=file,sheet=sheet_names[3],col_names=T))

  ### construct SummarizedExperiment ------
  D <- SummarizedExperiment(assay = assay,
                            rowData = rd,
                            colData = cd,
                            metadata = list(sessionInfo=utils::sessionInfo()))
  # add original metadata if exists
  if (!is.null(result$meta$results)) metadata(D)$results <- result$meta$results
  if (!is.null(result$meta$settings)) metadata(D)$settings <- result$meta$settings

  ### add status information ------
  funargs <- mti_funargs()
  D %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = glue::glue("Loaded assay from maplet-formatted Excel file: {basename(file)}, sheets: {sheet_text}")
      )

  D

}
