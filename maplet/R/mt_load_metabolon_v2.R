#' Load Metabolon-format data (new format)
#'
#' For Metabolon-format version used after 2020-10-29. Loads data from a NEW Metabolon-format Excel file. Needs to be in the original
#' "Client Data Table" new format that they deliver.
#'
#' @param D \code{SummarizedExperiment} input. Missing if first step in pipeline.
#' @param file Name of input Excel file.
#' @param data_sheet Name of sheet with raw data.
#' @param met_sheet Name of sheet with metabolite info.
#' @param samp_sheet Name of sheet with sample info.
#'
#' @return If first step in pipeline, creates \code{SummarizedExperiment} object. Populates empty assay, colData, and rowData data frames.
#'
#' @examples
#' \dontrun{D <-
#'   # load data
#'   mt_load_metabolon_v2(file=system.file("extdata", "example_data/sampledata.xlsx", package = "maplet"),
#'                        data_sheet="Raw Data",
#'                        met_sheet="Chemical Annotation",
#'                        samp_sheet="Sample Meta Data") %>%
#'   ...}
#'
#' @author RB
#'
#' @export
mt_load_metabolon_v2 <- function(D, file, data_sheet, met_sheet, samp_sheet) {

  # initialize outer result list
  result <- list()

  # validate arguments
  if (missing(file)) stop("file must be provided.")
  if (missing(data_sheet)) stop("sheet must be provided.")
  if (missing(met_sheet)) stop("met_sheet must be provided.")
  if (missing(samp_sheet)) stop("samp_sheet must be provided.")
  

  # get metadata from D if present
  if(!missing(D)){
    # validate SE
    if ("SummarizedExperiment" %in% class(D) == F) stop("D is not of class SummarizedExperiment")
    if (length(assays(D))!=0) stop("Passed SummarizedExperiment assay must be empty!")

    # get metadata
    result$meta <- metadata(D)
  }

  cols2discard <- c('PARENT_SAMPLE_NAME', 'SAMPLE_TYPE', 'CLIENT_IDENTIFIER',
                    'MATRIX_DESIGNATION')
  # using readxl package:
  raw = readxl::read_excel(path=file, sheet=data_sheet, col_names = T)
  met_info = readxl::read_excel(path=file, sheet=met_sheet, col_names = T)
  sample_info = readxl::read_excel(path=file, sheet=samp_sheet, col_names = T)

  result$data <- raw %>% select(-all_of(cols2discard))
  result$metinfo <- met_info[match(names(result$data), met_info$CHEM_ID), ]
  result$sampleinfo <- sample_info[match(raw$PARENT_SAMPLE_NAME, sample_info$PARENT_SAMPLE_NAME), ]

  # as column names, use "CHEMICAL_NAME", if available
  if ("CHEMICAL_NAME" %in% colnames(result$metinfo)) {
    colnames(result$data) = result$metinfo$CHEMICAL_NAME
  } else {
    colnames(result$data) = c()
  }

  # as row names, use "PARENT_SAMPLE_NAME", if available
  if ("PARENT_SAMPLE_NAME" %in% colnames(result$sampleinfo)) {
    rownames(result$data) = result$sampleinfo$PARENT_SAMPLE_NAME
  } else {
    rownames(result$data) = c()
  }

  # return SummarizedExperiment

  # add display name
  result$metinfo$name   <- result$metinfo$CHEMICAL_NAME
  # fix variable names
  colnames(result$data) <- result$metinfo$CHEMICAL_NAME %>% make.names()
  # generate summarized experiment
  D <- SummarizedExperiment(assay    = t(result$data),
                       colData  = result$sampleinfo,
                       rowData  = result$metinfo,
                       metadata = list(sessionInfo=utils::sessionInfo()))

  # add original metadata if exists
  if (!is.null(result$meta$results)) metadata(D)$results <- result$meta$results
  if (!is.null(result$meta$settings)) metadata(D)$settings <- result$meta$settings

  # ensure colnames and rownames exist
  if (is.null(colnames(D))) colnames(D) <- 1:ncol(D)
  if (is.null(rownames(D))) rownames(D) <- result$metinfo$CHEMICAL_NAME

  # add status information
  funargs <- maplet:::mti_funargs()
  D %<>% 
    maplet:::mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("loaded Metabolon file: %s, sheets: %s, %s, %s",
                       basename(file), data_sheet, met_sheet, samp_sheet)
    )

  # return
  D

}
