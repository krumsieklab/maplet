#' Load WCM core-format data
#'
#' Loads data from a WCM core format file (as they send it).
#'
#' @param D \code{SummarizedExperiment} input. Missing if first step in pipeline.
#' @param file Name of input Excel file.
#' @param sheet Name or number of sheet.
#' @param zero_to_na Replace zeros by NAs? Default: T.
#'
#' @return If first step in pipeline, creates \code{SummarizedExperiment} object. Populates empty assay, colData, and rowData data frames.
#'
#' @examples
#' \dontrun{D <-
#'   # load data in WCM format
#'   mt_load_wcm(file='DM369_Metabolites separated by nucleotides 12_6_18_altered.xlsx',
#'               sheet='peakHeight_metabolite_intensiti') %>%
#'   ...}
#'
#' @author JK
#'
#' @export
mt_load_wcm <- function(D, file, sheet, zero_to_na=T){

  # initialize outer result list
  result <- list()

  # validate arguments
  if (missing(file)) stop("file must be provided")
  if (missing(sheet)) stop("sheet must be provided")

  # save input information
  result$info$file <- file
  result$info$sheet <- paste(sheet, collapse = ", ")

  # get metadata from D if present
  if(!missing(D)){
    # validate SE
    if ("SummarizedExperiment" %in% class(D) == F) stop("D is not of class SummarizedExperiment")
    if (length(assays(D))!=0) stop("Passed SummarizedExperiment assay must be empty!")

    # get metadata
    result$meta <- metadata(D)
  }

  # load file in raw format
  raw = as.data.frame(readxl::read_excel(path=file, sheet=sheet)) %>% tibble::column_to_rownames("compound")
  names = rownames(raw)
  rownames(raw) = make.names(rownames(raw))
  # split off HMDB column if it exists
  if ("HMDB" %in% colnames(raw)) {
    metinfo = data.frame(name=names,HMDB=raw$HMDB)
    raw = raw %>% dplyr::select(-HMDB)
  } else {
    metinfo = data.frame(name=rownames(raw))
  }

  # zeros to NAs?
  if (zero_to_na) raw[raw==0] <- NA

  # generate summarized experiment
  D <- SummarizedExperiment(assay    = as.matrix(raw),
                            rowData  = metinfo,
                            colData  = data.frame(ID=colnames(raw)),
                            metadata = list(sessionInfo=utils::sessionInfo(), parseInfo=result$info))

  # add original metadata if exists
  if (!is.null(result$meta$results)) metadata(D)$results <- result$meta$results
  if (!is.null(result$meta$settings)) metadata(D)$settings <- result$meta$settings

  # add status information
  funargs <- mti_funargs()
  metadata(D)$results %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("loaded WCM file: %s, sheet: %s", basename(file), sheet)
    )

  # return
  D

}



