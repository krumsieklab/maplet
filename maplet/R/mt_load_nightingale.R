#' Load Nightingale-format data
#'
#' Loads data from a Nightingale format Excel file.
#'
#' @param D \code{SummarizedExperiment} input. Missing if first step in pipeline.
#' @param file Name of input excel file.
#' @param format_type Type of nightingale format: "single_sheet" or "multiple_sheets". Default: "multiple_sheets".
#' @param data_sheet If format_type is multiple sheets, name of sheet with data.
#' @param met_sheet If format type is multiple sheets, name of sheet with biomarker information. Default: "Biomarker annotations".
#' @param met_qc_sheet If format type is multiple sheets, name of sheet with biomarker qc tags. Default: "Tags per biomarker".
#' @param samp_qc_sheet If format type is multiple sheets, name of sheet with sample qc tags. Default: "Quality control tags and notes".
#' @param id_col OPTIONAL. Name of the column with sample IDs.
#'
#' @return Produces an initial \code{SummarizedExperiment}, with assay, colData, rowData, and metadata with first entry.
#'
#' @examples
#' \dontrun{D <-
#'   # load data
#'   mt_load_nightingale(file = =system.file("extdata", "example_data/sampledata.xlsx", package = "maplet")) %>%
#'   ...}
#'
#' @author RB
#'
#' @export
mt_load_nightingale <- function(D,
                                file,
                                data_sheet,
                                format_type = 'multiple_sheets',
                                met_sheet = 'Biomarker annotations',
                                met_qc_sheet = 'Tags per biomarker',
                                samp_qc_sheet = 'Quality control tags and notes',
                                id_col) {

  # initialize outer result list
  result <- list()

  # get metadata from D if present
  if(!missing(D)){
    # validate SE
    if ("SummarizedExperiment" %in% class(D) == F) stop("D is not of class SummarizedExperiment")
    if (length(assays(D))!=0) stop("Passed SummarizedExperiment assay must be empty!")

    # get metadata
    result$meta <- metadata(D)
  }

  if(format_type=='multiple_sheets'){
    if(missing(data_sheet)){data_sheet <- 'Results'}
    if(missing(id_col)){id_col <- 'Sample id'}
    D <- load_multiple_sheet_format(file=file,
                                    data_sheet=data_sheet, met_sheet=met_sheet,
                                    met_qc_sheet=met_qc_sheet, samp_qc_sheet=samp_qc_sheet, id_col=id_col)
    # add original metadata if exists
    if (!is.null(result$meta$results)) metadata(D)$results <- result$meta$results
    if (!is.null(result$meta$settings)) metadata(D)$settings <- result$meta$settings
  } else if (format_type=='single_sheet') {
    if(missing(data_sheet)){stop('data_sheet should be provided for single_sheet format!')}
    if(missing(id_col)){id_col <- 'sampleid'}
    D <- load_single_sheet_format(file=file,
                                    data_sheet=data_sheet, id_col=id_col)
    # add original metadata if exists
    if (!is.null(result$meta$results)) metadata(D)$results <- result$meta$results
    if (!is.null(result$meta$settings)) metadata(D)$settings <- result$meta$settings

  } else {stop(sprintf('Unknown format type, %s!', format_type))}

  # add status information
  funargs <- maplet:::mti_funargs()
  D %<>% 
    maplet:::mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("loaded Nightingale file: %s, sheet: %s", basename(file), data_sheet)
    )

  # return
  D

}

get_data <- function(mat, met_info, id_col, samp_qc_sheet=F){
  # find last metabolite column
  imetlast <- max(which(apply(is.na(mat),2,sum)<dim(mat)[1]))
  # find last sample row
  isamplast <- max(which(apply(is.na(mat),1,sum)<dim(mat)[2]))
  # find row with table header
  tab_header <- min(which(!is.na(mat[, 1]))) + 1
  # find the first row of the data
  data_start <- min(which(!is.na(mat[(tab_header+1):isamplast, 1]))) + (tab_header + 1)
  # subset to data rows and columns
  if(samp_qc_sheet){
    col_ids <- c(which(mat[tab_header, ] %in%id_col): imetlast)
    col_names <- unlist(c(id_col, mat[(tab_header-1), col_ids[-1]]))
    mat <- mat[data_start:isamplast, col_ids]
    names(mat) <- col_names
  } else{
    # find the names of the table
    col_names <- c(id_col, met_info[['Excel_column_name']])
    mat <- mat[data_start:isamplast, which(mat[tab_header, ] %in%col_names)]
    names(mat) <- col_names
  }
  return(mat)
}

load_multiple_sheet_format <- function(file,
                                       data_sheet, met_sheet, met_qc_sheet, samp_qc_sheet, id_col){
  # using readxl package:
  raw <- readxl::read_excel(path=file, sheet=data_sheet, col_names = F)
  met_info <- readxl::read_excel(path=file, sheet=met_sheet, col_names = T)
  qc_met <- readxl::read_excel(path=file, sheet=met_qc_sheet, col_names = F)
  qc_sample <- readxl::read_excel(path=file, sheet=samp_qc_sheet, col_names = F)
  # convert any spaces in the colnames to underscores
  colnames(met_info) <- gsub(" ", "_", colnames(met_info))

  result <- list()

  raw_data <- get_data(mat=raw, met_info, id_col)
  qc_met <- get_data(mat=qc_met, met_info, id_col) %>% t() %>% data.frame()
  names(qc_met) <- unlist(qc_met[1, ])
  qc_met <- qc_met[-1, ] %>% mutate('Excel_column_name' = rownames(qc_met[-1, ]))
  qc_sample <- get_data(mat=qc_sample, met_info, id_col, samp_qc_sheet = T)
  # add sample information
  result$sampleinfo <- data.frame(raw_data %>% select(id_col), stringsAsFactors = F, check.names = F) %>%
    left_join(qc_sample, by=id_col)
  # add metabolite information
  result$metinfo <- data.frame(met_info, check.names = F) %>%
    left_join(qc_met, by='Excel_column_name')
  # add data
  raw_data <- raw_data %>% select(-id_col) %>%
    mutate_all(as.matrix) %>% mutate_all(as.numeric)
  result$data <- data.frame(raw_data)

  # add display name
  result$metinfo$name   <- result$metinfo$Excel_column_name
  # fix variable names
  colnames(result$data) <- result$metinfo$CSV_column_name
  # generate summarized experiment
  D <- SummarizedExperiment(assay    = t(result$data),
                            colData  = result$sampleinfo,
                            rowData  = result$metinfo,
                            metadata = list(sessionInfo=utils::sessionInfo()))

  # ensure colnames and rownames exist
  if (is.null(colnames(D))) colnames(D) <- 1:ncol(D)
  if (is.null(rownames(D))) rownames(D) <- result$metinfo$Excel_column_name
  # return SummarizedExperiment
  return(D)
}
load_single_sheet_format <- function (file=file,
                         data_sheet=data_sheet, id_col=id_col){
  # using readxl package:
  raw <- readxl::read_excel(path=file, sheet=data_sheet, col_names = F)
  # find last metabolite column
  imetlast <- max(which(apply(is.na(raw),2,sum)<dim(raw)[1]))
  # find last sample row
  isamplast <- max(which(apply(is.na(raw),1,sum)<dim(raw)[2]))
  # find row with table header
  tab_header <- min(which(!is.na(raw[, 1]))) +1
  # find the first row of the data
  data_start <- 1 + (min(which(!is.na(raw[(tab_header+1):isamplast, 1]))) + (tab_header + 1))
  # subset to data rows and columns
  col_ids <- c(which(raw[tab_header, ] %in%id_col): imetlast)
  col_names <- unlist(c(id_col, raw[tab_header, col_ids[-1]]))
  raw_data <- raw[data_start:isamplast, col_ids]
  names(raw_data) <- col_names
  result <- list()
  # add sample information
  result$sampleinfo <- raw[data_start:isamplast, 1:col_ids[1]]
  names(result$sampleinfo) <- c(unlist(raw[(tab_header-1), 1:(col_ids[1]-1)]), id_col)
  result$sampleinfo %<>% select(id_col, everything())
  # add metabolite information
  result$metinfo <- data.frame(name=unlist(raw[tab_header, col_ids[-1]]),
                               fullname=unlist(raw[(tab_header-1), col_ids[-1]]),
                               check.names = F)
  # add data
  raw_data <- raw_data %>% select(-id_col) %>%
    mutate_all(as.matrix) %>% mutate_all(as.numeric)
  result$data <- data.frame(raw_data)

  # fix variable names
  colnames(result$data) <- result$metinfo$name
  # generate summarized experiment
  D <- SummarizedExperiment(assay    = t(result$data),
                            colData  = result$sampleinfo,
                            rowData  = result$metinfo,
                            metadata = list(sessionInfo=utils::sessionInfo()))

  # ensure colnames and rownames exist
  if (is.null(colnames(D))) colnames(D) <- 1:ncol(D)
  if (is.null(rownames(D))) rownames(D) <- result$metinfo$name
  # return SummarizedExperiment
  return(D)
}
