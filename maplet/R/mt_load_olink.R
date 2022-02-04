#' Load Olink-format data
#'
#' Loads data from an Olink-format Excel file.
#'
#'
#' @param D \code{SummarizedExperiment} input. Missing if first step in pipeline.
#' @param file Name of Excel file exported from NPX manger.
#'
#' @return If first step in pipeline, creates \code{SummarizedExperiment} object. Populates empty assay (note: 2^NPX), colData, and rowData data frames.
#'
#' @examples
#' \dontrun{D <-
#'   # load data
#'   mt_load_olink(file=system.file("extdata", "example_data/sampledata.xlsx", package = "maplet")) %>%
#'   ...}
#'
#' @author RB, KC
#'
#' @import readxl
#' @import tidyverse
#'
#' @export
mt_load_olink <- function(D, file){
  # columns to look for in Olink file
  REQ_COLS <- c("Panel", "Assay", "Uniprot ID", "OlinkID", "LOD") # required
  OPT_COLS <- c("Missing Data freq.", "Normalization") # optional
  # name of analyte data frame (rowData) column names corresponding to Olink row names
  ANA_COL_NAMES <- c("Panel"="Panel", "Uniprot ID"="UniprotID", "OlinkID"="OlinkID", "LOD"="LOD",
                     "Missing Data freq."="MissingDataFreq", "Assay"="name", "Normalization"="Normalization")

  # initialize result list
  result <- list()

  # get metadata from D if present
  if(!missing(D)){
    # validate SE
    if ("SummarizedExperiment" %in% class(D) == F) stop("D is not of class SummarizedExperiment")
    if (length(assays(D))!=0) stop("Passed SummarizedExperiment assay must be empty!")

    # get metadata
    result$meta <- metadata(D)
  }

  # read the Olink file
  mat <- read.xlsx(file)
  # find last column with some info
  imetlast <- max(which(apply(is.na(mat),2,sum)<dim(mat)[1]))
  # find last row with some info
  isamplast <- max(which(apply(is.na(mat),1,sum)<dim(mat)[2]))
  # find indices of Olink columns
  olink_cols <- c(REQ_COLS, OPT_COLS)
  col_index <- sapply(olink_cols, function(x){grep(x, mat[,1])}) %>% unlist()
  # throw an error if required columns missing
  missing_cols <- REQ_COLS[REQ_COLS %in% names(col_index)==F]
  if(length(missing_cols)>0) stop(paste0("The following required columns are missing:\n", paste0(missing_cols, collapse = ", ")))

  # check there are no duplicated sample IDs
  dup_samples <- mat[(col_index["OlinkID"]:col_index["LOD"]), 1] %>% duplicated()
  if(any(dup_samples)){
    dup_samples %<>% subset(mat[(col_index["OlinkID"]:col_index["LOD"]), 1], .) %>% unique()
    stop(paste0("Sample ID values must be unique (including controls). The following IDs are duplicated:\n",
                paste0(dup_samples, collapse = ", ")))
  }

  # assay data
  #   select rows beginning from and up to the indices associated with 'OlinkID' and 'LOD' in the first column
  #   select and set sample IDs as columns; select and set OlinkIDs as rownames
  assay_mat <- mat[(col_index["OlinkID"]):(col_index["LOD"]-1), which(is.na(mat[col_index["OlinkID"], ])==F)] %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames(var=colnames(.)[1]) %>%
    t() %>% as.data.frame() %>%
    dplyr::select(sort(colnames(.))) %>%
    dplyr::arrange(OlinkID) %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames(var="OlinkID") %>%
    mutate_all(as.matrix) %>%
    mutate_all(as.numeric)
  result$assay <- assay_mat

  # sample info
  #
  sample_info <- mat[(col_index["OlinkID"]+1):(col_index["LOD"]-1), which(is.na(mat[col_index["OlinkID"], ])==T)]
  # collapse Panel, Assay, and Uniprot ID columns to create sample column names
  tmp <- data.frame(unlist(mat[grep('Panel', mat[,1]), which(is.na(mat[col_index["OlinkID"], ])==T)]),
                    unlist(mat[grep('Assay', mat[,1]), which(is.na(mat[col_index["OlinkID"], ])==T)]),
                    unlist(mat[grep('Uniprot ID', mat[,1]), which(is.na(mat[col_index["OlinkID"], ])==T)]))
  names(sample_info) <- unlist(apply(tmp, 1, function(x) paste(na.omit(x), collapse = "--")))
  sample_info$sample_id <- mat[(col_index["OlinkID"]+1):(col_index["LOD"]-1), 1]
  rownames(sample_info) <- sample_info$sample_id
  result$sampleinfo <- sample_info %>% dplyr::arrange(sample_id) %>% dplyr::select(sample_id, tidyselect::everything())

  # protein info
  analyte_info <- mat %>% dplyr::slice(col_index) %>% column_to_rownames(colnames(mat)[1]) %>% t() %>%
    as.data.frame() %>% dplyr::filter(!is.na(OlinkID))
  colnames(analyte_info) <- ANA_COL_NAMES[colnames(analyte_info)]
  analyte_info %<>%
    dplyr::mutate(feature_id = OlinkID) %>%
    dplyr::mutate(BIOCHEMICAL = name) %>%
    dplyr::mutate(Panel_Version = sapply(Panel, function(x){regmatches(x, gregexpr("(?<=\\().*?(?=\\))", x, perl=T))[[1]]})) %>%
    dplyr::mutate(Panel = sapply(Panel, function(x){gsub("\\s*\\([^\\)]+\\)","",x)})) %>%
    dplyr::mutate(LOD = round(as.numeric(LOD), 5)) %>%
    dplyr::mutate(MissingFreq = round(as.numeric(MissingDataFreq), 2)) %>%
    dplyr::select(feature_id, LOD, MissingFreq, OlinkID, UniprotID, Panel, Panel_Version, BIOCHEMICAL, name)
  rownames(analyte_info) <- analyte_info$OlinkID
  result$metinfo <- analyte_info %>% dplyr::arrange(OlinkID)

  # generate summarized experiment
  D <- SummarizedExperiment(assay    = result$assay,
                            colData  = result$sampleinfo,
                            rowData  = result$metinfo,
                            metadata = list(sessionInfo=utils::sessionInfo()))

  # add original metadata if exists
  if (!is.null(result$meta$results)) metadata(D)$results <- result$meta$results
  if (!is.null(result$meta$settings)) metadata(D)$settings <- result$meta$settings

  # ensure colnames and rownames exist
  if (is.null(colnames(D))) colnames(D) <- 1:ncol(D)
  if (is.null(rownames(D))) rownames(D) <- result$metinfo$name

  # add status information
  funargs <- maplet:::mti_funargs()
  D %<>%
    maplet:::mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("Loaded Olink file: %s", basename(file))
    )

  # return
  D
}
