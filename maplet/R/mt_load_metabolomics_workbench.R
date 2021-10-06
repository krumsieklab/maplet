#' Loader for Metabolomics Workbench Studies
#'
#' Load SummarizedExperiment data frames using the mwtab file provided by studies on the metabolomics
#' data repository Metabolomics Workbench.\n\n Because of a formatting issue in the mwtab json objects
#'
#'
#' @param D \code{SummarizedExperiment} input. OPTIONAL.
#' @param study_id Identifier for the study of the form "ST######".
#' @param analysis_id Required if multiple datasets per study_id. Must be one of the form "AN######".
#' @param cache Use cached file to load data? Default: FALSE.
#' @param cache_file File for storing / retrieving cached data.
#'
#' @return If first step in pipeline, creates \code{SummarizedExperiment} object. Populates empty assay,
#' colData, and rowData data frames.
#'
#' @examples
#' \dontrun{D <-
#'   # load data
#'   D1 <- mt_load_metabolomics_workbench(study_id = "ST000991")
#'   # use analysis_id to distinguish between multiple datasets
#'   D2 <- mt_load_metabolomics_workbench(study_id = "ST001301", analysis_id = "AN002167")
#'   }
#'
#' @author KC
#'
#' @export
mt_load_metabolomics_workbench <- function(D, study_id, analysis_id, cache=FALSE, cache_file){

  # validate arguments
  if(missing(study_id)) stop("Value must be provided for argument \'study_id\'.")
  if(grepl("\\bST\\d+\\b", study_id)==FALSE) stop("Argument \'study_id\' must be of the form \'ST######\'.")
  if(!missing(analysis_id) && grepl("\\bAN\\d+\\b", analysis_id)==FALSE) stop("Argument \'analysis_id\' must be of the form \'AN######\'.")
  if(!missing(cache) && missing(cache_file)) stop("Value must be provided for \'cache_file\' when cache=TRUE.")

  # initialize result list
  result=list()

  # get metadata from D if present
  if(!missing(D)){
    # validate SE
    if ("SummarizedExperiment" %in% class(D) == F) stop("D is not of class SummarizedExperiment")
    if (length(assays(D))!=0) stop("Passed SummarizedExperiment assay must be empty!")

    # get metadata
    result$meta <- metadata(D)
  }

  # load cached file if available
  if(!missing(cache_file) && file.exists(cache_file)){
    mti_logstatus("Cached file found. Loading data from file.")
    D <- readRDS(cache_file)
    # download data from Metabolomics Workbench
  }else{
    # download data files
    mti_logstatus("Downloading data from Metabolomics Workbench...")
    # retrieve data from URL
    call <- glue::glue("https://www.metabolomicsworkbench.org/rest/study/study_id/{study_id}/mwtab")
    query_results <- httr::GET(call, timeout = 60, encode = "json")
    # if request not successful, crash
    call_status <- query_results$status_code
    if(call_status != 200) stop(glue::glue("Request error {call_status} for URL:\n{call}"))

    # convert call data to text and remove trailing garbage
    call_text <- httr::content(query_results,as="text", encoding = "UTF-8")
    call_text <- gsub("\\n", "", call_text)

    # try to parse mwtab json file
    #   if fails, call extract_json to find beginning and ending of objects
    #   before parsing
    study_list <- tryCatch({
      study_list <- list(jsonlite::fromJSON(call_text))
    },
    warning = function(cond){
      study_list <- extract_json(call_text)
      return(study_list)
    },
    error = function(cond){
      study_list <- extract_json(call_text)
      return(study_list)
    })

    if( length(study_list) > 1){
      # use analysis_type and ion_mode to determine which study to load
      if(missing(analysis_id)){
        stop("Multiple datasets detected. Argument \'analysis_id\' is required to select dataset to load.")
      }

      # get analysis IDs
      analysis_id_list <- study_list %>% purrr::map("METABOLOMICS WORKBENCH") %>% purrr::map("ANALYSIS_ID") %>% unlist
      analysis_idx <- which(analysis_id_list == analysis_id)
      if(length(analysis_idx) != 1) stop(glue::glue("Analysis ID {analysis_id} was not detected."))
      study <- study_list[[analysis_idx]]

    }else{
      study <- study_list[[1]]
    }

    data_idx <- which(endsWith(attributes(study)$names, "_METABOLITE_DATA"))
    if(length(data_idx) == 0){
      stop("No intensity data detected.")
    }else if(length(data_idx) > 1){
      stop("BUG: More than one intensity data frame detected for a single study!")
    }

    # extract data frames
    assay <- study[[data_idx]]$Data %>% tibble::column_to_rownames("Metabolite")
    rd <- study[[data_idx]]$Metabolites
    cd <- study$SUBJECT_SAMPLE_FACTORS %>% dplyr::bind_cols(study$SUBJECT)

    # check same number of metabolites / samples in rowData / colData as assay
    if(nrow(rd) > nrow(assay)) rd <- rd[rd$Metabolite %in% rownames(assay),]
    if(nrow(cd) > ncol(assay)) cd <- cd[cd$`Sample ID` %in% colnames(assay),]

    # assemble SummarizedExperiment object
    D <- SummarizedExperiment::SummarizedExperiment(assay = assay,
                                                    rowData = rd,
                                                    colData = cd,
                                                    metadata = list(sessionInfo=utils::sessionInfo()))

    # cache data
    if(cache) saveRDS(D, file = cache_file)

  }

  # add original metadata if exists
  if (!is.null(result$meta$results)) metadata(D)$results <- result$meta$results
  if (!is.null(result$meta$settings)) metadata(D)$settings <- result$meta$settings

  # ensure colnames and rownames exist
  if (is.null(colnames(D))) colnames(D) <- samp_names
  if (is.null(rownames(D))) rownames(D) <- rd$metabolite_identification

  # add status information
  logtxt = glue::glue("Loaded Data from Metabolomics Workbench for Study ID: {study_id}")
  if(!missing(analysis_id)) log_txt <- glue::glue(logtxt, ", Analysis ID: {analysis_id}")
  funargs <- mti_funargs()
  D %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = logtxt
    )

  # return
  D


}

# helper function - mwtab json pages not formatted correctly on Metabolomics Workbench REST pages.
#   Need to find beginning and end of the json object(s) before they can be parsed.
#   Returns all json objects found in mwtab json page.
extract_json <- function(call_text){
  start_idx_list <- gregexpr("\\{\"METABOLOMICS WORKBENCH", call_text)[[1]]
  num_studies <- length(start_idx_list)
  end_idx_list <- c(start_idx_list[2:num_studies]-1, nchar(call_text))
  study_list <- list()
  for(i in seq_len(num_studies)){
    start_idx = start_idx_list[i]
    end_idx = end_idx_list[i]
    tmp_str <- substr(call_text, start_idx, end_idx)
    study_list[[i]] <- jsonlite::fromJSON(tmp_str)
  }
  return(study_list)
}
