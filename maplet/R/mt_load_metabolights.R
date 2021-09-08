#' Loader for MetaboLights Studies
#'
#' Load SummarizedExperiment data frames using the metabolighteR package to retrieve publicly available
#' data files from the metabolomics data repository Metabolights. If cache=TRUE, the assembled
#' SummarizedExperiment object (with variable name "D_cache") will be saved to a local directory. This loader
#' will always check for the presence of a cached file corresponding to the study ID and ion mode (optional)
#' provided. If found, this file will be used and new data will not be downloaded.
#'
#' @param D \code{SummarizedExperiment} input. OPTIONAL.
#' @param study_id Identifier for the study of the form "MTBLS#" (where # is a number with 1 or more digits).
#' @param met_file Metabolite table file (prefix "m_") to download. Separates into assay and rowData data frames.
#'    Use metabolighteR::get_study_files(study_id) for a list of files available.
#' @param samp_file Sample table file to download (prefix "s_"). Data for colData data frame. Use
#'    metabolighteR::get_study_files(study_id) for a list of files available.
#' @param samp_col Name of column containing sample names in sample table. Default: "Sample.Name".
#' @param samp_suffix OPTIONAL. If sample names in metabolite table differ from names in samp_col by a suffix,
#'    provide the suffix.
#' @param cache Cache downloaded data to local directory? Default: FALSE.
#' @param cache_file File for storing / retrieving cached data.
#'
#' @return If first step in pipeline, creates \code{SummarizedExperiment} object. Populates empty assay,
#' colData, and rowData data frames.
#'
#' @examples
#' \dontrun{
#'   # load data
#'   D1 <-mt_load_metabolights(study_id = "MTBLS301",
#'                             met_file = "m_MTBLS301_LC-MS__metabolite_profiling_v2_maf.tsv",
#'                             samp_file = "s_MTBLS301.txt")
#'   D2 <- mt_load_metabolights(study_id = "MTBLS264",
#'                              met_file = "m_mtbls264_POS_mass_spectrometry_v2_maf.tsv",
#'                              samp_file = "s_MTBLS264.txt",
#'                              samp_suffix="_POS")}
#'
#' @author KC
#'
#' @export
mt_load_metabolights <- function(D,
                                 study_id,
                                 met_file,
                                 samp_file,
                                 samp_col="Sample.Name",
                                 samp_suffix,
                                 cache=FALSE,
                                 cache_file){

  # validate arguments
  if(missing(study_id)) stop("Value must be provided for argument \'study_id\'.")
  if(grepl("\\bMTBLS\\d+\\b", study_id)==FALSE) stop("Argument \'study_id\' must be of the form \'MTBLS#\'.")
  if(cache==TRUE && missing(cache_file)) stop("Argument \'cache_file\' must be provided when cache=TRUE.")

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
    maplet:::mti_logstatus("Cached file found. Loading data from file.")
    D <- readRDS(cache_file)
  # download data from Metabolights
  }else{
    # check met_file and samp_file provided
    studyFileList <- metabolighteR::get_study_files(study_id)
    if(missing(met_file) | missing(samp_file)) stop("To download data, values must be provided for both \'met_file\' and \'samp_file\'.")
    if(samp_file %in% studyFileList$file == F) stop(glue::glue("Unrecognized sample file name for study {study_id}: {samp_file}"))
    if(met_file %in% studyFileList$file == F) stop(glue::glue("Unrecognized metabolite file name for study {study_id}: {met_file}."))

    # download data files
    maplet:::mti_logstatus("Downloading data from MetaboLights...")
    # metabolite table - assay + rowData
    met_tab <- metabolighteR::download_study_file(study_id, met_file) %>% suppressMessages()

    # sample table - colData
    cd <- metabolighteR::download_study_file(study_id, samp_file) %>% suppressMessages()
    if(samp_col %in% colnames(cd)==F) stop(glue::glue("Column {samp_col} not found in sample table. Please check sample table file."))

    # separate metabolite table into assay and rowData
    samp_col_names <- cd[[samp_col]] %>% make.names()
    samp_names <- samp_col_names %>% intersect(colnames(met_tab))

    # if names not found, try adding a suffix
    if(length(samp_names)==0){
      if(!missing(samp_suffix)){
        samp_col_names <- cd[[samp_col]] %>% make.names() %>% paste0(samp_suffix)
        samp_names <- samp_col_names %>% intersect(colnames(met_tab))
        if(length(samp_names)==0) stop("Sample names don't match. Check sample names in file.")
      }else{
        stop("Sample names don't match. Check sample names in file.")
      }
    }

    assay <- met_tab %>% dplyr::select(samp_names)
    rd <- met_tab %>% dplyr::select(-samp_names)
    cd <- cd[samp_col_names %in% samp_names,]

    # create SummarizedExperiment object
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
  logtxt = glue::glue("Loaded Data from MetaboLights for Study ID: {study_id}")
  funargs <- maplet:::mti_funargs()
  D %<>%
    maplet:::mti_generate_result(
      funargs = funargs,
      logtxt = logtxt
    )

  # return
  D

}
