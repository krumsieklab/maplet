#' Load Metabolon-format lipidomics data
#'
#' Loads lipidomics data from a Metabolon-format Excel file. Needs to be in the original "Client Data Table" format that they deliver.
#'
#' @param D \code{SummarizedExperiment} input. Missing if first step in pipeline.
#' @param file Name of input Excel file.
#' @param sheet_list List with the sheet names to read.
#'
#' @return If first step in pipeline, creates \code{SummarizedExperiment} object. Populates empty assay, colData, and rowData data frames.
#'
#' @examples
#' \dontrun{
#' D <-
#'   # load data
#'   mt_load_files_metabolon_lipidomics(file=system.file("extdata", "example_data/sampledata.xls", package = "maplet"),
#'                                      sheet_list = c("Species Concentrations","Fatty Acid Concentrations")) %>%
#'   ...}
#'
#' @author EB
#'
#' @import readxl
#' @import stringr
#' @import tidyverse
#'
#' @export
mt_load_metabolon_lipidomics <- function(D, file, sheet_list) {

  # initialize outer result list
  result <- list()

  # validate arguments
  if (missing(file)) stop("file must be provided")
  if (missing(sheet_list)) stop("sheet must be provided")

  # get metadata from D if present
  if(!missing(D)){
    # validate SE
    if ("SummarizedExperiment" %in% class(D) == F) stop("D is not of class SummarizedExperiment")
    if (length(assays(D))!=0) stop("Passed SummarizedExperiment assay must be empty!")

    # get metadata
    result$meta <- metadata(D)
  }


  # using readxl package:
  raw <- lapply(sheet_list %>% {names(.)=.;.}, function(x){
    readxl::read_excel(path=file, sheet=x, col_names = F)
  })

  xx <- lapply(sheet_list %>% {names(.)=.;.}, function(x){

    result <- list()

    # find metabolite header row and last metabolite row
    imetheader = which(!is.na(raw[[x]][,1]))[3]
    imetlast = max(which(apply(is.na(raw[[x]]),1,sum)<dim(raw[[x]])[2]))
    # find sample header column and last sample row
    isampheader = min(which(!is.na(raw[[x]][5,])))
    isamplast = max(which(apply(is.na(raw[[x]]),2,sum)<dim(raw[[x]])[1]))

    # fix overlapping cell
    overl=gsub("\\s+", " ", stringr::str_trim(raw[[x]][imetheader,isampheader]))
    overl=strsplit(overl, " ")[[1]]
    overlmet = overl[2]
    overlsamp = overl[1]

    # extract metabolite information
    result$metinfo <- readxl::read_excel(path=file, sheet=x, col_names = T,
                                         range = readxl::cell_limits(ul = c(imetheader+1, 1),
                                                                     lr = c(imetlast+1 , isampheader)))
    result$metinfo <- as.data.frame(result$metinfo)

    # convert any spaces in the colnames to underscores
    colnames(result$metinfo) <- gsub(" ", "_", colnames(result$metinfo))
    # fix last one
    colnames(result$metinfo)[ncol(result$metinfo)] = overlmet

    # extract sample information
    result$sampleinfo = data.frame(t(raw[[x]][5:imetheader,(isampheader+1):isamplast]), stringsAsFactors = F)
    #colnames(result$sampleinfo) = as.list(raw[1:imetheader-1,isampheader])[[1]] # dirty hack, had something to do with the output format of read_excel
    colnames(result$sampleinfo) = c(as.vector(as.matrix(raw[[x]][6:imetheader-1,isampheader])),overlsamp) # dirty hack, had something to do with the output format of read_excel
    rownames(result$sampleinfo) = c()
    result$sampleinfo %<>% dplyr::mutate_all(parse_guess)
    # convert any spaces in the colnames to underscores
    colnames(result$sampleinfo) <- gsub(" ", "_", colnames(result$sampleinfo))

    # extract data
    result$data <- t(raw[[x]][(imetheader+1):imetlast, (isampheader+1):isamplast])
    result$data <- as.data.frame(apply(result$data,2, as.numeric))

    # as column names, use "Name", if available
    if ("Name" %in% colnames(result$metinfo)) {
      colnames(result$data) = result$metinfo$Name
    } else {
      colnames(result$data) = c()
    }

    # as row names, use "CLIENT_IDENTIFIER", if available
    if ("CLIENT_IDENTIFIER" %in% colnames(result$sampleinfo)) {
      rownames(result$data) = result$sampleinfo$CLIENT_IDENTIFIER
    } else {
      rownames(result$data) = c()
    }
    # add extra column for later merging
    result$data %<>% as.data.frame
    result$data$mergeby <- rownames(result$data)

    # adding variable for matching later on
    result$sampleinfo$id = rownames(result$data)

    # add display name
    result$metinfo$name <- result$metinfo$Name

    result

  })

  ## merge data and annotations from the different data sheets
  # join all data
  result$data <- xx[[sheet_list[1]]]$data
  if(length(sheet_list)>1) {
    for (i in 2:length(sheet_list)) {
      # throw an error if the sample names are different across data sheets
      if (any(xx[[sheet_list[1]]]$data$mergeby!=xx[[sheet_list[i]]]$data$mergeby)) {
        stop(sprintf("Sample names of %s and %s are not the same.", sheet_list[1], sheet_list[i]))
      }
      result$data <- result$data %>%
        dplyr::full_join(xx[[sheet_list[i]]]$data, by = "mergeby")
    }
  }
  rownames(result$data) <- result$data$mergeby
  result$data$mergeby <- NULL

  # join all sample annotations
  result$sampleinfo <- xx[[sheet_list[1]]]$sampleinfo
  if(length(sheet_list)>1) {
    for (i in 2:length(sheet_list)) {
      # throw a warning if sample annotations are different across data sheets
      if(any(colnames(xx[[sheet_list[1]]]$sampleinfo)!=colnames(xx[[sheet_list[i]]]$sampleinfo))) {
        mti_logwarning(sprintf("The sample annotations of %s and %s are not the same. Sample annotations from %s used.", sheet_list[1], sheet_list[i], sheet_list[1]))
      }
    }
  }

  # join all metabolite annotations
  result$metinfo <- xx[[sheet_list[1]]]$metinfo
  if(length(sheet_list)>1) {
    for (i in 2:length(sheet_list)) {
      result$metinfo <- result$metinfo %>%
        rbind(xx[[sheet_list[i]]]$metinfo)
    }
  }

  # sort metabolites and samples in annotations to match order in data
  result$sampleinfo <- result$sampleinfo[match(rownames(result$data),result$sampleinfo$id),]
  result$metinfo <- result$metinfo[match(colnames(result$data), result$metinfo$name),]

  # rownames of assay (right now colnames) must be R-compatible names
  colnames(result$data) <- make.names(colnames(result$data))

  # construct SummarizedExperiment
  D <- SummarizedExperiment(assay    = t(result$data),
                            colData  = result$sampleinfo,
                            rowData  = result$metinfo,
                            metadata = list(sessionInfo=utils::sessionInfo()))

  # add original metadata if exists
  if (!is.null(result$meta$results)) metadata(D)$results <- result$meta$results
  if (!is.null(result$meta$settings)) metadata(D)$settings <- result$meta$settings

  # add status information
  funargs <- mti_funargs()
  D %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("loaded Metabolon lipidomics file: %s, sheets: %s", basename(file), paste(sheet_list, collapse = ", "))
    )

  # return
  D
}
