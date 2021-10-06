#' Load Metabolon-format data (old format)
#'
#' For Metabolon-format version used prior to 2020-10-29. Loads data from a Metabolon-format Excel file. Needs to be in the original
#' "Client Data Table" format that they deliver.
#'
#' @param D \code{SummarizedExperiment} input. Missing if first step in pipeline.
#' @param file Name of input Excel file.
#' @param sheet Name of sheet.
#' @param copy_nan_sheet If given, which sheet to copy the NA pattern from. Default: '' (i.e. None).
#'
#' @return If first step in pipeline, creates \code{SummarizedExperiment} object. Populates empty assay, colData, and rowData data
#' frames.
#'
#' @examples
#' \dontrun{D <-
#'   # load data
#'   mt_load_metabolon_v1(file=system.file("extdata", "example_data/sampledata.xlsx", package = "maplet"),
#'                        sheet="OrigScale") %>%
#'   ...}
#'
#' @author JK
#'
#' @export
mt_load_metabolon_v1 <- function(D, file, sheet, copy_nan_sheet='') {

  # initialize result list
  result=list()

  # validate arguments
  if (missing(file)) stop("file must be provided")
  if (missing(sheet)) stop("sheet must be provided")

  # get metadata from D if present
  if(!missing(D)){
    # validate SE
    if ("SummarizedExperiment" %in% class(D) == F) stop("D is not of class SummarizedExperiment")
    if (length(assays(D))!=0) stop("Passed SummarizedExperiment assay must be empty!")

    # get metadata
    result$meta <- metadata(D)
  }

  # using readxl package:
  raw = readxl::read_excel(path=file, sheet=sheet, col_names = F)

  # find metabolite header row and last metabolite row
  imetheader = min(which(!is.na(raw[,1])))
  imetlast = max(which(apply(is.na(raw),1,sum)<dim(raw)[2]))
  # find sample header column and last sample row
  isampheader = min(which(!is.na(raw[1,])))
  isamplast = max(which(apply(is.na(raw),2,sum)<dim(raw)[1]))

  # fix overlapping cell
  overl=stringr::str_replace(gsub("\\s+", " ", stringr::str_trim(raw[imetheader,isampheader])), "B", "b")
  overl=strsplit(overl, " ")[[1]]
  overlmet = overl[2]
  overlsamp = overl[1]

  # extract metabolite information
  result$metinfo <- readxl::read_excel(path=file, sheet=sheet, col_names = T,
                                       range = readxl::cell_limits(ul = c(imetheader, 1),
                                                                   lr = c(imetlast  , isampheader)))
  result$metinfo <- as.data.frame(result$metinfo)

  # convert any spaces in the colnames to underscores
  colnames(result$metinfo) <- gsub(" ", "_", colnames(result$metinfo))
  # fix last one
  colnames(result$metinfo)[ncol(result$metinfo)] = overlmet
  # also add Metabolon ID in M00000 format
  result$metinfo$COMP_IDstr = sapply(sprintf('%05d',as.numeric(result$metinfo$COMP_ID)), function(x)paste0('M',x))

  # extract sample information
  result$sampleinfo = data.frame(t(raw[1:imetheader,(isampheader+1):isamplast]), stringsAsFactors = F)
  #colnames(result$sampleinfo) = as.list(raw[1:imetheader-1,isampheader])[[1]] # dirty hack, had something to do with the output format of read_excel
  colnames(result$sampleinfo) = c(as.vector(as.matrix(raw[1:imetheader-1,isampheader])),overlsamp) # dirty hack, had something to do with the output format of read_excel
  rownames(result$sampleinfo) = c()
  result$sampleinfo %<>% dplyr::mutate_all(readr::parse_guess)

  # extract data
  result$data <- t(raw[(imetheader+1):imetlast, (isampheader+1):isamplast])
  result$data <- as.data.frame(apply(result$data,2, as.numeric))

  # as column names, use "BIOCHEMICAL", if available
  if ("BIOCHEMICAL" %in% colnames(result$metinfo)) {
    colnames(result$data) = result$metinfo$BIOCHEMICAL
  } else {
    colnames(result$data) = c()
  }

  # as row names, use "SAMPLE_NAME", if available
  if ("SAMPLE_NAME" %in% colnames(result$sampleinfo)) {
    rownames(result$data) = result$sampleinfo$SAMPLE_NAME
  } else {
    rownames(result$data) = c()
  }

  # copy NanN from another sheet?
  if (nchar(copy_nan_sheet)>0) {
    # recursive call
    nandf = parseMetabolonFile(file=file, sheet=copy_nan_sheet)
    # sanity checks
    if (!(all.equal(colnames(result$data), colnames(nandf$data)))==T) {
      # figure out which are different
      l1 = colnames(result$data)
      l2 = colnames(nandf$data)
      dummy=sapply(1:length(l1), function(i){if(l1[i]!=l2[i]){sprintf('"%s" vs. "%s"\n',l1[i],l2[i])}})
      stop('some metabolites are different between data sheet and NaN sheet');
    }
    if (!all.equal(rownames(result$data), rownames(nandf$data)))
      stop('some sample names are different between data sheet and NaN sheet');
    # copy over NaNs
    result$data[is.na(nandf$data)]=NA

  }

  # return SummarizedExperiment

  # add display name
  result$metinfo$name   <- result$metinfo$BIOCHEMICAL
  # fix variable names
  colnames(result$data) <- result$metinfo$BIOCHEMICAL %>% make.names()

  # construct SummarizedExperiment
  D <- SummarizedExperiment(assay    = t(result$data),
                            colData  = result$sampleinfo,
                            rowData  = result$metinfo,
                            metadata = list(sessionInfo=utils::sessionInfo()))

  # add original metadata if exists
  if (!is.null(result$meta$results)) metadata(D)$results <- result$meta$results
  if (!is.null(result$meta$settings)) metadata(D)$settings <- result$meta$settings


  # ensure colnames and rownames exist
  if (is.null(colnames(D))) colnames(D) <- 1:ncol(D)
  if (is.null(rownames(D))) rownames(D) <- result$metinfo$BIOCHEMICAL

  # add status information
  funargs <- mti_funargs()
  D %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("loaded Metabolon file: %s, sheet: %s", basename(file), sheet)
    )

  # return
  D

}
