#' Load Olink-format data
#'
#' @description
#' Loads data from an Olink-format Excel file.
#' Uses Olink's R code read_NPX taken from
#' \href{https://github.com/Olink-Proteomics/OlinkRPackage/tree/master/OlinkAnalyze}{https://github.com/Olink-Proteomics/Olink
#' RPackage/tree/master/OlinkAnalyze}.
#'
#' @description
#' In case the Olink file is not in XLSX format, but CSV or TSV text format:
#' \itemize{
#'    \item you may need to remove all Ctrl data columns
#'    \item save file as xlsx (make sure to rename SEPT9 back from Excel date)
#'    \item don't keep the Excel file open while running R - this throws a file access denied error
#' }
#'
#' @param D \code{SummarizedExperiment} input. Missing if first step in pipeline.
#' @param file Name of NPX file exported from NPX manger.
#'
#' @return If first step in pipeline, creates \code{SummarizedExperiment} object. Populates empty assay (note: 2^NPX), colData, and rowData data frames.
#'
#' @examples
#' \dontrun{D <-
#'   # load data
#'   mt_load_olink(file=system.file("extdata", "example_data/sampledata.xlsx", package = "maplet")) %>%
#'   ...}
#'
#' @author KS
#'
#' @import readxl
#' @import tidyverse
#'
#' @export
mt_load_olink <- function(D, file){

  # initialize outer result list
  result <- list()

  # validate arguments
  if (missing(file)) stop("file must be provided")

  # save input information
  result$info$file <- file

  # get metadata from D if present
  if(!missing(D)){
    # validate SE
    if ("SummarizedExperiment" %in% class(D) == F) stop("D is not of class SummarizedExperiment")
    if (length(assays(D))!=0) stop("Passed SummarizedExperiment assay must be empty!")

    # get metadata
    result$meta <- metadata(D)
  }

  # read the Olink file
  odata = read_NPX(file)

  # convert to wide format
  wdata = reshape2::dcast(odata,formula = " SampleID ~ OlinkID", value.var = "NPX")
  xdata = simplify2array(wdata[,-1])
  rownames(xdata) = wdata$SampleID

  # convert to SummarizedExperiment
  # the assay is called "exprs" to be compatible with autonomics
  D = SummarizedExperiment(assay = list(exprs = t(xdata)),
                           colData = DataFrame(sample_id = rownames(xdata)),
                           rowData = DataFrame(feature_id = colnames(xdata),
                                               LODdata = reshape2::dcast(odata,formula = " SampleID ~ OlinkID", value.var = "LOD")[1,-1] %>% unlist() %>% as.numeric() %>% unname(),
                                               MissingFreq = reshape2::dcast(odata,formula = " SampleID ~ OlinkID", value.var = "MissingFreq")[1,-1] %>% unlist() %>% as.numeric() %>% unname(),
                                               OlinkID = reshape2::dcast(odata,formula = " SampleID ~ OlinkID", value.var = "OlinkID")[1,-1] %>% unlist() %>% unname(),
                                               UniProt = reshape2::dcast(odata,formula = " SampleID ~ OlinkID", value.var = "UniProt")[1,-1] %>% unlist() %>% unname(),
                                               Assay = reshape2::dcast(odata,formula = " SampleID ~ OlinkID", value.var = "Assay")[1,-1] %>% unlist() %>% unname(),
                                               Panel = reshape2::dcast(odata,formula = " SampleID ~ OlinkID", value.var = "Panel")[1,-1] %>% unlist() %>% unname(),
                                               Panel_Version = reshape2::dcast(odata,formula = " SampleID ~ OlinkID", value.var = "Panel_Version")[1,-1] %>% unlist() %>% unname(),
                                               PlateID = reshape2::dcast(odata,formula = " SampleID ~ OlinkID", value.var = "PlateID")[1,-1] %>% unlist() %>% unname()
                           ),
                           metadata = list(sessionInfo=utils::sessionInfo(), parseInfo=result$info))

  # add original metadata if exists
  if (!is.null(result$meta$results)) metadata(D)$results <- result$meta$results
  if (!is.null(result$meta$settings)) metadata(D)$settings <- result$meta$settings

  # do some sanity checks
  stopifnot (length(which(colnames(xdata) != rowData(D)$OlinkID)) == 0)
  MissingFreq2 = reshape2::dcast(odata,formula = " SampleID ~ OlinkID", value.var = "MissingFreq")[2,-1] %>% unlist() %>% as.numeric() %>% unname()
  stopifnot (length(which(MissingFreq2 != rowData(D)$MissingFreq)) == 0)

  # as column names, use "BIOCHEMICAL" for autonomics and "name" for maplet
  rowData(D)$BIOCHEMICAL = rowData(D)$Assay
  rowData(D)$name = rowData(D)$Assay

  # as row names, use "SAMPLE_NAME", if available
  D$SAMPLE_NAME = D$sample_id

  # add status information
  funargs <- mti_funargs()
  D %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("loaded Olink file: %s", basename(file))
    )

  # return
  D

}
###########################################################################
# KS: note - the function read_NPX has been copied out of the OlinkAnalyze package.
#     it coud also be obtained by installing the OlinkAnalyze package instead, by using the following:
#       devtools::install_github(repo ='Olink-Proteomics/OlinkRPackage/OlinkAnalyze')
#     however, it appears that installing that package creates problems with maplet,
#     so for now we just copy the function here
###########################################################################
read_NPX <- function (filename)
{
  NORM_FLAG <- F
  meta_dat <- readxl::read_excel(
    filename,
    skip = 2,
    n_max = 4,
    col_names = F,
    .name_repair = "minimal"
  )
  meta_dat[4, 1] <- "SampleID"
  NR_DEVIATIONS <-
    sum(stringr::str_detect(meta_dat[2,], "QC Deviation from median"))
  nr_col <- ncol(meta_dat)
  names(meta_dat) <- as.character(1:nr_col)
  meta_dat <- meta_dat %>% rename(Name = `1`)
  dat <- readxl::read_excel(
    filename,
    skip = 6,
    col_names = F,
    .name_repair = "minimal",
    col_types = c("text")
  )
  nr_col <- ncol(dat)
  names(dat) <- as.character(1:nr_col)
  dat <- dat %>% rename(Name = `1`)
  missfreq <-
    dat %>% filter(stringr::str_detect(Name, "Missing Data freq."))
  LOD <- dat %>% filter(stringr::str_detect(Name, "LOD"))
  norm_method <-
    dat %>% filter(stringr::str_detect(Name, "Normalization"))
  if (nrow(norm_method) == 0) {
    dat <- dat[c(-1 * (nrow(dat) - 2):nrow(dat)),]
  }
  else {
    dat <- dat[c(-1 * (nrow(dat) - 3):nrow(dat)),]
    NORM_FLAG <- T
  }
  meta_dat <- rbind(meta_dat, missfreq, LOD, norm_method)
  nr_panel <- (ncol(meta_dat) - 1 - NR_DEVIATIONS) / 94
  SampleID <- dat$Name
  Index_nr <- c(1:length(SampleID))
  panel_data <- list()
  QC_list <- list()
  meta_data_list <- list()
  panel_list <- list()
  assay_name_list <- list()
  panel_list_long <- list()
  deviations_list <- list()
  for (i in 1:nr_panel) {
    panel_data[[i]] <- dat[, (2 + ((i - 1) * 92)):(93 + ((i -
                                                            1) * 92))]
    if (NR_DEVIATIONS == 0) {
      QC_list[[i]] <- dat[, c((2 + ((nr_panel) * 92) +
                                 (i - 1)), (2 + ((nr_panel) * 92) + (i - 1)) +
                                nr_panel)]
      meta_data_list[[i]] <- meta_dat[, c((2 + ((i - 1) *
                                                  92)):(93 + ((i - 1) * 92)), (2 + ((nr_panel) *
                                                                                      92) + (i - 1)), (2 + ((nr_panel) * 92) + (i -
                                                                                                                                  1)) + nr_panel)]
    }
    else {
      QC_list[[i]] <- dat[, c((2 + ((nr_panel) * 92) +
                                 (i - 1)),
                              (2 + ((nr_panel) * 92) + (i - 1)) +
                                nr_panel,
                              (2 + ((nr_panel) * 92) + (i - 1)) +
                                2 * nr_panel,
                              (2 + ((nr_panel) * 92) + (i - 1)) +
                                3 * nr_panel)]
      meta_data_list[[i]] <- meta_dat[, c((2 + ((i - 1) *
                                                  92)):(93 + ((i - 1) * 92)),
                                          (2 + ((nr_panel) *
                                                  92) + (i - 1)),
                                          (2 + ((nr_panel) * 92) + (i -
                                                                      1)) + nr_panel,
                                          (2 + ((nr_panel) * 92) + (i -
                                                                      1)) + 2 * nr_panel,
                                          (2 + ((nr_panel) * 92) +
                                             (i - 1)) + 3 * nr_panel
      )]
      meta_data_list[[i]][4, 95] <- "QC Deviation Inc Ctrl"
      meta_data_list[[i]][4, 96] <- "QC Deviation Det Ctrl"
    }
    meta_data_list[[i]][4, 93] <- meta_data_list[[i]][2,
                                                      93]
    meta_data_list[[i]][4, 94] <- meta_data_list[[i]][2,
                                                      94]
    panel_list[[i]] <- cbind(panel_data[[i]], QC_list[[i]])
    colnames(panel_list[[i]]) <- unlist(meta_data_list[[i]][4,])
    panel_list[[i]] <-
      panel_list[[i]][,!is.na(stringr::str_detect(colnames(panel_list[[i]]),
                                                  "SampleID|OID[0-9]{5}"))]
    panel_list[[i]][, c(-93,-94)] <- lapply(panel_list[[i]][,
                                                            c(-93,-94)], function(x)
                                                              as.numeric(stringr::str_replace_all(x,
                                                                                                  c(
                                                                                                    `#` = "",
                                                                                                    `,` = ".",
                                                                                                    `No Data` = NA
                                                                                                  ))))
    assay_name_list[[i]] <- tibble(
      ID = c(t(meta_data_list[[i]][4,])),
      Name = c(t(meta_data_list[[i]][2,])),
      UniProt = c(t(meta_data_list[[i]][3,])),
      Panel = c(t(meta_data_list[[i]][1,])),
      MissingFreq = c(t(meta_data_list[[i]][5,])),
      LOD = as.numeric(c(t(
        meta_data_list[[i]][6,]
      )))
    )
    if (NORM_FLAG == T) {
      assay_name_list[[i]] <- bind_cols(assay_name_list[[i]],
                                        Normalization = c(t(meta_data_list[[i]][7,])))
    }
    panel_list_long[[i]] <-
      panel_list[[i]] %>% mutate(SampleID = SampleID) %>%
      mutate(Index = Index_nr) %>% gather(
        Assay,
        NPX,
        -SampleID,-`QC Warning`,
        -`Plate ID`,
        -Index,
        -matches("*Inc Ctrl*"),-matches("*Det Ctrl*")
      ) %>% left_join(assay_name_list[[i]],
                      by = c(Assay = "ID")) %>% select(
                        SampleID,
                        Index,
                        Assay,
                        UniProt,
                        Name,
                        MissingFreq,
                        Panel,
                        `Plate ID`,
                        `QC Warning`,
                        LOD,
                        NPX,
                        matches("Normalization"),
                        matches("*Inc Ctrl*"),
                        matches("*Det Ctrl*")
                      ) %>%
      rename(PlateID = `Plate ID`) %>% rename(QC_Warning = `QC Warning`) %>%
      rename(OlinkID = Assay, Assay = Name)
  }
  bind_rows(panel_list_long) %>% filter(!is.na(SampleID)) %>%
    tbl_df %>% mutate(Panel_Version = gsub(".*\\(", "", Panel)) %>%
    mutate(Panel_Version = gsub("\\)", "", Panel_Version)) %>%
    mutate(Panel = gsub("\\(.*\\)", "", Panel)) %>% select(
      SampleID,
      Index,
      OlinkID,
      UniProt,
      Assay,
      MissingFreq,
      Panel,
      Panel_Version,
      PlateID,
      QC_Warning,
      LOD,
      NPX,
      matches("Normalization"),
      matches("*Inc Ctrl*"),
      matches("*Det Ctrl*")
    )
}
