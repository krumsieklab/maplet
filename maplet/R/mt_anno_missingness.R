#' Add missingness annotation column
#'
#' Adds a rowData or colData column representing the missingness of features or samples, respectively.
#'
#' @param D \code{SummarizedExperiment} input.
#' @param anno_type Either "samples" (colData) or "features" (rowData).
#' @param out_col Name of new missingness column. Default: "missingness".
#'
#' @return rowData or colData: New annotation column added.
#'
#' @examples
#' \dontrun{
#'    # load data
#'    mt_load_metabolon(file=file, sheet="data") %>%
#'    # add feature missingness column to rowData
#'    mt_anno_missingness(anno_type = "features", out_col = "Feat_Missing")
#'    # add sample missingness column to colData
#'    mt_anno_missingness(anno_type = "samples", out_col = "Samp_Missing")
#' }
#'
#' @author KC
#'
#' @export
mt_anno_missingness <- function(D, anno_type, out_col="missingness"){

  # helper function
  # form mt_plots_qc_missingness - should be moved to mt_internal_helpers
  missingness <- function(X)apply(is.na(X),2,sum)/dim(X)[1]

  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))

  if (!(anno_type %in% c("samples","features"))) stop("anno_type must be either 'samples' or 'features'")
  if(length(anno_type) > 1) stop("anno_type must be either 'samples' OR 'features' - NOT both")

  X <- t(assay(D))

  if(anno_type == "features"){
    # get missingness of features
    miss_col <- missingness(X)

    # add missingness column to rowData
    rowData(D)[[out_col]] <- miss_col

    logtxt = sprintf("Added missingness annotation for %i features", length(miss_col))

  }else if(anno_type == "samples"){
    # get missingness of samples
    miss_col <- missingness(t(X))

    if(missing(out_col)){
      out_col <- "missingness"
    }
    colData(D)[[out_col]] <- miss_col

    logtxt = sprintf("Added missingness annotation for %i samples", length(miss_col))

  }else{
    stop("Bug")
  }

  # add status information
  funargs <- mti_funargs()
  D %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = logtxt
    )

  # return
  D

}
