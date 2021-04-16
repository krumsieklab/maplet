#' Filter by missingness
#'
#' Filters either samples or features by fraction of missing values. Won't do both in one call.
#'
#' @param D \code{SummarizedExperiment} input.
#' @param feat_max Maximum fraction of missing features (between 0 and 1.0). Default: NA.
#' @param samp_max Maximum fraction of missing samples (between 0 and 1.0). Default: NA.
#' @param group_col Name of of a colData sample annotation column; feat_max will be applied to each group separately,
#'    feature must have at most feat_max in any of the groups. Default: NA.
#' @param report_filtered Write list of filtered features into the status log text? Default: F.
#' @param report_sample_col Required if report_filtered=T for sample filtering. Specifies which column of colData() to
#'    output into the status log. Default: ''.
#'
#' @return assay: Rows or columns will be filtered.
#' @return $results$output: Logical vector of kept features/samples.
#'
#' @examples
#' \dontrun{# first remove samples with >10% missingness, then features with >20% missingness
#' ... %>%
#'   mt_pre_filter_missingness(samp_max=0.1) %>%
#'   mt_pre_filter_missingness(feat_max=0.2) %>%
#' ...}
#'
#' @author JK
#'
#' @export
mt_pre_filter_missingness <- function(D,
                                      feat_max=NA,
                                      samp_max=NA,
                                      group_col = NA,
                                      report_filtered = F,
                                      report_sample_col='') {

  # check arguments, SummarizedExperiment, and exactly one cutoff argument must be non-NA
  stopifnot("SummarizedExperiment" %in% class(D))
  stopifnot(!is.na(feat_max) || !is.na(samp_max))
  stopifnot(!(is.na(feat_max) && is.na(samp_max)))
  stopifnot(!(!is.na(feat_max) && (feat_max<0 || feat_max>1)))
  stopifnot(!(!is.na(samp_max) && (samp_max<0 || samp_max>1)))

  # perform filtering
  if (!is.na(feat_max)) {
    NA.mat = is.na(assay(D))
    # feature
    lst_filtered = c()
    if(is.na(group_col)){ # if group feat_max
      na.stat = rowSums(NA.mat)
      metKeep = na.stat/ncol(D) <= feat_max
      # store list of filtered
      lst_filtered = D %>% rowData() %>% .$name %>% .[!metKeep]
      # filter
      D=D[metKeep,]
      na.stat[metKeep]
    }else{
      stopifnot(group_col %in% (colData(D) %>% colnames))
      xgroup_col = colData(D)[,group_col]
      na.stat = xgroup_col %>% unique %>% {. = as.character(.); names(.) = .; .} %>%
        sapply(function(x) rowSums(NA.mat[, xgroup_col == x])/sum(xgroup_col == x))
      metKeep = rowSums( na.stat <= feat_max ) == ncol(na.stat)
      # store list of filtered
      lst_filtered = D %>% rowData() %>% .$name %>% .[!metKeep]
      # filter
      D=D[metKeep,]
      na.stat[metKeep, ]
    }

    # generate logtxt
    logtxt = sprintf('features filtered, %.2f%%, %d of %d removed', round(feat_max*100),sum(!metKeep),length(metKeep))
    if (report_filtered) {
      logtxt %<>% paste0(sprintf("  -  list of filtered features: %s", paste0(lst_filtered, collapse=" / ")))
    }

    # add status information
    funargs <- mti_funargs()
    D %<>% 
      mti_generate_result(
        funargs = funargs,
        logtxt = logtxt,
        output = list(kept=as.vector(metKeep), na.stat = na.stat, na.mat = NA.mat[metKeep,])
      )

  } else {

    # sample
    sampleKeep = apply(is.na(assay(D)),2,sum)/nrow(D) <= samp_max

    # generate logtxt
    logtxt = sprintf('samples filtered, %.2f%%, %d of %d removed', samp_max*100,sum(!sampleKeep),length(sampleKeep))
    if (report_filtered) {
      if (nchar(report_sample_col)==0) stop("If report_filtered=T and filtering samples, report_sample_col must be specified.")
      lst_filtered = D %>% colData() %>% .[[report_sample_col]] %>% .[!sampleKeep]
      logtxt %<>% paste0(sprintf("  -  list of filtered features: %s", paste0(lst_filtered, collapse=" / ")))
    }

    # filter
    D=D[,sampleKeep]

    # add status information
    funargs <- mti_funargs()
    D %<>% 
      mti_generate_result(
        funargs = funargs,
        logtxt = logtxt,
        output = list(kept=as.vector(sampleKeep))
      )

  }

  # throw error if filtering caused empty dataset
  if (ncol(D)==0 || nrow(D)==0) stop("Filtering resulted in empty dataset.")

  # return
  D
}


