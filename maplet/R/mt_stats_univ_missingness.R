#' Perform missingness significance analysis
#'
#' This function will determine if NAs significantly accumulate in one of the sample groups. It is recommended that this function
#' is run without prior missing value filtering.
#'
#' @param D \code{SummarizedExperiment} input.
#' @param in_col Sample annotation (colData) column to compare against.
#' @param stat_name Name of comparison for later reference.
#' @param samp_filter Sample filter term to restrict to.
#'
#' @return $results$output: Statistics object.
#'
#' @examples
#' \dontrun{# run on sample field 'Group', name output stats object 'miss'
#' ... %>% mt_stats_univ_missingness(in_col = 'Group', stat_name='miss') %>% ...
#' }
#'
#' @author JK
#'
#' @export
mt_stats_univ_missingness <- function(D, in_col, stat_name, samp_filter) {

  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  stopifnot(length(in_col)==1)

  # merge data with sample info
  Ds <- D

  ## FILTER SAMPLES
  if(!missing(samp_filter)) {

    filter_q <- dplyr::enquo(samp_filter)
    num_samp <- ncol(Ds)
    samples.used <-  Ds %>% mti_format_se_samplewise() %>%
      mti_filter_samples(filter_q, num_samp)
    Ds <- Ds[,samples.used]

  } else {
    samples.used = rep(T, ncol(Ds))
  }

  # get variable to compare to
  if (!(in_col %in% colnames(colData(Ds)))) stop(sprintf("'%s' not found in sample annotations.", in_col))
  fixorder = function(x){o= unique(as.character(x)); gdata::reorder.factor(x, new.order=o)}
  vc = fixorder(as.factor(colData(Ds)[[in_col]]))
  if (length(levels(vc))<2) stop(sprintf("'%s' has less than 2 factor levels",in_col))




  # run models
  rawres <- sapply(1:nrow(Ds), function(i){
    # for (i in 1:nrow(Ds)) {

    # get feature
    m <- assay(Ds)[i,]

    # construct table
    tab <- table(is.na(m), vc)
    # fix table
    if (!("TRUE" %in% rownames(tab))) {
      tab %<>% rbind(rep(0,ncol(tab)))
      rownames(tab)[2] = "TRUE"
    }
    if (!("FALSE" %in% rownames(tab))) {
      tab %<>% rbind(rep(0,ncol(tab)))
      rownames(tab)[2] = "FALSE"
    }
    # run test
    test <- stats::fisher.test(tab)
    if (!("estimate" %in% names(test))) test$estimate=NA
    # construct extra fields
    ex <- gdata::unmatrix(tab)
    names(ex) = gsub("TRUE","missing",gsub("FALSE","present", names(ex)))
    # return
    as.data.frame(cbind(data.frame(var=rownames(Ds)[i], statistic=test$estimate, p.value=test$p.value),t(as.data.frame(ex))))

    # }
  },simplify=F)
  # combine
  res <- do.call(rbind, rawres)
  rownames(res) <- NULL

  ## add status information & results
  funargs <- mti_funargs()
  D %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("missingness analysis with variable %s", as.character(in_col)),
      output = list(
        table   = res,
        name    = stat_name,
        samples.used = samples.used,
        outcome = in_col
      )
    )

  ## return
  D

}



