#' Composite score models
#'
#' Compute composite scores for outcome ~ features.
#'
#' @param D \code{SummarizedExperiment} input.
#' @param in_col Name of the column in colData to compare with features.
#' @param id_col Name of the column in colData having patient_ID.
#' @param stat_name Name under which this comparison will be stored, must be unique to all other statistical results.
#' @param samp_filter Term which samples to filter to first... e.g. used if the data contains >2 groups but the user wants to run a two-group comparison.
#'
#' @return $result$output: statistics object
#'
#' @import survival
#' @import MatrixEQTL
#'
#' @examples
#' \donttest{
#'  mt_stats_univ_cs(
#'    in_col     = Group,
#'    id_col = Subject_ID,
#'    stat_name         = "composite_score",
#'    samp_filter = (Group %in% c("grp1","grp2")),
#'  ) %>% ...
#'  }
#'
#' @author RB
#'
#' @export
mt_stats_univ_cs <- function(D,
                             in_col,
                             id_col,
                             stat_name,
                             samp_filter) {

  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))

  # make sure name does not exist yet
  if (stat_name %in% unlist(maplet::mtm_res_get_entries(D, "stats") %>% purrr::map("output") %>% purrr::map("name"))) stop(sprintf("stat element with name '%s' already exists",stat_name))

  # merge data with sample info
  Ds <- D %>% maplet:::mti_format_se_samplewise() # NOTE: No explosion of dataset size, no gather() - 8/17/20, JK

  ## FILTER SAMPLES
  if(!missing(samp_filter)) {

    filter_q <- dplyr::enquo(samp_filter)
    num_samp <- ncol(D)
    samples.used <- maplet:::mti_filter_samples(Ds, filter_q, num_samp)
    Ds <- Ds[samples.used,]

  } else {
    samples.used = rep(T, ncol(D))
  }
  # check that outcome is either binary or numerical and run test accordingly
  feats <- rownames(D)
  outvec <- Ds[[in_col]]
  cl <- outvec %>% class()

  if (("character" %in% cl) || ("factor" %in% cl)) {
    stop("Composite scores must be numeric!")
  } else {
    outgroups <- NULL
    # run concordance
    # we want to model feature ~ outcome
    cstest <- lapply(feats, function(x){
      co <- survival::concordance(Ds[[in_col]] ~ Ds[[x]] + cluster(Ds[[id_col]]))
      z <- (co$concordance-0.5)/sqrt(co$var)
      p <- -expm1(pnorm(abs(z),lower.tail = T,log.p = T ))
      stat <- unname(co$concordance)
      res <- data.frame("estimate"=stat, "statistic"=z, "p.value"=p, "method"='concordance')
      return(res)
    }) %>% do.call(rbind,.) %>% data.frame()
  }

  # add columns with feature names and y variable
  cstest %<>% mutate(var=feats, term=in_col)
  # reorder columns
  cstest %<>% select(var, term, estimate, statistic, p.value)
  # rearrange back to original feature order
  o <- match(rownames(D),cstest$var)
  stopifnot(!any(is.na(o))) # sanity check
  cstest <- cstest[o,]
  # make sure that NAs in the outcome are set to FALSE in the list of used samples
  samples.used[is.na(Ds[[in_col]])] <- F

  ## add status information & results
  funargs <- maplet:::mti_funargs()
  metadata(D)$results %<>%
    maplet:::mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("Composite score analysis, %s", in_col),
      output = list(
        table = cstest,
        #formula = as.formula(glue::glue("~ {in_col}")),
        name    = stat_name,
        groups = outgroups,
        #lstobj = NULL,
        samples.used = samples.used,
        outcome = in_col
      )
    )

  ## return
  D

}
