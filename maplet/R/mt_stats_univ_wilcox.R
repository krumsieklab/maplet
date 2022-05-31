#' Wilcox test
#'
#' Perform Wilcoxon test.
#'
#' @param D \code{SummarizedExperiment} input.
#' @param in_col Name of the column in colData to compare with feature.
#' @param stat_name Name under which this comparison will be stored, must be unique to all other statistical results.
#' @param samp_filter Term which samples to filter to first (e.g. used if the data contains >2 groups but the user
#'    wants to run a two-group comparison).
#' @param paired OPTIONAL. Set the paired flag in wilcox.test function, would be applied for numeric in_col. Default: F.
#' @param exact OPTIONAL. Set the exact flag in wilcox.test function. Default: NULL.
#'
#' @return $results$output: Statistics object.
#'
#' @examples
#' \donttest{# run lm with no confounders, "Group" as outcome
#' # filter to groups "Li_2" and "Li_5"
#' # name the comparison "Li's"
#' ... %>%
#'  mt_stats_univ_wilcox(
#'    in_col     = Group,
#'    samp_filter = (Group %in% c("Li_2","Li_5")),
#'    stat_name         = "Li's"
#'  ) %>% ...
#'  }
#'
#' @author RB
#'
#' @import glue
#'
#' @export
mt_stats_univ_wilcox <- function(D, in_col, stat_name, samp_filter, paired=F, exact=NULL) {

  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  # check that in_col is in the colData
  if (!(in_col %in% colnames(colData(D)))) stop(sprintf("There is no column called %s in the colData", in_col))
  # make sure name does not exist yet
  if (stat_name %in% unlist(maplet::mtm_res_get_entries(D, "stats") %>% purrr::map("output") %>% purrr::map("name"))) stop(sprintf("stat element with name '%s' already exists",stat_name))

  # merge data with sample info
  Ds <- D %>% mti_format_se_samplewise() # NOTE: No explosion of dataset size, no gather() - 8/17/20, JK

  ## FILTER SAMPLES
  if(!missing(samp_filter)) {

    filter_q <- dplyr::enquo(samp_filter)
    num_samp <- ncol(D)
    samples.used <- mti_filter_samples(Ds, filter_q, num_samp)
    Ds <- Ds[samples.used,]

  } else {
    samples.used = rep(T, ncol(D))
  }

  # check that outcome is either binary or numerical and run test accordingly
  mets <- rownames(D)
  outvec <- Ds[[in_col]]
  cl <- outvec %>% class()

  if (("character" %in% cl) || ("factor" %in% cl)) {
    if ((outvec %>% as.factor() %>% levels() %>% length()) != 2) {
      stop("If outcome is a factor, it must have exactly two levels")
    }
    # save groups
    outgroups <- outvec %>% as.factor() %>% levels()
    # run wilcox
    # we want to model feature ~ outcome
    wt <- lapply(mets, function(x){
      input <- split(Ds[,x], Ds[[in_col]])
      res  <- wilcox.test(input[[1]], input[[2]], alternative = "two.sided", exact=exact)
      res <- data.frame("statistic"=res$statistic, "p.value"=res$p.value, "method"=res$method)
      return(res)
    }) %>% do.call(rbind,.) %>% data.frame()

  } else {
    outgroups <- NULL
    # run wilcox
    # we want to model feature ~ outcome
    wt <- lapply(mets, function(x){
      res  <- wilcox.test(Ds[,x], Ds[[in_col]], paired=paired, alternative = "two.sided", exact=exact)
      res <- data.frame("statistic"=res$statistic, "p.value"=res$p.value, "method"=res$method)
      return(res)
    }) %>% do.call(rbind,.) %>% data.frame()
  }

  # add columns with feature names and in_col variable
  wt %<>% mutate(var=mets, term=in_col)
  # reorder columns
  wt %<>% select(var, term, statistic, p.value)
  # rearrange back to original feature order
  o <- match(rownames(D),wt$var)
  stopifnot(!any(is.na(o))) # sanity check
  wt <- wt[o,]
  # make sure that NAs in the outcome are set to FALSE in the list of used samples
  samples.used[is.na(Ds[[in_col]])] <- F

  ## add status information & results
  funargs <- mti_funargs()
  D %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("Wilcox rank sum test, %s", in_col),
      output = list(
        table = wt,
        #formula = as.formula(glue::glue("~ {in_col}")),
        formula = sprintf("~ %s", in_col),
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
