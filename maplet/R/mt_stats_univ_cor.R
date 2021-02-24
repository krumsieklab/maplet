#' Computes correlation to a given phenotype
#'
#' Test for association between paired samples. If present, NAs will be omitted.
#'
#' @param D \code{SummarizedExperiment} input.
#' @param stat_name Name under which this comparison will be stored, must be unique to all other statistical results.
#' @param in_col Name of the colData variable to use for the correlation calculation. If method="kendall", class(D[[in_col]]) needs
#'    to be numeric.
#' @param method Correlation method to use. Can be any among "pearson", "kendall", "spearman".
#' @param samp_filter OPTIONAL. Sample filter condition.
#' @param exact OPTIONAL. Set the exact flag in cor.test function. Default: NULL.
#'
#' @return $results$output: List of Kendall's correlation coefficients and pvalues, as well as the corresponding variable names.
#'
#' @examples
#' \dontrun{... %>%
#'   mt_stats_univ_cor(in_col = "Stage", samp_filter = (GROUP %in% "Tumor"), stat_name = "tau", method = "tau") %>%
#' ...
#' }
#'
#' @author EB, RB
#'
#' @export
mt_stats_univ_cor <- function(D, stat_name, in_col, method, samp_filter, exact=NULL) {

  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  # check method
  stopifnot(method %in% c("pearson", "kendall", "spearman"))
  # check that "in_col" is in the colData
  if (!(in_col %in% colnames(colData(D)))) stop(sprintf("There is no column called %s in the colData", in_col))
  # "in_col" must be numeric
  if (class(D[[in_col]])!="numeric") stop(sprintf("For Kendall's correlation, %s must be numeric",in_col))

  # make sure name does not exist yet
  if (stat_name %in% unlist(maplet::mtm_res_get_entries(D, "stats") %>% purrr::map("output") %>% purrr::map("stat_name"))) stop(sprintf("stat element with stat_name '%s' already exists",stat_name))

  # merge data with sample info
  Ds <- D %>% maplet:::mti_format_se_samplewise() # NOTE: No explosion of dataset size, no gather() - 6/2/20, JK

  ## FILTER SAMPLES
  if(!missing(samp_filter)) {

    filter_q <- dplyr::enquo(samp_filter)
    num_samp <- ncol(D)
    samples.used <- maplet:::mti_filter_samples(Ds, filter_q, num_samp)
    Ds <- Ds[samples.used,]

  } else {
    samples.used = rep(T, ncol(D))
  }

  met <- colnames(Ds)[(length(colnames(colData(D)))+2):length(colnames(Ds))]
  # compute association to the phenotype
  rr <- lapply(met, function(x){
    d=stats::cor.test(Ds[,x], Ds[[in_col]], method=method, alternative = "two.sided", exact=exact)
    list("statistic"=d$estimate, "p.value"=d$p.value, "method"=d$method)
  })
  names(rr) <- met

  # revert list structure
  revert_list_str <- function(ls) {
    # get sub-elements in same order
    x <- lapply(ls, `[`, names(ls[[1]]))
    # stack and reslice
    apply(do.call(rbind, x), 2, as.list)
  }
  rr_reverse <- revert_list_str(rr)

  # arrange results in dataframe
  tab <- cbind.data.frame(as.data.frame(do.call(rbind, rr_reverse$statistic)),
                         as.data.frame(do.call(rbind, rr_reverse$p.value)),
                         as.data.frame(do.call(rbind, rr_reverse$method)) )
  colnames(tab) <- c("statistic","p.value","method")

  # add term column with ordinal variable
  tab$term <- rep(in_col, dim(tab)[1])
  # add column with names
  tab$var <- rownames(tab)

  # reorder columns
  cc <- c("var","term","statistic","p.value")
  tab <- tab[,match(cc,colnames(tab))]

  ## construct output groups variable
  outgroups <- unique(Ds[[in_col]])

  # add status information
  funargs <- maplet:::mti_funargs()
  metadata(D)$results %<>%
    maplet:::mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("%s correlation to %s", method, in_col),
      output = list(
        table = tab,
        name = stat_name,
        lstobj = NULL,
        groups = outgroups,
        samples.used = samples.used,
        outcome = in_col
      )
    )

  # return
  D

}
