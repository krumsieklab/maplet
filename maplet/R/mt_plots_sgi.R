# KC TO-DO
# - next thing to implement is 'tree pruning' (if needed)
# - perform analysis as long as one provided outcome is present in colData (with warning)?

#' SGI analysis plot
#'
#' Wrapper function for the \code{plot_overview} function from the \code{sgi} package. Performs a
#' simple SGI analysis and plots all parts in one combined overview plot. Refer to
#' \code{\link[sgi]{plot_overview}} for details.
#'
#' @param D A \code{SummarizedExperiment} object.
#' @param outcomes A vector of column names from colData. The outcomes to be tested by SGI.
#'    See \code{\link[sgi]{sgi_init}} for details.
#' @param distance_function The distance method to be passed to the \code{\link[stats]{dist}}
#'    function. Refer to link for details. Default: "euclidean".
#' @param hclust_method The clustering method to be passed to the \code{\link[stats]{hclust}}
#'    function. Refer to link for details. Default: "ward.D2".
#' @param min_size Minimum size for a cluster to be considered. If no value provided, set to
#'    5\% of sample size. See \code{\link[sgi]{sgi_init}} for details.
#' @param padj_thresh Adjusted p-value threshold to filtering out significant enrichments.
#'
#' @return $result$output: list containing ggplot object
#'
#' @examples
#' \donttest{
#'  D %>%
#'  mt_plots_sgi(outcomes = c("AGE", "SEX", "BMI")
#'               distance_function = "euclidean",
#'               hclust_method = "ward.D2",
#'               min_size = 18,
#'               padj_thresh = 0.05
#'              ) %>% ...
#'  }
#'
#' @author KC
#'
#' @export
mt_plots_sgi <- function(D,
                         outcomes,
                         distance_function = "euclidean",
                         hclust_method = "ward.D2",
                         min_size,
                         padj_thresh = 0.05){

  require(sgi)

  ### argument check ------
  if("SummarizedExperiment" %in% class(D)==F) stop("D must be a SummarizedExperiment object.")
  if(missing(outcomes)) stop("A value must be provided for argument \'outcomes\'.")
  if(padj_thresh > 1 || padj_thresh < 0) stop("\'padj_thresh\' must be a value between 0 and 1")


  ### prepare data and parameters ------
  # sgi expects samples in rows and features in columns
  data = t(assay(D))
  cd = as.data.frame(colData(D))

  # check if data contrains NAs
  if (any(is.na(data))) stop("Data matrix for sgi overview plot cannot contain NAs.")

  # check if any metabolites have zero variance
  if(any(data %>% apply(2, sd)==0)) stop("Some features have zero variance!")

  # if min cluster size not provided, use 5% sample size
  if(missing(min_size)) min_size = ceiling(nrow(data) * 0.05)

  # extract data frame of sample outcomes
  missing_cols <- outcomes[outcomes %in% colnames(cd) == F]
  if(length(missing_cols) > 0) stop(paste0("The following outcomes were not found in colData: ", paste0(missing_cols, collapse = ", ")))
  outcomes_df = cd %>% dplyr::select(all_of(outcomes))


  ### perform sgi analysis ------
  # hierarchical clustering
  hc = hclust(dist(data, method = distance_function), method = hclust_method)
  # initialize SGI structure
  sg = sgi_init(hc, minsize = min_size, outcomes = outcomes_df)
  # run sgi
  as = sgi_run(sg)

  ### generate tree and overview plot ------
  # subgroup tree plot
  gg_tree = plot(as, padj_th = padj_thresh)
  # overview plot - tree plot with sample outcome and feature heat maps
  p = plot_overview( gg_tree = gg_tree,
                     as = as,
                     outcomes = outcomes_df,
                     xdata = data )

  ### add status information & plot ------
  funargs <- mti_funargs()
  D %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = paste0("SGI overview plot with outcomes: ", paste0(outcomes, collapse=", ")),
      output = list(p)
    )

  D
}
