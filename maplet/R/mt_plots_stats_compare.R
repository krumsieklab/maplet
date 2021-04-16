#' Comparative plot between two statistical results
#'
#' Produces a plot that compares the directed -log10 p-values between two previously executed statistics.
#'
#' @param D1 First \code{SummarizedExperiment} input to compare; the one in the pipeline.
#' @param stat1 Name of statistical result in the first dataset.
#' @param filter1 Filter term, defining which features to label from first statistical result (can use elements of stats table).
#' @param D2 Second \code{SummarizedExperiment} input. If not given, will be the same as the first. Default: D1.
#' @param stat2 Name of statistical result in second dataset.
#' @param filter2 Filter term, defining which features to label from second statistical result (can use elements of stats table).
#' @param filter_op  If AND -> two colors, one for those where both stats match the criterion, and one where they don't.
#'                  If OR -> three colors, a third one where only one stat matches the criterion. Default: "AND".
#' @param title OPTIONAL. Title of plot. Default: "".
#' @param label_col OPTIONAL. Name of column in the statistical results data frame to use for labeling points. Default: "name".
#' @param point_size Size of the points on the ggplot. Default: 1.5.
#' @param return_plot_only Return only the plot object. WARNING: setting this to true makes the function non-MT
#'    pipeline compatible. Default: F.
#' @param data_outfile Name of Excel file to output statistical result data to. Default: NULL.
#' @param use_estimate Use estimate for comparison, instead of statistic. Default: F.
#'
#' @return $result$output: plot, p-value histogram
#' @return $result$output2: statistical result data frame
#'
#' @examples
#' \dontrun{## compare two stats from inside the same pipeline
#' ... %>%
#' mt_plots_stats_compare(stat1='WT',
#'   filter1= p.adj<0.1,
#'   stat2='KO',
#'   filter2= p.adj<0.1,
#'   filter_op = 'OR'
#' ) %>% ...
#'
#' ## compare two stats from different pipelines, as part of the pipeline of the second
#' # 'comp' is a string that contains the name of a comparison (here both SEs have the same comparison on two datasets)
#' .. %>% mt_plots_stats_compare(
#'   stat1 = comp, filter1 = p.adj<0.1,
#'   D2 = firstPipeSE, stat2 = comp, filter2 = p.adj<0.1,
#'   filter_op = "OR") %>% ...
#'
#' ## compare two stats from different pipelines, output as plot object
#' ## not part of the actual MT pipelines, but separate call
#' # 'comp' is a string that contains the name of a comparison (here both SEs have the same comparison on two datasets)
#' gg <- mt_plots_stats_compare(
#'   D1 = D1, stat1 = comp, filter1 = p.adj<0.1,
#'   D2 = D2, stat2 = comp, filter2 = p.adj<0.1,
#'   filter_op = "OR", return_plot_only=T)
#' }
#'
#' @author JK
#'
#' @import ggplot2
#'
#' @export
mt_plots_stats_compare <- function(D1,
                                   stat1,
                                   filter1,
                                   D2 = D1,
                                   stat2,
                                   filter2,
                                   filter_op="AND",
                                   title = "",
                                   label_col = "name",
                                   point_size = 1.5,
                                   return_plot_only=F,
                                   data_outfile = NULL,
                                   use_estimate = F) {

  ## check input
  stopifnot("SummarizedExperiment" %in% class(D1))
  stopifnot("SummarizedExperiment" %in% class(D2))
  if (missing(filter1)) stop("Must provide 'filter1'")
  if (missing(filter2)) stop("Must provide 'filter2'")

  filter1q <- dplyr::enquo(filter1)
  filter2q <- dplyr::enquo(filter2)

  if (!(filter_op %in% c("AND","OR"))) stop("filter_op must be 'AND' or 'OR'")

  ## obtain the two stats structures
  s1 <- maplet:::mtm_get_stat_by_name(D1, stat1, fullstruct=T)
  s1t <- s1$table
  s2 <- maplet:::mtm_get_stat_by_name(D2, stat2, fullstruct=T)
  s2t <- s2$table

  # if use_estimate==T, check that estimate column exists
  if(use_estimate){
    if("estimate" %in% colnames(s1t) == F){stop("Column estimate was not found in stat table 1")}
    if("estimate" %in% colnames(s2t) == F){stop("Column estimate was not found in stat table 2")}
  }

  ## add directed p-value, filter, and merge
  s1t$dp1 <- if(use_estimate) -log10(s1t$p.value) * sign(s1t$estimate) else -log10(s1t$p.value) * sign(s1t$statistic)
  s1t$filtered1 <- s1t$var %in% (s1t %>% dplyr::filter(!!filter1q))$var
  s2t$dp2 <- if(use_estimate) -log10(s2t$p.value) * sign(s2t$estimate) else -log10(s2t$p.value) * sign(s2t$statistic)
  s2t$filtered2 <- s2t$var %in% (s2t %>% dplyr::filter(!!filter2q))$var
  st <- merge(s1t, s2t, by='var')
  st <- merge(st, rowData(D1), by.x="var", by.y='row.names', all.x=T) # add names


  # combine filters
  if (filter_op=="AND")
    # AND
    st$filtered = as.numeric(st$filtered1 & st$filtered2)
  else if (filter_op=="OR")
    # OR
    st$filtered = as.numeric(st$filtered1) + as.numeric(st$filtered2)
  else
    stop("bug")

  ## create axis labels
  if ("groups" %in% names(s1) && length(s1$groups)==2) {
    xlabel = sprintf("%s high <--   dir. log10(p-value)   --> %s high", s1$groups[1], s1$groups[2])
  } else {
    xlabel = 'directed log10(p-value)'
  }
  if ("groups" %in% names(s2) && length(s2$groups)==2) {
    ylabel = sprintf("%s high <--   dir. log10(p-value)   --> %s high", s2$groups[1], s2$groups[2])
  } else {
    ylabel = 'directed log10(p-value)'
  }


  ## plot
  st <- as.data.frame(st)
  p <- st %>%
    ggplot2::ggplot(ggplot2::aes(x=dp1,y=dp2,color=as.factor(filtered))) +
    ggplot2::geom_point(size = point_size) +
    ggplot2::labs(color='filtered') +
    ggrepel::geom_text_repel(data=dplyr::filter(st, filtered>0), ggplot2::aes_string(label=label_col), size=3, colour = "black",
                             max.overlaps = Inf) +
    ggplot2::xlab(xlabel) + ggplot2::ylab(ylabel)

  if (title != "") {
    p <- p + ggplot2::ggtitle(title)

  }

  ## export to file?
  if (!is.null(data_outfile)) {
    # can't handle list columns, drop those
    keep = !sapply(st, is.list)
    openxlsx::write.xlsx(x=st[,keep], file=data_outfile, asTable=F)
  }

  if (!return_plot_only) {
    ## add status information & plot
    funargs <- maplet:::mti_funargs()
    D1 %<>%
      maplet::: mti_generate_result(
        funargs = funargs,
        logtxt = sprintf("comparison plot between '%s' and '%s'", stat1, stat2),
        output = list(p),
        output2 = st
      )
    ## return
    D1
  } else {
    p
  }

}
