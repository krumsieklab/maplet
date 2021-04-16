#' Creates a volcano plot
#'
#' Create a volcano plot using data from a statistical result.
#'
#' @param D \code{SummarizedExperiment} input.
#' @param x Name of column in statistics table containing values to be plotted on the x-axis. Default: fc.
#' @param stat_name Name of the statistics object to plot.
#' @param feat_filter If given, filter will be applied to data and remaining variables will be labeled in plot
#' @param xlabel Label for the x-axis. Default: gsub("~","",as.character(x)).
#' @param vline Where to draw vertical line (for fold-change), has to be single value. Default: NA.
#' @param hline Where to draw horizontal line (for p-values), has to be an expression such as 'p.adj < 0.1'.
#' @param ggadd Further elements/functions to add (+) to the ggplot object. Default: NULL.
#' @param ... Additional expression directly passed to aes() of ggplot, can refer to colData.
#'
#' @return $result$output: plot, volcano
#'
#' @examples
#' \dontrun{# Volcano plot as overview of results with a result already in 'comp'
#' ... %>%
#' mt_plots_volcano(stat_name     = "comp",
#'  feat_filter = p.adj < 0.1,
#'  color       = p.value < 0.05) %>%
#'  ...}
#'
#' @author JZ, JK
#'
#' @import ggplot2
#'
#' @export
mt_plots_volcano <- function(D,
                             x = fc,
                             stat_name,
                             feat_filter = p.value < 0.05,
                             xlabel=gsub("~","",as.character(x)),
                             vline=NA,
                             hline,
                             ggadd=NULL,
                             ...){
  x <- dplyr::enquo(x)

  ## check input
  stopifnot("SummarizedExperiment" %in% class(D))
  if(missing(stat_name))
    stop("stat_name must be given for volcano plot")

  # get argument names from dots
  n <- sapply(as.list(substitute(list(...)))[-1L], deparse)
  dot_args <- names(n)
  
  # check for defunct argument names
  if ("statname" %in% dot_args) stop("You used the old MT naming convention statname Should be: stat_name")
  if ("metab_filter" %in% dot_args) stop("You used the old MT naming convention metab_filter Should be: feat_filter")
  
  ## rowData
  rd <- rowData(D) %>%
    as.data.frame() %>%
    dplyr::mutate(var = rownames(D))
  # remove rows not needed for plotting
  vars <- c(maplet:::mti_extract_variables(c(dplyr::enquo(x), dplyr::enquo(feat_filter), quos(...))),"var","name")
  rd <- rd[,colnames(rd) %in% vars,drop=F]


  ## stat
  data_plot <- maplet:::mtm_get_stat_by_name(D, stat_name)
  if(quo_name(x) %in% colnames(data_plot)==F) stop(glue::glue("Column {quo_name(x)} not found in stat table."))
  data_plot %<>% dplyr::inner_join(rd, by = "var") %>%
    dplyr::mutate(xxx = !!x)

  ## SCALE -log10
  reverselog_trans <- function (base = exp(1)){
    trans <- function(x) -log(x, base)
    inv <- function(x) base^(-x)
    scales::trans_new(paste0("reverselog-", format(base)), trans, inv,
                      scales::log_breaks(base = base),
                      domain = c(1e-100, Inf))
  }

  ## determine if and where to draw hline
  if (!missing(hline)) {
    hliney <- maplet:::mtm_get_stat_by_name(D, stat_name) %>%
      dplyr::inner_join(rd, by = "var") %>%
      dplyr::mutate(xxx = !!x) %>% dplyr::filter(!!dplyr::enquo(hline)) %>% .$p.value %>% max()
  } else {
    hliney <- NA
  }

  ## sanity check that there is something to plot
  if (all(is.na(data_plot$p.value))) stop("All p-values for Volcano plot are NA")

  ## CREATE PLOT
  p <- data_plot %>%
    ## do plot
    ggplot(aes(x = xxx, y = p.value)) +
    # vline?
    (if(!is.na(vline)){geom_vline(xintercept = c(-vline, vline), linetype='dashed', color='#F8766D')}else{NULL}) +
    # hline?
    (if(!is.na(hliney)){geom_hline(yintercept = hliney, linetype='dashed', color='#F8766D')}else{NULL}) +
    # points
    geom_point(aes(...)) +
    scale_y_continuous(trans = reverselog_trans(10),
                       breaks = scales::trans_breaks("log10", function(x) 10^x),
                       labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    labs(x = xlabel, y = "p-value") +
    ggtitle(stat_name)

  ## ADD FEATURE LABELS
  if(!missing(feat_filter)){
    maplet:::mti_logstatus("add label")
    feat_filter_q <- dplyr::enquo(feat_filter)
    data_annotate <- data_plot %>%
      dplyr::filter(!!feat_filter_q)
    p <- p + ggrepel::geom_text_repel(data = data_annotate,
                                      aes(label = name), max.overlaps = Inf)
  }

  ## ADD AXIS GROUPS
  d <- maplet:::mtm_get_stat_by_name(D, stat_name, fullstruct=T)
  if ("groups" %in% names(d) && length(d$groups)==2) {
    p <- maplet:::mti_add_leftright_gg(p, paste0(d$groups[1],' high'), paste0(d$groups[2],' high'))
  }

  # add custom elements?
  if (!is.null(ggadd)) p <- p+ggadd

  ## add status information & plot
  funargs <- maplet:::mti_funargs()
  D %<>% 
    maplet:::mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("volcano plot, aes: %s", maplet:::mti_dots_to_str(...)),
      output = list(p)
    )
  ## return
  D
}

