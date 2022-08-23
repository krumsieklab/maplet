# To Do List
# verify patients all have exactly one of each time point
# box plots or I plots
# time + group - coloring or grouping


#' Special Box Plots
#'
#' Create specialized box plots for paired samples (requires values for id_col and time_col) or time-resolved plots (requires values
#' for id_col, time_col, and group_col).
#'
#' @param D \code{SummarizedExperiment} input.
#' @param stat_name Index of the entry in metadata(D)$results that contains statistic object.
#' @param id_col Name of colData column containing sample IDs. Value for mapping sample pairs.
#' @param time_col Name of colData column containing factors or ordinal numbers representing time points.
#' @param group_col OPTIONAL. Name of colData column for grouping samples. Value on x-axis.
#' @param use_se For three or more time-points, use standard error (TRUE) or standard deviation (FALSE)? Default: TRUE.
#' @param feat_filter If given, filter will be applied to data and remaining variables will be labeled in plot. Default: p.value<0.05.
#' @param feat_sort If given, arrange will be applied to data variables will be sorted. Default: p.value.
#' @param annotation If given adds annotation to plot. Default: p.value (see Usage for formatting).
#' @param text_size Text size of the annotations. Default: 3.88.
#' @param restrict_to_used_samples Whether to filter to the samples that were used in the statistical test. Default: T.
#' @param ylabel Label for y-axis. Default: NULL.
#' @param ... Additional expression directly passed to aes() of ggplot, can refer to colData.
#'
#' @return $results$output: plot, special box plot
#'
#' @examples
#' \dontrun{mt_plots_repeating_test(stat_name = "comp2",
#' time_col           = Group,
#' group_col          = Sex,
#' id_col             = Individual,
#' fill               = Group,
#' feat_filter        = p.adj<0.2,
#' feat_sort          = fc)}
#'
#' @author KC
#'
#' @import ggplot2
mt_plots_box_special <- function(D,
                                 stat_name,
                                 id_col,
                                 time_col,
                                 group_col,
                                 use_se = TRUE,
                                 feat_filter = p.value < 0.05,
                                 feat_sort = p.value,
                                 annotation = "{sprintf('P-value: %.1e', p.value)}",
                                 text_size = 3.88,
                                 restrict_to_used_samples = T,
                                 ylabel=NULL,
                                 ...){

  if("SummarizedExperiment" %in% class(D)==F) stop("D is not a SummarizedExperiment object.")
  if(missing(id_col)) stop("The argument id_col is required.")
  if(missing(time_col)) stop("The argument time_col is required.")
  if(missing(stat_name)) stop("The argument stat_name is required.")

  # enquo variables
  time_col <- dplyr::enquo(time_col)
  id_col <- dplyr::enquo(id_col)

  n <- sapply(as.list(substitute(list(...)))[-1L], deparse)
  dot_args <- names(n)

  # create dummy SE so original not changed
  Ds <- D

  ## KC NOTE: Ask Jan if include correct confounders

  ## extract rowData dataframe
  rd <- rowData(Ds) %>%
    as.data.frame() %>%
    dplyr::mutate(var = rownames(Ds))

  ## extract stat table
  if(!missing(stat_name)){
    stat <- maplet::mtm_get_stat_by_name(Ds, stat_name) %>%
      dplyr::inner_join(rd, by = "var")
  }else{
    stat <- rd
    restrict_to_used_samples <- F # not dependend on a stat
  }

  ## FILTER FEATURES
  if(!missing(feat_filter)){
    feat_filter_q <- dplyr::enquo(feat_filter)
    stat <- stat %>%
      dplyr::filter(!!feat_filter_q)
    mti_logstatus(glue::glue("filter features: {feat_filter_q} [{nrow(stat)} remaining]"))
  }

  ## SORT FEATURES
  if(!missing(feat_sort)){
    feat_sort_q <- dplyr::enquo(feat_sort)
    stat <- stat %>%
      dplyr::arrange(!!feat_sort_q) %>%
      ## sort according to stat
      dplyr::mutate(name = factor(name, levels = unique(name)))
    mti_logstatus(glue::glue("sorted features: {feat_sort_q}"))
  }

  ## CREATE PLOT DATA FRAME
  dummy <- Ds %>%
    mti_format_se_samplewise() %>% # NOTE: No explosion of dataset size due to active restriction - 6/2/20, JK
    tidyr::gather(var, value, dplyr::one_of(rownames(Ds)))
  ## filter to groups?
  if (restrict_to_used_samples) {
    filterto <- maplet::mtm_get_stat_by_name(Ds, stat_name, fullstruct=T)$samples.used
    dummy <- dummy[filterto,]
  }

  # check x is a column in dataset
  mainvar <- time_col %>% dplyr::quo_name()
  if(mainvar %in% colnames(dummy) == F) stop(glue::glue("No column in plot data frame with name \"{mainvar}\"."))

  # filter down only to the variables needed for plotting
  # need to parse x and ... list
  vars <- time_col %>% dplyr::quo_name()
  vars <- c(vars, id_col %>% dplyr::quo_name())
  q <- dplyr::quos(...)
  if (length(q) > 0) {
    vars <- c(vars, q %>% lapply(function(x){x %>% as.character() %>% gsub("~","",.)}) %>% unlist() %>% as.vector())
  }
  vars <- unique(vars)
  mti_logstatus(glue::glue("these are the vars: {vars}"))

  plottitle <- ifelse(missing(stat_name),"",stat_name)
  # make sure the main outcome variable x is a factor
  mainvar <-time_col %>% dplyr::quo_name()
  dummy[[mainvar]] <- as.factor(dummy[[mainvar]])

  p <- dummy %>%
    dplyr::select(dplyr::one_of(c("var","value", vars))) %>%
    ## add feature names, but only restricted subset from statistics table
    dplyr::inner_join(stat[,dplyr::intersect(colnames(stat),c('var','statistic','p.value','p.adj','name'))], by = "var") %>%
    ggplot() +
    geom_boxplot(aes(x=as.factor(!!time_col), y=value, ...)) +
    geom_point(aes(x=as.factor(!!time_col), y=value, ...),colour="black", size=2, alpha=0.5) +
    geom_line(aes(x=as.factor(!!time_col), y=value, group=!!id_col, ...), colour="black")

  p <- p + facet_wrap(.~name, scales = "free_y", ncol=2)
  p <- list(p)
  output2 <- ceiling(length(unique(stat$name))/2)
  print(output2)

  # add status information & plot
  funargs <- mti_funargs()
  D %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("Feature ",ifelse("boxplots, aes: %s", mti_dots_to_str(...))),
      output = p,
      output2 = output2
    )

  # return
  D

}
