#' Plot Box or Scatter Plots
#'
#' Creates one box plot or scatter plot per feature based on given sample annotations.
#'
#' @param D \code{SummarizedExperiment} input.
#' @param x Name of column from colData containing the phenotype to plot on x axis.
#' @param plot_type Either "box" or "scatter".
#' @param stat_name Index of the entry in metadata(D)$results that contains statistic object.
#' @param correct_confounder Confounders to adjust for before plotting, formula notation.
#' @param feat_filter If given, filter will be applied to data and remaining variables will be labeled in plot. Default: p.value<0.05.
#' @param feat_sort If given, arrange will be applied to data variables will be sorted. Default: p.value.
#' @param annotation If given adds annotation to plot. Default: p.value (see Usage for formatting).
#' @param full_info Add full information of all sample annotations and statistics results to plot table data.frame?
#'    Makes plotting more flexible but can render SE objects huge. Default: F
#' @param text_size Text size of the annotations. Default: 3.88.
#' @param jitter Add geom_jitter ("jitter") or geom_beeswarm ("beeswarm") to boxplot. Exclude if NULL.  Default: "beeswarm".
#' @param restrict_to_used_samples Whether to filter to the samples that were used in the statistical test. Default: T.
#' @param ylabel Label for y-axis. Default: NULL.
#' @param fit_line Add fit line? Default: T.
#' @param fit_line_se Add standard error range? Default: T.
#' @param ggadd Further elements/functions to add (+) to the ggplot object. Default: NULL.
#' @param pages Whether to divide plots into pages (repeats legend and y-axis label). Default F.
#' @param ... Additional expression directly passed to aes() of ggplot, can refer to colData.
#'
#' @return if plot_type = "box", $result$output: plot, box plot
#' @return if plot_type = "scatter", $result$output: plot, scatter plot
#'
#' @examples
#' \dontrun{# boxplots as overview of results with a result already in 'comp'
#' # color by "Group" variable in colData
#' mt_plots_box_scatter(x                  = Group,
#'                      stat_name          = "comp",
#'                      plot_type          = "box",
#'                      correct_confounder = ~BATCH_MOCK,
#'                      feat_filter       = p.value<0.01,
#'                      feat_sort         = p.value,
#'                      annotation         = "{sprintf('P-value: %.1e', p.value)}\nStatistic: {sprintf('%.2f', statistic)}",
#'                     ) %>%
#'                  ...
#'  }
#'
#' @author JZ, KC
#'
#' @import ggplot2
#'
#' @export
mt_plots_box_scatter <- function(D,
                                 x,
                                 stat_name,
                                 plot_type,
                                 correct_confounder,
                                 feat_filter = p.value < 0.05,
                                 feat_sort = p.value,
                                 annotation = "{sprintf('P-value: %.1e', p.value)}",
                                 full_info = F,
                                 text_size = 3.88,
                                 jitter = "beeswarm",
                                 restrict_to_used_samples = T,
                                 ylabel=NULL,
                                 fit_line = T,
                                 fit_line_se = T,
                                 ggadd = NULL,
                                 pages = F,
                                 ...){

  stopifnot("SummarizedExperiment" %in% class(D))
  if(missing(x)) stop("x must be provided")
  x <- dplyr::enquo(x)

  # get argument names from dots
  n <- sapply(as.list(substitute(list(...)))[-1L], deparse)
  dot_args <- names(n)

  # check for defunct argument names
  if ("metab_filter" %in% dot_args) stop("You used the old MT naming convention metab_filter. Should be: feat_filter")
  if ("metab_sort" %in% dot_args) stop("You used the old MT naming convention metab_sort. Should be: feat_sort")
  if ("rows" %in% dot_args) stop("\"rows\" is no longer an acceptecd argument.")
  if ("cols" %in% dot_args) stop("\"cols\" is no longer an accepted argument.")
  if ("manual_ylab" %in% dot_args) stop("You used the old MT naming convention manual_ylab. Should be: ylabel")
  if ("manual_ylabel" %in% dot_args) stop("You used the old MT naming convention manual_ylabel. Should be: ylabel")
  if ("fitline" %in% dot_args) stop("You used the old MT naming convention fitline. Should be: fit_line")
  if ("fitline_se" %in% dot_args) stop("You used the old MT naming convention fitline_se. Should be: fit_line_se")

  # create dummy SE so original not changed
  Ds <- D

  ## CONFOUNDER
  if(!missing(correct_confounder)){
    mti_logstatus(glue::glue("correcting for {correct_confounder}"))
    Ds <- mti_correctConfounder(Ds, correct_confounder)
  }

  ## rowData
  rd <- rowData(Ds) %>%
    as.data.frame() %>%
    dplyr::mutate(var = rownames(Ds))

  ## stat
  if(!missing(stat_name)){
    stat <- maplet::mtm_get_stat_by_name(Ds, stat_name) %>%
      dplyr::inner_join(rd, by = "var")
  }else{
    stat <- rd
    ### KC: ONLY IN BOXPLOT (Why?)
    restrict_to_used_samples <- F # not dependend on a stat
  }

  ## FILTER FEATURES
  ### KC: feat_filter is never not missing, should there be no default?
  if(!missing(feat_filter)){
    feat_filter_q <- dplyr::enquo(feat_filter)
    stat <- stat %>%
      dplyr::filter(!!feat_filter_q)
    mti_logstatus(glue::glue("filter features: {feat_filter_q} [{nrow(stat)} remaining]"))
  }

  ## SORT FEATURES
  ### KC: feat_sort is never not missing, should there be no default?
  if(!missing(feat_sort)){
    feat_sort_q <- dplyr::enquo(feat_sort)
    stat <- stat %>%
      dplyr::arrange(!!feat_sort_q) %>%
      ## sort according to stat
      dplyr::mutate(name = factor(name, levels = unique(name)))
    mti_logstatus(glue::glue("sorted features: {feat_sort_q}"))
  }

  ## CREATE PLOT
  dummy <- Ds %>%
    mti_format_se_samplewise() %>% # NOTE: No explosion of dataset size due to active restriction - 6/2/20, JK
    tidyr::gather(var, value, dplyr::one_of(rownames(Ds)))
  ## filter to groups?
  if(plot_type=="box"){
    if (restrict_to_used_samples) {
      filterto <- maplet::mtm_get_stat_by_name(Ds, stat_name, fullstruct=T)$samples.used
      dummy <- dummy[filterto,]
    }
  }

  # check x is a column in dataset
  mainvar <- x %>% dplyr::quo_name()
  if(mainvar %in% colnames(dummy) == F) stop(glue::glue("No column in plot data frame with name \"{mainvar}\"."))

  if(!full_info){
    # filter down only to the variables needed for plotting
    # need to parse x and ... list
    vars <- x %>% dplyr::quo_name()
    q <- dplyr::quos(...)
    if (length(q) > 0) {
      vars <- c(vars, q %>% lapply(function(x){x %>% as.character() %>% gsub("~","",.)}) %>% unlist() %>% as.vector())
    }
    vars <- unique(vars)

    plottitle <- ifelse(missing(stat_name),"",stat_name)
    if(plot_type=="box"){
      # make sure the main outcome variable x is a factor
      mainvar <-x %>% dplyr::quo_name()
      dummy[[mainvar]] <- as.factor(dummy[[mainvar]])

      p <- dummy %>%
        dplyr::select(dplyr::one_of(c("var","value", vars))) %>%
        ## add feature names, but only restricted subset from statistics table
        dplyr::inner_join(stat[,dplyr::intersect(colnames(stat),c('var','statistic','p.value','p.adj','name'))], by = "var") %>%
        dplyr::select(-var) %>%
        ## do plot
        ggplot() +
        geom_boxplot(aes(x = as.factor(!!x), y = value, ...), outlier.shape = ifelse(jitter, NA, 19)) +
        labs(x = NULL, y = NULL) +
        ggtitle(plottitle)
    }else{

      p <- dummy %>%
        dplyr::select(dplyr::one_of(c("var","value", vars))) %>%
        ## add feature names, but only restricted subset from statistics table
        dplyr::inner_join(stat[,c('var','statistic','p.value','p.adj','name')], by = "var") %>%
        dplyr::select(-var) %>%
        ## do plot
        ggplot() +
        ## add fit line?
        {if (fit_line) geom_smooth(aes(x = !!x, y = value), method = "lm", se=fit_line_se, color = "black") else NULL} +
        geom_point(aes(x = !!x, y = value, ...)) +
        labs(x = dplyr::quo_name(x), y="feature") +
        ggtitle(plottitle)

    }

  }else{

    if(plot_type=="box"){
      # leave full info in
      # can create huge data.frames

      plottitle <- ifelse(missing(stat_name),"",stat_name)
      p <- dummy %>%
        ## add feature names
        dplyr::inner_join(stat, by = "var") %>%
        ## do plot
        ggplot() +
        geom_boxplot(aes(x = !!x, y = value, ...), outlier.shape = ifelse(jitter, NA, 19)) +
        labs(x = NULL, y = NULL) +
        ggtitle(plottitle)

    }else{
      # leave full info in
      # can create huge data.frames

      #
      plottitle <- ifelse(missing(stat_name),"",stat_name)
      p <- dummy %>%
        ## add feature names
        dplyr::inner_join(stat, by = "var") %>%
        ## do plot
        ggplot() +
        ## add fit line?
        {if (fit_line) geom_smooth(aes(x = !!x, y = value), method = "lm", se=fit_line_se, color = "black") else NULL} +
        geom_point(aes(x = !!x, y = value, ...)) +
        labs(x = dplyr::quo_name(x), y="feature") +
        ggtitle(plottitle)

    }

  }

  ### BOX PLOT SPECIFIC
  if(plot_type=="box"){
    ## add ylabel
    if (!is.null(ylabel)) {
      p <- p + ylab(ylabel)
    } else {
      # add label if this is logged data
      r <- Ds %>% maplet::mtm_res_get_entries(c("pre","trans","log"))
      if (length(r)>0) {
        p <- p + ylab(r[[1]]$logtxt) # log text contains e.g. "log2"
      }
    }

    ## ADD JITTER
    if(jitter=="beeswarm"){
      p <- p +
        ggbeeswarm::geom_beeswarm(aes(x = !!x, y = value, ...))
    }else if(jitter=="jitter"){
      p <- p +
        geom_jitter(aes(x = !!x, y = value, ...))
    }
  }


  ### COMMON TO BOTH PLOTS
  ## ADD ANNOTATION
  if(!missing(annotation)){
    data_annotate <- stat %>%
      dplyr::mutate(annotate = glue::glue(annotation)) %>%
      dplyr::distinct(name, annotate)
    p <- p + geom_text(data = data_annotate,
                       aes(label = annotate),
                       x = -Inf, y = Inf, hjust = -0.05, vjust = 1.05, size=text_size )
  }

  if (!is.null(ggadd)) p <- p+ggadd

  ## SPLIT TO MULTIPLE PAGES
  # if there is no plot, create a single empty page
  if (length(unique(stat$name))==0) {
    p <- list(ggplot() + geom_text(aes(x=0, y=0, label='no plots'), size=10))
    output2 <- NULL
  } else if(pages) {

    npages <- ceiling(length(unique(stat$name))/8)
    tmp <- lapply(1:npages, function(x){
      p + ggforce::facet_wrap_paginate(.~name, scales = "free", ncol = 2, nrow = 4, page=x)
      })
    p <- ggpubr::ggarrange(plotlist = tmp, ncol=1)
    p <- ggpubr::annotate_figure(p, top = ggpubr::text_grob(plottitle))
    output2 <- ceiling(length(unique(stat$name))/2)
  }else {
    p <- p + facet_wrap(.~name, scales = "free", ncol=2)
    p <- list(p)
    output2 <- ceiling(length(unique(stat$name))/2)
  }

  ## add status information & plot
  funargs <- mti_funargs()
  D %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("Feature ",ifelse(plot_type=="box", "boxplots", "scatter plots"),", aes: %s", mti_dots_to_str(...)),
      output = p,
      output2 = output2
    )
  ## return
  D

}

