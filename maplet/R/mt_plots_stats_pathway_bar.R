#' Creates a pathway comparison bar plot
#'
#' Create a pathway comparison bar plot using a list of statistical results.
#'
#' @param D \code{SummarizedExperiment} input.
#' @param stat_list List of names of the statistics objects to be used for filtering
#' @param feat_filter Filter will be applied to data and remaining variables will be used to create plot
#' @param group_col The rowData variable used to aggregate variables.
#' @param color_col OPTIONAL. A rowData variable used to color barplot. Default NULL.
#' @param y_scale Plot percentage or frequency of variables. Values c("fraction","count"). Default "fraction".
#' @param sort_by_y Sort pathways in plot according to y_scale. Default FALSE.
#' @param assoc_sign_col OPTIONAL. Parameter to discriminate between positive and negative associations. Needs to be the name
#'    of a column in the statistical results indicated by stat_list.
#' @param add_empty BOOLEAN. If TRUE adds also empty pathways to the barplot.
#' @param keep_unmapped BOOLEAN. If TRUE keep NULL values (specified as "Unmapped"). Default: FALSE.
#' @param outfile OPTIONAL. Excel filename to save data to
#' @param ggadd Further elements/functions to add (+) to the ggplot object.
#' @param ... Additional expression directly passed to aes() of ggplot, can refer to colData.
#'
#' @return $result$output: plot, barplot
#'
#' @examples
#' \dontrun{# Barplot as overview of results with a result already in 'comp'
#' ... %>%
#' mt_plots_stats_pathway_bar(stat_list     = "comp",
#'  feat_filter = p.adj < 0.05,
#'  group_col    = "SUB_PATHWAY",
#'  color_col      = "SUPER_PATHWAYdevtools::install(codes.makepath("MT/maplet"))",
#'  y_scale       = "count",
#'  sort_by_y         = TRUE) %>%
#'  ...}
#'
#' @author Elisa Benedetti
#'
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#' @import openxlsx
#' @importFrom tidyselect any_of
#'
#' @export

mt_plots_stats_pathway_bar <- function(D,
                               stat_list,
                               feat_filter = p.value < 1,
                               group_col = "SUB_PATHWAY",
                               color_col = NULL,
                               y_scale = "fraction",
                               sort_by_y = FALSE,
                               assoc_sign_col,
                               add_empty = FALSE,
                               keep_unmapped = FALSE,
                               outfile = NULL,
                               ggadd = NULL,
                               ...){

  # get argument names from dots
  n <- sapply(as.list(substitute(list(...)))[-1L], deparse)
  dot_args <- names(n)

  # check for defunct argument names
  if ("stat_name" %in% dot_args) stop("You used the old MT naming convention stat_name. Should be: stat_list.")
  if ("metab_filter" %in% dot_args) stop("You used the old MT naming convention metab_filter. Should be: feat_filter.")
  if ("aggregate" %in% dot_args) stop("You used the old MT naming convention aggregate Should be: group_col.")
  if ("colorby" %in% dot_args) stop("You used the old MT naming convention colorby Should be: color_col.")
  if ("sort" %in% dot_args) stop("You used the old MT naming convention sort. Should be: sort_by_y.")
  if ("yscale" %in% dot_args) stop("You used the old MT naming convention yscale. Should be: y_scale.")
  if ("assoc_sign" %in% dot_args) stop("You used the old MT naming convention assoc_sign. Should be: assoc_sign_col.")
  if ("keep.unmapped" %in% dot_args) stop("You used the old MT naming convention keep.unmapped. Should be: keep_unmapped.")
  if ("output.file" %in% dot_args) stop("You used the old MT naming convention output.file. Should be: outfile.")

  ## check input
  stopifnot("SummarizedExperiment" %in% class(D))
  if(missing(stat_list) & !missing(feat_filter))
    stop("stat_list must be given for feat_filter to work.")
  if(missing(stat_list) & !missing(assoc_sign_col))
    stop("stat_list must be given for assoc_sign_col to work.")
  if(!(group_col %in% colnames(rowData(D))))
    stop(sprintf("group_col column '%s' not found in rowData", group_col))
  if(!is.null(color_col))
    if(!(color_col %in% colnames(rowData(D))))
      stop(sprintf("color_col column '%s' not found in rowData", color_col))

  ## rowData
  rd <- rowData(D) %>%
    as.data.frame() %>%
    dplyr::mutate(var = rownames(D))
  # set the nulls to unknown
  if(keep_unmapped){
    rd[[group_col]][which(rd[[group_col]]=="NULL")] <- "Unmapped"
  } else{
    rd <- rd[which(rd[[group_col]]!="NULL"), ]
  }
  perc <- rd[[group_col]] %>%
    unlist %>% table(exclude = NULL) %>% as.data.frame()
  colnames(perc) <- c("name","count")

  flag_filter <- ifelse((!missing(feat_filter)), T,F)
  flag_sign <- ifelse((!missing(assoc_sign_col)), T,F)

  data <- lapply(stat_list %>% {names(.)=.;.}, function(ss){
    ## subselect variables
    if(flag_filter) {
      feat_filter_q <- dplyr::enquo(feat_filter)
      sel <- maplet::mtm_get_stat_by_name(D=D,name=ss) %>%
        dplyr::filter(!!feat_filter_q) %>%
        dplyr::filter(var %in% rd$var)
      rd <- rd %>%
        dplyr::filter(var %in% sel$var)

    }

    # if filtering gives an empty matrix, produce an empty df
    if(nrow(rd)==0) {
      data_plot <- data.frame(name=as.character(),
                              count = as.numeric())
      anno <- data.frame()

    } else {
      # if assoc_sign_col given, include in data
      if(flag_sign){
        if(!(assoc_sign_col %in% colnames(sel))) {
          stop(sprintf("Could not find column called %s in the statistical results called %s", assoc_sign_col, stat_list))
        } else {
          # reorder sel according to rd
          sel <- sel[match(sel$var,rd$var),] %>%
            dplyr::mutate(association=ifelse(sign(!!sym(assoc_sign_col))>0, "positive", "negative")) %>%
            dplyr::select(var,association)

          data_plot <- data.frame(name=rd[[group_col]] %>% unlist %>% as.vector,
                                  association=rep(sel$association, times= (rd[[group_col]] %>% sapply(length)))) %>%
            table(exclude = NULL) %>% as.data.frame()
          colnames(data_plot) <- c("name","association","count")

        }
      } else {
        # reorder sel according to rd
        sel <- sel[match(sel$var,rd$var),] %>%
          dplyr::select(var)

        data_plot <- rd[[group_col]] %>%
          unlist %>% table(exclude = NULL) %>% as.data.frame()
        colnames(data_plot) <- c("name","count")

      }

      if(add_empty){
        # check which group_col entries are not included
        agg <- rowData(D) %>% as.data.frame() %>% .[[group_col]] %>% unlist %>% unique
        agg_empty <- agg[which(!(agg %in% unique(as.character(data_plot$name))))]

        # create data frame
        empty <- data.frame(name = agg_empty, count = rep(0, times=length(agg_empty)))

        if("association" %in% colnames(data_plot)){
          empty$association <- "positive"
        }

        # add to data
        data_plot %<>% dplyr::full_join(empty, by=colnames(data_plot))
      }

      # add number of features in each pathway
      perc <- data_plot %>% dplyr::select(name) %>%
        dplyr::left_join(perc, by="name")
      data_plot <- data_plot %>%
        # add fraction variable
        dplyr::mutate(fraction= count/perc$count)

      # add color column if not given
      if(is.null(color_col)) {
        color_col <- paste(group_col,"color", collapse = "_")
        rd[[color_col]] <- "pathway"
      }

      # create dictionary between group_col and color_col variables
      dict <- rd %>% dplyr::select(!!sym(group_col),!!sym(color_col)) %>% tidyr::unnest_longer(col=group_col) %>% as.data.frame()
      dict <- dict[!duplicated(dict[[group_col]]),]

      # add color to data_plot
      data_plot <- data_plot %>%
        dplyr::left_join(dict, by=c("name"=group_col)) %>%
        dplyr::rename(color=sym(color_col))

      # create annotation data
      anno <- data.frame(name = rep(rd$name, times=sapply(rd[[group_col]], length) %>% as.vector()),
                         var = rep(rd$var, times=sapply(rd[[group_col]], length) %>% as.vector()),
                         pathway = unlist(rd[[group_col]]),
                         color = if(!is.null(color_col)){rep(rd[[color_col]], times=sapply(rd[[group_col]], length) %>% as.vector())}else{"pathway"}) %>%
        dplyr::left_join(maplet::mtm_get_stat_by_name(D=D,name=ss) , by="var") %>%
        dplyr::select(-var)

      # if pathway mapping exists in the metadata, use the names provided there
      x <- D %>% metadata
      if ("pathways" %in% names(x)){
        if (group_col %in% names(x$pathways)) {
          # add pathway names to dataframe
          data_plot %<>%
            dplyr::left_join(x$pathways[[group_col]][,c("ID","pathway_name")], by=c("name"="ID"))
          anno %<>%
            dplyr::left_join(x$pathways[[group_col]][,c("ID","pathway_name")], by=c("pathway"="ID")) %>%
            dplyr::select(name,pathway_name,pathway,color,everything()) %>%
            dplyr::rename(pathway_id=pathway, pathway=pathway_name)

          # set Unknown pathway names to Unknown
          if(length(which(is.na(data_plot$pathway_name)))>0){
            data_plot$pathway_name[which(is.na(data_plot$pathway_name))] <- "Unknown"
          }
          # substitute codes for names and remove extra column
          data_plot %<>%
            dplyr::mutate(name=pathway_name) %>%
            dplyr::select(-pathway_name)
        } else{
          mti_logwarning(sprintf("%s field not found in the metadata",group_col))
        }
      }
      # create labels for plotting
      data_plot %<>% dplyr::mutate(label=sprintf("%s [%d]", name, perc$count))

      # convert labels to factor to sort alphabetically
      data_plot$label <- as.factor(data_plot$label)

      # add comparison name to df
      data_plot$comp <- ss
      anno$comp <- ss

    }
    list(dt = data_plot, anno = anno)
  })

  # function to revert string structure
  revert_list_str_4 <- function(ls) {
    # get sub-elements in same order
    x <- lapply(ls, `[`, names(ls[[1]]))
    # stack and reslice
    apply(do.call(rbind, x), 2, as.list)
  }

  data <- revert_list_str_4(data)

  # if there is at least one result, produce plot, otherwise output empty plot
  if((sapply(data$dt, function(ss){dim(ss)[1]}) %>% sum()) >0) {
    # merge list into a single dataframe
    data_plot <- do.call(rbind, data$dt) %>% as.data.frame()
    # get common colnames in stat table
    tt <- sapply(data$anno, colnames) %>% unlist %>% table %>% as.data.frame
    colnames(tt)[1] <- "var"
    tt %<>%
      dplyr::filter(Freq == max(Freq)) %>%
      dplyr::pull(var) %>%
      as.character
    anno <- lapply(data$anno, function(x){
      x %>%
        dplyr::select(tidyselect::any_of(tt))
    }) %>% {do.call(rbind, .)} %>% as.data.frame()

    # optional sorting (only for single statistical results)
    if (sort_by_y){
      data_plot$label <- reorder(data_plot$label, -data_plot[[y_scale]])
    }
    # sort comp so that facets appear in the same order given by the user
    data_plot$comp <- factor(data_plot$comp,levels=stat_list)

    # convert count to numeric
    data_plot$count %<>% as.numeric

    ## CREATE PLOT
    p <- ggplot(data_plot, aes(label)) +
      (if("association" %in% colnames(data_plot)) {geom_bar(data = subset(data_plot, association == "positive"), aes(y = !!sym(y_scale), fill = color), stat = "identity", position = "dodge", color="black", size=0.4)}) +
      (if("association" %in% colnames(data_plot)) {geom_bar(data = subset(data_plot, association == "negative"), aes(y = -!!sym(y_scale), fill = color), stat = "identity", position = "dodge", color="black", size=0.4)} else{geom_bar(aes(x=label, y=!!sym(y_scale), fill=color), stat = "identity", color="black", size=0.4)}) +
      (if(y_scale=="fraction") {ggtitle(sprintf("Fraction of pathway affected, %s", gsub("~", "", rlang::expr_text(dplyr::enquo(feat_filter)))))}else{ggtitle(sprintf("Number of hits per pathway, %s", gsub("~", "", rlang::expr_text(enquo(feat_filter)))))}) +
      (if(y_scale=="count" & "association" %in% colnames(data_plot)) {expand_limits(y=c(-max(data_plot$count, na.rm = T)*1.7, max(data_plot$count, na.rm = T)*1.7))}) +
      (if(y_scale=="count" & !("association" %in% colnames(data_plot))) {expand_limits(y=c(0, max(data_plot$count, na.rm = T)*1.7))}) +
      (if(y_scale=="fraction" & "association" %in% colnames(data_plot)) {expand_limits(y=c(-1, 1))}) +
      geom_hline(yintercept = 0,colour = "black", size=0.4) +
      labs(x="",fill = color_col) +
      theme(plot.title = element_text(hjust = 0.4)) +
      scale_x_discrete(limits = rev(levels(data_plot$label)))

    # add phenotype labels to x axis
    if("association" %in% colnames(data_plot) & length(stat_list)==1){
      d <- maplet::mtm_get_stat_by_name(D, stat_list, fullstruct=T)
      if ("groups" %in% names(d) && length(d$groups)==2) {
        # get breaks
        ggbld <- ggplot2::ggplot_build(p)
        yticks = ggbld$layout$panel_params[[1]]$y$minor_breaks # using minor_breaks because sometimes breaks would not work
        # edit labels to include groups
        ytlabs = yticks
        ytlabs[1] <- sprintf("%s\n%s", yticks[1], sprintf("high in %s", d$groups[1]))
        ytlabs[length(ytlabs)] <- sprintf("%s\n%s", yticks[length(yticks)], sprintf("high in %s", d$groups[2]))
        # apply new labels
        p <- p +
          scale_y_continuous(breaks = yticks, labels = ytlabs)
      }
    }

    # flip axes and add annotations on bars
    p <- p +
      coord_flip() +
      (if(y_scale=="count" & !("association" %in% colnames(data_plot))) {geom_text(data=data_plot, aes(label, !!sym(y_scale), label= sprintf("%.2f%%", fraction*100)),
                                                                               position = position_dodge(width=0.9), hjust = -0.1, size=2.5)}) +
      (if(y_scale=="count" & "association" %in% colnames(data_plot)) {geom_text(data=dplyr::filter(data_plot, association=="positive"), aes(label, !!sym(y_scale), group= association,label= sprintf("%.2f%%", fraction*100)),
                position = position_dodge(width=0.9), hjust = -0.1, size=2.5)}) +
      (if(y_scale=="count" & "association" %in% colnames(data_plot)) {geom_text(data=data_plot %>% dplyr::filter(association=="negative"), aes(label, -!!sym(y_scale), group= association,label= sprintf("%.2f%%", fraction*100)),
                position = position_dodge(width=0.9), hjust = 1.1, size=2.5)}) +
      facet_wrap(~comp)

    # add custom elements?
    if (!is.null(ggadd)) p <- p + ggadd

    # save plot parameters to be passed to the html generator for dynamical plot height
    re <- p %>%
      ggplot2::ggplot_build() %>%
      magrittr::extract2('layout') %>%
      magrittr::extract2('layout')

    nr <- data_plot$name %>% unique %>% length # number of pathways
    ncol <- re$COL %>% max() # number of panel columns
    nrow <- re$ROW %>% max() # number of panel rows

  } else {

    p <- ggplot() +
      geom_text(aes(x=0,y=0, label="No significant results"), size=10)

    # save plot parameters to be passed to the html generator for dynamical plot height
    nr <- 0 # number of pathways
    ncol <- NULL
    nrow <- NULL

  }

  if(!is.null(outfile)){
    if(exists("data_plot")){
      wb = openxlsx::createWorkbook()
      sheet = openxlsx::addWorksheet(wb, "Parameters")
      if(is.null(color_col)){color_col <- 'none'}
      openxlsx::writeData(wb, sheet=sheet, list(comparisons = stat_list, feat_filter = gsub("~", "", rlang::expr_text(enquo(feat_filter))), group_col = group_col, coloredby = color_col))
      sheet = openxlsx::addWorksheet(wb, "AggregatedPathways")
      openxlsx::writeData(wb, sheet=sheet, data_plot, rowNames = F, colNames = T)
      sheet = openxlsx::addWorksheet(wb, "IndividualResults")
      openxlsx::writeData(wb, sheet=sheet, anno, rowNames = F, colNames = T)
      openxlsx::saveWorkbook(wb, outfile, overwrite = T)
    } else {
      warning("mt_plots_statsbarplot: No significant results. outfile ignored.")
    }
  }

  ## add status information & plot
  funargs <- mti_funargs()
  D %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = ifelse(exists("stat_list"), sprintf("bar plot for comparison %s, by %s, filtered for %s, using %s", paste(stat_list,collapse = ", "), group_col, gsub("~", "", rlang::expr_text(enquo(feat_filter))), y_scale),
                      sprintf("bar plot by %s using %s", group_col, y_scale)),
      output = list(p),
      output2 = list(nr = nr, npancol = ncol, npanrow = nrow)
    )
  ## return
  D
}

