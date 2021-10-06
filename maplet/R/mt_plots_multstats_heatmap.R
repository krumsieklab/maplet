#' Heatmap of Multiple Statistical Objects
#'
#' Generate a heatamp for a list of statistical results using \code{pheatmap::pheatmap}.
#'
#' @param D \code{SummarizedExperiment} input.
#' @param stat_list List of stat names to plot. Default: NA (for all stat objects in D).
#' @param cutoff The p.adj value significance cutoff.
#' @param feat_anno_col The rowData() column name to use as column annotation.
#' @param signif_mark Marker used to indicate significant p-values. Default: "•".
#' @param mark_size Size of the marker. Default: 20.
#' @param show_mark Show significance marker? Default: T.
#' @param color_signif Boolean. Color only significant values? Default: F.
#' @param filter_signif Show only features significant in one or more results? Default: F.
#' @param cluster_rows Cluster by rows? Default: F.
#' @param cluster_cols Cluster by cols? Default: F.
#' @param show_colnames Should column (feature) labels be included in the plot? Default: F.
#' @param main Main title of pheatmap. Default: "Stat Results Overview".
#'
#' @return results$output: heatmap of statistical results
#'
#' @examples
#' \dontrun{D <- D %>% mt_plots_multstats_heatmap(feat_anno_col = "SUPER_PATHWAY",
#'                                                cluster_cols = T,
#'                                                cluster_rows = T)}
#'
#' @author KC
#'
#' @import ggplot2
#'
#' @export
mt_plots_multstats_heatmap <- function(D,
                                       stat_list=NA,
                                       cutoff,
                                       feat_anno_col,
                                       signif_mark = "•",
                                       mark_size = 20,
                                       show_mark = T,
                                       color_signif = F,
                                       filter_signif = F,
                                       cluster_rows = F,
                                       cluster_cols = F,
                                       show_colnames = F,
                                       main = "Stat Results Overview"){

  pheat_arg <- list()
  pheat_arg$cluster_rows <- cluster_rows
  pheat_arg$cluster_cols <- cluster_cols
  pheat_arg$show_colnames <- show_colnames
  pheat_arg$main <- main

  ## NTS: add checks for all parameters
  # validate and extract arguments
  stopifnot("SummarizedExperiment" %in% class(D))

  if(color_signif==T){
    if(missing(cutoff)){
      stop("color_signif is set to True, but no values provided for cutoff")
    }else if(length(cutoff) != 1){
      stop("cutoff must consist of two and only two elements")
    }else if(!is.numeric(cutoff)){
      stop("cutoff must contain only numeric values")
    }
  }

  cutoff <- -log10(cutoff)
  left_bound <- -cutoff
  right_bound <- cutoff

  # if NA, stat_list is all stat results
  if(is.na(stat_list)){
    stat_list <- maplet::mtm_res_get_entries(D, "stats")%>%
      purrr::map_chr(~.x$output$name)
  }

  if(length(stat_list) == 0){
    stop("There are no statistical results.")
  }

  # extract stat results
  stat_res <- stat_list %>% lapply(function(i){mtm_get_stat_by_name(D, i) %>%
      dplyr::mutate(stat_name = i)
  })

  # get stat names
  stat_names <- stat_list %>% unname()

  # determine which features are significant if signif_expr provided

  # generate matrix to plot
  color_signif_matrix <- stat_res %>% lapply(function(x){
    stat_name = x$stat_name[1]
    x <- x %>%  dplyr::mutate(color = sign(statistic) * -log10(p.adj)) %>%
      dplyr::mutate(signif_col = ifelse(color > left_bound & color < right_bound , 0, 1)) %>%
      dplyr::select(., var, color, signif_col)
    colnames(x)[2] <- stat_name
    colnames(x)[3] <- paste0(stat_name, "_signif_col")
    x
  }) %>% purrr::reduce(dplyr::full_join, by="var")

  # filter significant features
  if(filter_signif == T){
    color_signif_matrix <- color_signif_matrix %>% dplyr::filter_at(dplyr::vars(dplyr::ends_with("_signif_col")), dplyr::any_vars(. == 1))
  }

  color_matrix <- color_signif_matrix %>% dplyr::select(., var, tidyselect::all_of(stat_names))
  signif_matrix <- color_signif_matrix %>% dplyr::select(., var, dplyr::ends_with("_signif_col"))

  # generate significance label table
  if(show_mark==T){
    signif_label <- signif_matrix %>% t()
    colnames(signif_label) <- signif_label[1,,drop=F]
    signif_label <- signif_label[-1,,drop=F]
    signif_label <- ifelse(signif_label==1,signif_mark,"") %>% as.data.frame()
    indx <- sapply(signif_label, is.factor)
    signif_label[indx] <- lapply(signif_label[indx], function(x) as.character(x))
    rownames(signif_label) <- gsub("_signif_col", "", rownames(signif_label))
    pheat_arg$display_numbers <- signif_label
    pheat_arg$fontsize_number <- mark_size
  }


  # add annotations?
  if(!missing(feat_anno_col)){
    rd <- rowData(D) %>% as.data.frame()
    tmp <- rd[,c("name", feat_anno_col)] %>% as.data.frame()
    colnames(tmp)[1] <- "var"
    tmp$var <- make.names(tmp$var)
    tmp <- tmp %>% dplyr::semi_join(color_matrix, by="var")
    rownames(tmp) <- tmp$var
    tmp <- tmp %>% dplyr::select(all_of(feat_anno_col))

    pheat_arg$annotation_col <- tmp
  }

  # format data for plotting
  color_matrix <- color_matrix[,c("var", stat_names),drop=F] %>% t()
  colnames(color_matrix) <- color_matrix[1,]
  color_matrix <- color_matrix[-1,,drop=F] %>% as.data.frame()
  indx <- sapply(color_matrix, function(x){is.factor(x) | is.character(x)})
  color_matrix[indx] <- lapply(color_matrix[indx], function(x) as.numeric(as.character(x)))
  # if there's just one row, make sure row clustering is deactivated
  if (nrow(color_matrix)==1) pheat_arg$cluster_rows=F
  # set plot argument
  pheat_arg$mat <- color_matrix



  # maximum absolut value
  range <- max(abs(color_matrix), na.rm=T)
  # check if anything is significant at all
  if (range>cutoff) {

    # if coloring by significance, need to calculate breaks
    if(color_signif == T){

      # generate colors and color breaks
      colors <- gplots::bluered(100)
      n = length(colors)

      color_breaks <- seq(-range, range, length.out = n + 1)

      # find breaks that represent significant values
      graycolor <- "#CCCCCC"
      left_side <- {color_breaks < -cutoff} %>% which() %>% max() # values between -cutoff and cutoff will have the same constant color
      right_side <- {color_breaks > cutoff} %>% which() %>% min()


      color_breaks <- c(color_breaks[1:left_side], -cutoff, cutoff, color_breaks[right_side:length(color_breaks)])
      # rearrange colors
      colors <- c(colors[1:(left_side)], graycolor, colors[(right_side-3):length(colors)])

      pheat_arg$color <- colors
      pheat_arg$breaks <- color_breaks


    }

    # plot pheatmap
    re <- do.call(pheatmap::pheatmap, pheat_arg)

  } else {
    # empty plot
    re <- (ggplot() + ggtitle(main)+ geom_text(aes(x=0, y=0, label='no significant results'), size=8))
  }

  # NTS: add option to convert to ggplot

  # add status information & plot
  funargs <- mti_funargs()
  D %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("Pheatmap of stat results.\n Colors represent (sign(statistic) * -log10(p.adj)))"),
      output = list(re)
    )

  # return
  D


}
