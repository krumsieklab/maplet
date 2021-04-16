#' 'Equalizer' plots
#'
#' Creates a nested plot based on feature annotations, e.g. of super-/ and sub-pathways, or of sub-pathways and features.
#'
#' @param D1 \code{SummarizedExperiment} input 1, the coarse one.
#' @param stat1  Name of first comparison output to take arguments from, the coarse one [first one has to be the less
#'    granular one (e.g. D1 super, D2 sub)].
#' @param D2 \code{SummarizedExperiment} input 2, the fine one.
#' @param stat2 Name of second comparison output to take arguments from, the fine one.
#' @param legend_fine Fine legend to be plotted.
#' @param legend_coarse Coarse legend to be plotted. Default: NULL.
#' @param vline_fine Filter expression where to draw the red, dashed line, for fine. Default: p.adj < 0.05.
#' @param vline_coarse Filter expression where to draw the red, dashed line, for coarse. Default: p.adj < 0.05.
#' @param color_list Colors for fine and coarse. Default: c("#9494FF","red") (light blue and red).
#'
#' @return $result$output: plot, equalizer
#'
#' @examples
#' \dontrun{# super-pathway / sub-pathway equalizer
#' # sub-pathway analysis must already be stored in D_sub, and this is part of the super-pathway pipeline, with a result already in 'comp'
#'  ... %>%
#'  mt_plots_equalizer(stat1='comp',
#'                     D2=D_sub,
#'                     stat2=='comp',
#'                     legend_fine="sub pathway",
#'                     legend_coarse='super pathway',
#'                     vline_fine = p.adj < 0.1,
#'                     vline_coarse = p.adj < 0.1) %>%
#'  ...}
#'
#' @author JK, MB
#'
#' @import ggplot2
#'
#' @export
mt_plots_equalizer <- function(D1,
                               stat1,
                               D2,
                               stat2,
                               legend_fine,
                               legend_coarse = NULL,
                               vline_fine = p.adj < 0.05,
                               vline_coarse = p.adj < 0.05,
                               color_list = c("#9494FF","red")) {

  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D1))
  stopifnot("SummarizedExperiment" %in% class(D2))
  # stopifnot(stat1!=stat2)

  # get results
  res1 <- mtm_get_stat_by_name(D1, stat1)
  res2 <- mtm_get_stat_by_name(D2, stat2)

  ##### Mustafa's code starts here

  # rm "name" duplicates
  frm_ <- function(x){
    xnames = dplyr::setdiff(colnames(x),"name")
    x[,c("name", xnames[!sapply(x[,xnames,drop =F],identical, x$name)]),drop =F]
  }

  # shortcuts
  rd1 = rowData(D1) %>% as.data.frame %>% frm_
  rd2 = rowData(D2) %>% as.data.frame %>% frm_


  # find field in first that contains all the second (e.g. all subpathways in the feature rowData)
  col2 <- rd2 %>% sapply(function(x) all((x %>% na.omit()) %in% rd1[["name"]]) )

  if (sum(col2)>1) stop(sprintf("Multiple columns in first rowData map to names of second: %s", paste0(colnames(rd2)[col2],collapse=", ")))
  if (sum(col2)==0) stop(sprintf("No columns in first rowData map to names of second. Specified correct rowData frames?"))


  col2 = col2 %>% which %>% names
  colnames(rd2)[colnames(rd2) == col2] = "COARSE"
  # if legend_coarse not given
  if(is.null(legend_coarse)) legend_coarse = col2

  colnames(rd1)[1] = "COARSE"
  colnames(rd2)[1] = "FINE"



  # df: data frame includes columns: "SUB_PATHWAY", "SUPER_PATHWAY", "statistic", "p.value"
  # name.df: primary key(column) name in df
  # df2: data.frame for super pathways, columns: SUPER_PATHWAY", "statistic", "p.value"
  # name.df2: foreign key(column) name in df2
  # th: log10(p.value) threholds for red dashed lines
  # colors: color_list for sub and super pathways

  mti_plot_equalizer_gg <- function(df, name.df="SUB", df2=NULL, name.df2="SUPER" ){

    # create x-axis
    df$x = abs(log10(df$p.value)) * sign(df$statistic)
    # x axis limits
    a = max(abs(df$x))

    # find x coordinates for cutoff lines
    xfine <- res2 %>% dplyr::filter(!!dplyr::enquo(vline_fine)) %>% .$p.value %>% max()
    xcoarse <- res1 %>% dplyr::filter(!!dplyr::enquo(vline_coarse)) %>% .$p.value %>% max()

    # main facetted plot
    gg<-
      ggplot(df, aes(x = x, y = FINE)) +
      geom_vline(xintercept = 0, color ="gray") +
      geom_vline(xintercept = c(-log10(xfine),log10(xfine)), color=color_list[1], alpha=0.4) +
      geom_vline(xintercept = c(-log10(xcoarse),log10(xcoarse)), color=color_list[2], alpha=0.4) +
      geom_point(pch = 22, fill = color_list[1], size = 3) +
      facet_grid(COARSE~. , scales = "free_y", space = "free_y") +
      theme(strip.background =element_rect(fill=NA),
            strip.text = element_text(colour = 'black', face = "bold"),
            strip.text.y = element_text(angle = 0, hjust = 0),
            panel.grid.major.y = element_line(color ="gray"),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            panel.background = element_rect(fill=NA, color ="black")) +
      scale_x_continuous(limits = c(-a,a))

    # add super pathways
    df.super = dplyr::summarise(dplyr::group_by(df,COARSE),
                                yy = mean(as.numeric(factor(FINE))),
                                xx = mean(x))

    if(!is.null(df2)){
      # if NA exist in super but in df2
      if(any(is.na(df.super$COARSE)) & !any(is.na(df2$COARSE))){
        df2 = rbind(df2, rep(NA,ncol(df2)))
      }

      df.super = dplyr::inner_join(df.super, df2, "COARSE")
      # create x-axis
      df.super$xx = abs(log10(df.super$p.value)) * sign(df.super$statistic)
    }


    gg = gg + geom_point(data = df.super, aes(x = xx,y = yy),pch = 22,
                         fill = color_list[2], size = 5, alpha = 0.7)

    # add legend
    df.legend = data.frame(x = rep(1,2), y = rep(NA,2), class = factor(c(name.df2, name.df),levels = c(name.df, name.df2)))
    gg + geom_point(data = df.legend,aes(x= x,y=y, fill =class),pch = 22, size = 4) +
      labs(fill="", x = expression(paste("directed log10(p)")), y = "") +
      scale_fill_manual(values = color_list) +
      theme(legend.position = "top", legend.key = element_blank(),
            legend.direction = "vertical", legend.justification = c(0,0))

  }



  p =  mti_plot_equalizer_gg(df = data.frame(rd2, res2), name.df = legend_fine,
                             df2 = data.frame(rd1, res1), name.df2 = legend_coarse )

  ## ADD AXIS GROUPS
  d <- mtm_get_stat_by_name(D1, stat1, fullstruct=T)
  if ("groups" %in% names(d) && length(d$groups)==2) {
    xlabel <- sprintf("%s high <--     directed log10(p)     --> %s high", d$groups[1], d$groups[2])
    p <- p + xlab(xlabel)
  }

  # add status information & plot
  funargs <- mti_funargs()
  D1 %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("equalizer plot of '%s' and '%s'", stat1, stat2),
      output = list(p)
    )

  # return
  D1


}









