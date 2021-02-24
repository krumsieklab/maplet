#' Box Plot of samples
#'
#' Create a box plot of the samples. Can be colored by factor.
#'
#' @param D \code{SummarizedExperiment} input.
#' @param title Title of box plot. Default: "Sample boxplot".
#' @param show_legend Show legend? Default: T.
#' @param ylabel The y-axis label. Default: "Feature concentrations".
#' @param plot_logged Show plot logged? Note: plot will still be logged if data have been logged before in pipeline. Default: F.
#' @param ggadd Further elements/functions to add (+) to the ggplot object. Default: NULL.
#' @param ...  Additional arguments directly passed to aes() of ggplot.
#'
#' @return $results$output: plot, sample box plot
#'
#' @examples
#' \dontrun{## sample boxplot, color by colData 'group' variable, with specific title, on log scale,
#' ... %>% mt_plots_sample_boxplot(color=group, title='after quotient normalization', plot_logged=T) %>% ...
#' }
#'
#' @author JK
#'
#' @import ggplot2
#'
#' @export
mt_plots_sample_boxplot <- function(D,
                                    title="Sample boxplot",
                                    show_legend=T,
                                    ylabel = "Feature concentrations",
                                    plot_logged=F,
                                    ggadd=NULL,
                                    ...) {

  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))

  # get argument names from dots
  n <- sapply(as.list(substitute(list(...)))[-1L], deparse)
  dot_args <- names(n)
  
  # check for defunct argument names
  if ("plottitle" %in% dot_args) stop("You used the old MT naming convention plottitle. Should be: title")
  if ("legend" %in% dot_args) stop("You used the old MT naming convention legend. Should be: show_legend")
  if ("manual_ylab" %in% dot_args) stop("You used the old MT naming convention manual_ylab. Should be: ylabel")
  if ("manual_ylabel" %in% dot_args) stop("You used the old MT naming convention manual_ylabel. Should be: ylabel")
  if ("logged" %in% dot_args) stop("You used the old MT naming convention logged. Should be: plot_logged")
  
  # plot_logged?
  Dplot = D
  if (plot_logged) {
    assay(Dplot) <- log2(assay(Dplot))
    ylabel = sprintf("%s [log2]", ylabel)
  }

  # merge with sample annotations, only keep the ones that were actually used
  cd <- Dplot %>% colData() %>% as.data.frame() %>% tibble::rownames_to_column("merge.primary")
  keep <- c(mti_extract_variables(quos(...)), "merge.primary")
  cd <- cd[,colnames(cd) %in% keep,drop=F]
  df <- cbind(cd, t(assay(Dplot)))
  # generate ggplot
  p <- df %>%  tidyr::gather(metab, value, dplyr::one_of(rownames(Dplot))) %>%
    ggplot(aes(x = merge.primary, y = value, ...)) +
    geom_boxplot() +
    ylab(ylabel) +
    ggtitle(title) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  # todo add rowname

  # remove legend?
  if (!show_legend) p = p + theme(legend.position="none")

  # add custom elements?
  if (!is.null(ggadd)) p <- p+ggadd

  # add status information & plot
  funargs <- mti_funargs()
  metadata(D)$results %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("sample boxplot, aes: %s", mti_dots_to_str(...)),
      output = list(p)
    )

  # return
  D


}
