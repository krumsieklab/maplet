#' Plot dilution factors
#'
#' @description
#' Correlate variable with dilution factors from quotient normalization.
#'
#' @description
#' \itemize{
#' \item for factors: will produce boxplot
#' \item for quantitative variable: will produce scatter plot
#' \item requires exactly one previous quotient normalization call in the pipeline
#' }
#'
#' @param D \code{SummarizedExperiment} input.
#' @param in_col Name of colData column (sample annotation) to correlation dilution factors with.
#' @param boxplot  Produce boxplot (TRUE), or beeswarm plot (FALSE). Only relevant for factor comparison. Default: T.
#' @param ggadd Further elements / functions to add (+) to the ggplot object. Default: NULL.
#'
#' @return $result$output: plot, comparison with dilution factors
#'
#' @examples
#' \dontrun{%>% mt_plots_dilution_factor(in_col="group") %>% # compare with 'group' sample annotation
#' }
#'
#' @author JK
#'
#' @import ggplot2
#'
#' @export
mt_plots_dilution_factor <- function(D,
                              in_col,
                              boxplot=T,
                              ggadd=NULL) {

  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  stopifnot(length(in_col)==1)

  # validate that there has been exactly one quotient normalization step
  q <- D%>% mtm_res_get_entries(c("pre","norm","quot"))
  if (length(q)>1) stop("There has been more than one quotient normalization call.")
  if (length(q)==0) stop("No quotient normalization performed.")
  # get dilution factors
  vd = q[[1]]$output$dilution

  # get variable to compare to
  if (!(in_col %in% colnames(colData(D)))) stop(sprintf("'%s' not found in sample annotations.", in_col))
  vc = colData(D)[[in_col]]

  # generate data frame for plotting
  dfplot <- data.frame(
    dilution.factor = vd,
    y = vc
  )
  colnames(dfplot)[2] <- in_col

  # either produce boxplot or scatter plot
  if (is.character(vc) || is.factor(vc)) {
    # ensure factor
    dfplot[[in_col]] = as.factor(dfplot[[in_col]])
    # boxplot
    p <- dfplot %>%
      ggplot(aes_string(x=in_col,y="dilution.factor",color=in_col)) +
      ggtitle("quotient normalization dilution factors")
    if (boxplot)
      p <- p+geom_boxplot()
    if (!boxplot)
      p <- p+ggbeeswarm::geom_quasirandom()

  } else {
    if (!is.numeric(vc)) stop(sprintf("'%s' has to be character, factor, or numeric.", in_col))
    # scatter
    p <- dfplot %>% ggplot() +
      geom_point(aes_string(x=in_col,y="dilution.factor")) +
      ggtitle("quotient normalization dilution factors")
  }

  if (!is.null(ggadd)) p <- p+ggadd

  # add status information & plot
  funargs <- mti_funargs()
  D %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("dilution factor plot, '%s'",in_col),
      output = list(p)
    )

  # return
  D

}

