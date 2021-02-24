#' 2D UMAP of samples
#'
#' Create a UMAP plot of samples. Can be colored by any variable in colData.
#'
#' @param D \code{SummarizedExperiment} input.
#' @param title Title of the plot. Default: "UMAP".
#' @param scale_data  Scale data before plotting? Mean 0, std 1. Default: F.
#' @param label_col Column to use as label. Default: ''.
#' @param text_repel Try to avoid all text overlaps when labeling? Default:T
#' @param store_matrices Store scores in result structure? Default: F
#' @param n_neighbors The size of the local neighborhood UMAP will look at when attempting to learn the meanifold structure
#'    of the data. Default: 15.
#' @param ggadd Further elements/functions to add (+) to the ggplot object. Default: NULL.
#' @param ... Additional expression directly passed to aes() of ggplot. Can refer to colData.
#'
#' @return result output: plot(s)
#' @return result output2: scores and loadings matrix
#
#' @examples
#' \dontrun{## UMAP on scale_data, color and shape by "Group" variable in colData
#' ... $>$ mt_plots_umap(scale_data=T, color=Group, shape=Group, title="UMAP - scaled data") %>% ...
#' }
#'
#' @author JK
#'
#' @import ggplot2
#'
#' @export
mt_plots_umap <- function(D,
                          title="UMAP",
                          scale_data=F,
                          label_col='',
                          text_repel=T,
                          store_matrices=F,
                          n_neighbors=15,
                          ggadd=NULL,
                          ...) {

  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))

  # get argument names from dots
  n <- sapply(as.list(substitute(list(...)))[-1L], deparse)
  dot_args <- names(n)
  
  # check for defunct argument names
  if ("scaledata" %in% dot_args) stop("You used the old MT naming convention scaledata Should be: scale_data")
  if ("labelby" %in% dot_args) stop("You used the old MT naming convention labelby Should be: label_col")
  if ("textrepel" %in% dot_args) stop("You used the old MT naming convention textrepel Should be: text_repel")
  if ("store.matrices" %in% dot_args) stop("You used the old MT naming convention store.matrices Should be: store_matrices")
  
  # helper function to combine two aesthetics, e.g. from aes() and aes_string()
  combine_aes <- function(...) {
    v <- c(...)
    class(v) <- "uneval"
    v
  }

  # extract data and verify that there are no NA values
  X = t(assay(D))
  if (any(is.na(X))) stop("Data matrix for UMAP cannot contain NAs")

  # scale?
  if (scale_data) X <- scale(X) #By default, the scale R-function: mean-centers and scales to unit variance the X matrix

  # UMAP
  um <- umap::umap(d=as.matrix(X), n_neighbors=n_neighbors)

  # assemble data frame, two components and sample info
  df = data.frame(x = um$layout[,1], y = um$layout[,2], colData(D)) # scores and colData
  colnames(df)[1:2] <- c("comp1","comp2")

  # plot
  p <- ggplot(data=df) +
    geom_point(combine_aes(aes_string(x="comp1",y="comp2"),aes(...))) +
    xlab("comp 1") + ylab("comp 2") + ggtitle(title)
  # add text?
  if (nchar(label_col)>0) {
    if (text_repel) p <- p + ggrepel::geom_text_repel(aes_string(x="comp1",y="comp2",label=label_col,...))
    else p <- p + geom_text(combine_aes(aes_string(x="comp1",y="comp2",label=label_col),aes(...)))
  }

  # add custom elements?
  if (!is.null(ggadd)) p <- p+ggadd

  # prep output matrices
  if (store_matrices) {
    output2 <- list(layout=umap$layout)
  } else {
    output2 <- NULL
  }

  # add status information & plot
  funargs <- mti_funargs()
  metadata(D)$results %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("UMAP, label_col: %s, aes: %s", label_col,  mti_dots_to_str(...)),
      output = list(p),
      output2 = output2
    )

  # return
  D


}
