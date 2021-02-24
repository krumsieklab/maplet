#' 2D PCA of samples
#'
#' Create PCA plot of samples. Can be colored by any variable in colData.
#'
#' @param D \code{SummarizedExperiment} input.
#' @param title Title of the plot. Default: "PCA".
#' @param scale_data  Scale data before plotting? Scale to mean 0, std 1. Default: F.
#' @param pc1 First PC to plot. Default: 1 (PC1).
#' @param pc2 Second PC to plot. Default: 2 (PC2).
#' @param data_type Data of type 'scores' or 'loadings'. Default: 'scores'.
#' @param label_col Column to use for labels. Default: ''.
#' @param text_repel Try to avoid all text overlaps when labeling? Default:T.
#' @param ellipse Confidence interval for ellipse. Default: NA (no ellipse).
#' @param exp_var_plot Add explained variance plot? Default: F.
#' @param store_matrices Store scores and loadings matrices in result structure? Default: F.
#' @param ggadd Further elements/functions to add (+) to the ggplot object. Default: NULL.
#' @param ... Additional expression directly passed to aes() of ggplot, can refer to colData.
#'
#' @return results$output: plot(s)
#' @return results$output2: scores and loadings matrix
#
#' @examples
#' \dontrun{## PCA on scale_data, color and shape by "Group" variable in colData
#' ... $>$ mt_plots_pca(scale_data=T, color=Group, shape=Group, title="PCA - scaled data") %>% ...
#' ## PCA scores plot on non-scaled data, with ellipse and extra explained variance plot, and two ggadds (white background and centering of title)
#' mt_plots_pca(title="PCA scores", data_type = 'scores', scale_data=F, pc1=1, pc2=2, ellipse=0.95, exp_var_plot=T, ggadd = theme_bw() + theme(plot.title=element_text(hjust=0.5)))
#' }
#'
#' @author JK
#'
#' @import ggplot2
#'
#' @export
mt_plots_pca <- function(D,
                         title="PCA",
                         scale_data=F,
                         pc1=1,
                         pc2=2,
                         data_type='scores',
                         label_col='',
                         text_repel=T,
                         ellipse=NA,
                         exp_var_plot=F,
                         store_matrices=F,
                         ggadd=NULL,
                         ...) {

  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))

  # get argument names from dots
  n <- sapply(as.list(substitute(list(...)))[-1L], deparse)
  dot_args <- names(n)
  
  # check for defunct argument names
  if ("scaledata" %in% dot_args) stop("You used the old MT naming convention scaledata Should be: scale_data")
  if ("PCa" %in% dot_args) stop("You used the old MT naming convention PCa Should be: pc1")
  if ("PC1" %in% dot_args) stop("You used the old MT naming convention PC1 Should be: pc1")
  if ("PCb" %in% dot_args) stop("You used the old MT naming convention PCb Should be: pc2")
  if ("PC2" %in% dot_args) stop("You used the old MT naming convention PC2 Should be: pc2")
  if ("showscore" %in% dot_args) stop("You used the old MT naming convention showscore Should be: data_type")
  if ("show" %in% dot_args) stop("You used the old MT naming convention show Should be: data_type")
  if ("labelby" %in% dot_args) stop("You used the old MT naming convention labelby Should be: label_col")
  if ("textrepel" %in% dot_args) stop("You used the old MT naming convention textrepel Should be: text_repel")
  if ("expvarplot" %in% dot_args) stop("You used the old MT naming convention expvarplot Should be: exp_var_plot")
  if ("store.matrices" %in% dot_args) stop("You used the old MT naming convention store.matrices Should be: store_matrices")
  
  # helper function to combine two aesthetics, e.g. from aes() and aes_string()
  combine_aes <- function(...) {
    v <- c(...)
    class(v) <- "uneval"
    v
  }

  # extract data and verify that there are no NA values
  X = t(assay(D))
  if (any(is.na(X))) stop("Data matrix for PCA cannot contain NAs")
  # check that data_type is either "scores" or "loadings"
  if (!(data_type %in% c("scores","loadings"))) stop("Show must be either 'scores' or 'loadings'")

  # scale?
  if (scale_data) X <- scale(X) #By default, the scale R-function: mean-centers and scales to unit variance the X matrix

  # PCA
  pca <- stats::prcomp(x=as.matrix(X), center=F, scale=F)
  expvar = (pca$sdev)^2 / sum(pca$sdev^2)
  # assemble data frame, two PCS and sample info
  if (data_type=='scores') df = data.frame(x = pca$x[,pc1], y = pca$x[,pc2], colData(D)) # scores and colData
  else  df = data.frame(x = pca$rotation[,pc1], y = pca$rotation[,pc2], rowData(D)) # loadings and rowData
  colnames(df)[1:2] <- c(sprintf("PC%d",pc1),sprintf("PC%d",pc2))

  # prep
  pc1name = sprintf('PC%d (%.1f%%)', pc1, expvar[pc1]*100)
  pc2name = sprintf('PC%d (%.1f%%)', pc2, expvar[pc2]*100)

  # plot
  p <- ggplot(data=df, combine_aes(aes_string(x=sprintf("PC%d",pc1),y=sprintf("PC%d",pc2)),aes(...))) +
    geom_point() +
    xlab(pc1name) + ylab(pc2name) + ggtitle(title)
  # add text?
  if (nchar(label_col)>0) {
    if (text_repel) p <- p + ggrepel::geom_text_repel(aes_string(x=sprintf("PC%d",pc1),y=sprintf("PC%d",pc2),label=label_col,...))
    else p <- p + geom_text(combine_aes(aes_string(x=sprintf("PC%d",pc1),y=sprintf("PC%d",pc2),label=label_col),aes(...)))
  }
  # add ellipse?
  if (!is.na(ellipse)) p <- p + stat_ellipse(level=ellipse)
  # add custom elements?
  if (!is.null(ggadd)) p <- p+ggadd


  # add explained variance plot?
  plotlist <- list(p)
  if (exp_var_plot) {
    nlimit <- 10
    expdf <- data.frame(n=1:length(expvar), expvar=expvar*100)
    # cut to first nlimit at most
    expdf <- expdf[1:min(nlimit,nrow(expdf)),]
    newp <- ggplot(data=expdf, aes(x=n, y=expvar)) +
      geom_bar(stat="identity", fill="steelblue")+
      ylab("") + xlab("PC") + ggtitle(sprintf("Explained variance (first %d PCs)", nrow(expdf))) +
      geom_text(aes(label=sprintf('%.1f%%',expdf$expvar)), vjust=-0.3, size=3.5) +
      scale_x_discrete(limits=1:nrow(expdf))
    # add to list of plots for results
    plotlist[[2]] <- newp
  }

  # prep output matrices
  if (store_matrices) {
    scores=pca$x
    loadings=pca$rotation
    rownames(loadings) <- D %>% rowData %>% .$name
    output2 <- list(scores=scores, loadings=loadings)
  } else {
    output2 <- NULL
  }

  # add status information & plot
  funargs <- mti_funargs()
  metadata(D)$results %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("PCA, %s, label_col: %s, aes: %s", data_type, label_col,  mti_dots_to_str(...)),
      output = plotlist,
      output2 = output2
    )

  # return
  D


}
