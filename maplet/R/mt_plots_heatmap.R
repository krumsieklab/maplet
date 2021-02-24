#' Heatmap plot by pheatmap
#'
#' Creates a heatmap using the function \code{pheatmap] from the \code{pheatmap] package. All \code{pheatmap::pheatmap}
#' arguments can be passed \code{\href{https://github.com/raivokolde/pheatmap}{https://github.com/raivokolde/pheatmap}}.
#'
#' @param D \code{SummarizedExperiment} object.
#' @param scale_data Scaling the data. Default: T.
#' @param sym_zero Make color scale symmetric around 0? (should only be used for scaled data). Default: F.
#' @param fun Function to transform/scale \code{t(assay(D))}, ie \code{mat = fun(t(assay(D)))} will be plotted.
#' @param silent Don't draw the table? A pheatmap argument, maplet uses a different default. Default: T.
#' @param return_gg Should pheatmap object be converted to gg object. Default: T.
#' @param gg_scale Scaling of plot to be converted to gg object.
#' @param annotation_col Data frame that specifies the annotations shown on the top of the heatmap.
#' @param annotation_row Data frame that specifices the annotations shown on the left side of the heatmap.
#' @param gg_ymin Minimum coordinate for ggplot y-axis. Default: 1-gg_scale.
#' @param gg_xmin Minimum coordinate for ggplot x-axis. Default: 1-gg_scale.
#' @param gg_xmax Maximum coordinate for ggplot x-axis. Default: gg_scale.
#' @param gg_ymax Maximum coordinate for ggplot y-axis. Default: gg_scale.
#' @param color Vector of colors used in heatmap.
#' @param ggadd  Further elements/functions to add (+) to the ggplot object.
#' @param \dots  See \code{pheatmap::pheatmap} for pheatmap arguments.
#'
#' @return $results$output: plot, heatmap
#'
#' @author MB
#'
#' @examples
#' \dontrun{D %>%
#' mt_plots_heatmap(annotation_row = c("SUPER_PATHWAY", "PLATFORM", "RI"),
#'                  annotation_col = c("GROUP_DESC","BATCH_MOCK","gender"),
#'                  fun = function(x) scale(exp(scale(x))),
#'                  clustering_distance_cols =  "correlation",
#'                  clustering_distance_rows = "minkowski")
#'  }
#'
#' @import ggplot2
#' @import pheatmap
#'
#' @export
mt_plots_heatmap <- function(D,
                             scale_data=F,
                             sym_zero=F,
                             fun = function(x){ if(scale_data) return(scale(x)); x},
                             silent = TRUE,
                             return_gg = T,
                             annotation_col = NA,
                             annotation_row = NA,
                             gg_scale = 1,
                             gg_ymin = 1 - gg_scale,
                             gg_xmin = 1 - gg_scale,
                             gg_xmax = gg_scale,
                             gg_ymax = gg_scale,
                             ggadd=NULL,
                             color = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100),
                             ...){

  # upon Jan's comment annotation_col and annotation_row are swapped for compatibility with SummarizedExperiment

  # get all inputs
  aa = c(as.list(environment()), list(...))
  
  # get argument names from dots
  n <- sapply(as.list(substitute(list(...)))[-1L], deparse)
  dot_args <- names(n)
  
  # check for defunct argument names
  if ("scaledata" %in% dot_args) stop("You used the old MT naming convention scaledata. Should be: scale_data")
  if ("sym0" %in% dot_args) stop("You used the old MT naming convention sym0. Should be: sym_zero")
  if ("fD" %in% dot_args) stop("You used the old MT naming convention fD. Should be: fun")
  if ("gg.scale" %in% dot_args) stop("You used the old MT naming convention gg.scale. Should be: gg_scale")
  if ("gg.ymin" %in% dot_args) stop("You used the old MT naming convention gg.ymin. Should be: gg_ymin")
  if ("gg.ymax" %in% dot_args) stop("You used the old MT naming convention gg.ymax. Should be: gg_ymax")
  if ("gg.xmin" %in% dot_args) stop("You used the old MT naming convention gg.xmin. Should be: gg_xmin")
  if ("gg.xmax" %in% dot_args) stop("You used the old MT naming convention gg.xmax. Should be: gg_xmax")
  if ("return.gg" %in% dot_args) stop("You used the old MT naming convention return.gg. Should be: return_gg")

  # fun(t(assay(D))) will be heatmapped
  x = t(assay(D))
  if (any(is.na(x))) stop("Data matrix for heatmap cannot contain NAs")

  # if x doesn't have rownames, use numbers 1:nrow
  if (is.null(rownames(x))) rownames(x) <- 1:nrow(x)
  # store
  x.colnames = colnames(x)
  x.rownames = rownames(x)
  aa$mat = fun(x)
  aa$color = color
  
  # keep only pheatmap::pheatmap parameters
  aa = aa[!(names(aa) %in% c("D","scale_data","fun","return_gg", "gg_scale", "gg_ymin", "gg_xmin", "gg_xmax", "gg_ymax"))]
  # annotations will be added later
  aa = aa[!(names(aa) %in% c("annotation_col", "annotation_row"))]

  # deprecated pheatmap parameter 'annotation', see pheatmap::pheatmap
  aa$annotation = NA

  # annotate rows with given variables in attr(D, "colData")
  if(!is.na(annotation_col[1])){
    annotation_col = D %>% colData %>% as.data.frame %>% `[`(,annotation_col,drop = F)
    rownames(annotation_col) = x.rownames
    aa$annotation_row = annotation_col
  }

  # annotate columns with given variables in attr(D, "rowData")
  if(!is.na(annotation_row[1])){
    annotation_row = D %>% rowData %>% as.data.frame %>% `[`(,annotation_row,drop = F)
    rownames(annotation_row) = x.colnames
    aa$annotation_col = annotation_row
  }

  # symmetric around zero?
  if (sym_zero) {
    cap <- max(abs(aa$mat)) # from -cap to +cap
    cs = length(aa$color) # number of color steps
    aa$breaks = seq(-cap, cap, 2*cap/cs) # technical stuff
  }

  # plot pheatmap
  re <- do.call(pheatmap::pheatmap, aa)

  # cast pheatmap object to gg object
  if(return_gg){
    re <- ggplot(data.frame(x = 0:1, y = 0:1), aes_(x = ~x, y = ~y)) +
      geom_blank() + scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
      scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
      annotation_custom(re$gtable, xmin = gg_xmin, xmax = gg_xmax, ymin = gg_ymin, ymax = gg_ymax) +
      theme_void()

    # add custom elements?
    if (!is.null(ggadd)) re <- re+ggadd
  }


  # add status information & plot
  funargs <- mti_funargs()
  metadata(D)$results %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("Heatmap of assay data."),
      output = list(re)
    )

  # return
  D
}


