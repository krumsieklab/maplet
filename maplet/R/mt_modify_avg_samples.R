#' Averages duplicate samples
#'
#' Averages the values of samples (columns) with the same values in specified colData column.
#'
#' @param D \code{SummarizedExperiment} input.
#' @param group_col Name of colData column (sample annotation) by which duplicates can be identified.
#'
#' @return assay: Duplicate samples (columns) combined.
#'
#' @examples
#' \dontrun{... %>% mt_modify_avg_samples(group_col = "RID") %>% ...}
#'
#' @author AS
#'
#' @export
mt_modify_avg_samples <- function(D, group_col) {

  stopifnot("SummarizedExperiment" %in% class(D))
  if(missing(group_col))
    stop("group_col can't be empty")

  # TODO: check that coldata for duplicates are the same, throw mti_logwarning if not
  X <- t(assay(D))
  ave_col <- colData(D)[[group_col]]
  u.col <- sort(unique(ave_col))
  counts <- as.numeric(table(ave_col))
  dup <- unique(ave_col[ave_col%in%u.col[counts > 1]])
  to_remove <- c()
  for (i in 1:length(dup)) {
    index <- which(ave_col == dup[i])
    dat.slice <- X[index, ]
    X[index, ] <- sapply(colMeans(dat.slice, na.rm = TRUE), rep, length(index))
    to_remove <- c(to_remove,index[2:length(index)])
  }
  assay(D) <- t(X)
  D<- D[,-to_remove]

  ## add status information
  funargs <- mti_funargs()
  D %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf(' %d duplicate samples combined', length(to_remove))
    )

  ##return

  D
}
