#' Calculate coefficient of variation (CV)
#'
#' Use QC samples to calculate coefficient of variation (CV) and add values as new column to rowData.
#'
#' @param D \code{SummarizedExperiment} input.
#' @param qc_samples Logical expression. Can use columns from colData.
#' @param out_col Name of new rowData column to store the cv values.
#' @param replicates OPTIONAL. Flag to indicate if the selected samples are replicates. Default: F.
#' @param id_col OPTIONAL. If replicates is T, it is the name of the column in colData containing sample IDs. Default: NULL.
#'
#' @return rowData: New annotation column added.
#'
#' @examples
#' \dontrun{... %>% mt_pre_cv(qc_samples=="PQC", out_col = "PQC_cv") %>% ...}
#'
#' @author AS, RB
#'
#' @export
mt_pre_cv <- function(D, qc_samples, out_col, replicates=F, id_col=NULL){

  stopifnot("SummarizedExperiment" %in% class(D))

  if(missing(qc_samples)) stop("qc_samples can't be empty")

  ## APPLY FILTER TO ROW DATA
  qc_samples_q <- dplyr::enquo(qc_samples)
  cd <- colData(D) %>%
    data.frame(row.names = 1:ncol(D)) %>%
    tibble::rownames_to_column("colnames") %>%
    dplyr::filter(!!qc_samples_q)

    ## SUBSET SAMPLES
  D_cv <- D[, as.numeric(as.matrix(cd$colnames))]
  # calc_cv function
  calc_cv <- function(x){sd(x, na.rm = T)/mean(x, na.rm = TRUE)}
  replace_zeros <- function(x){ x[which(x==0)] <- NA; return(x)}
  ## Calculate feature cv scores
  if (replicates) {
    if(missing(id_col)) stop(sprintf('id_col must be provided if cv is to be calcualted on replicates!'))
    # get the assay of duplicates with id column
    cv_data <- dplyr::bind_cols(D_cv %>% assay() %>% t() %>% data.frame(),
         D_cv %>% colData() %>% data.frame() %>% dplyr::select(!!id_col) %>% data.frame())
    # mean per ID
    rowData(D)[[out_col]] <- cv_data %>%
      dplyr::group_by_at(vars(starts_with(!!id_col))) %>% # calculate cv per duplicated ID
      dplyr::summarise_at(vars(-starts_with(!!id_col)), calc_cv) %>% select(-!!id_col) %>%
      summarise_all(replace_zeros) %>% colMeans(na.rm = T) %>% unlist()

  } else {
    rowData(D)[[out_col]] <- apply(assay(D_cv), 1, calc_cv)
  }

  ## add status information
  funargs <- maplet:::mti_funargs()
  metadata(D)$results %<>%
    maplet:::mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("Added QC cv to rowData")
      )
  ## return
  D
}
