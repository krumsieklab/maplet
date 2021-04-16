#' Correct by reference
#'
#' Correct data by a reference (e.g. NIST).
#'
#' @param D \code{SummarizedExperiment} input
#' @param qc_samples Logical expression. Can use fields from \code{colData()}
#' @param plate_col name of the colData column representing plate number
#'
#' @return D QCed
#'
#' @examples
#' \dontrun{... %>% mt_pre_reference_correction(qc_samples=Sample.Identification==102,
#' plate_col='Plate.Bar.Code') %>% ...}
#'
#' @author RB
#'
#' @export
mt_pre_reference_correction <- function(D, qc_samples, plate_col){

  stopifnot("SummarizedExperiment" %in% class(D))
  if(missing(qc_samples)) stop("qc_samples can't be empty!")
  if(missing(plate_col)) stop("plate_col can't be empty!")

  ## APPLY FILTER TO ROW DATA
  qc_samples_q <- dplyr::enquo(qc_samples)
  cd <- colData(D) %>%
    data.frame(row.names = 1:ncol(D)) %>%
    tibble::rownames_to_column("colnames") %>%
    dplyr::filter(!!qc_samples_q)

  ## QC and non-QC part
  D1 <- D[, as.numeric(as.matrix(cd$colnames))]
  D <- D[, -as.numeric(as.matrix(cd$colnames))]

  ## Compute correction on qc samples
  my_mean <- function(x){mean(x, na.rm = T)}
  replace_nan <- function(x){ x[is.nan(x)] <- 1; return(x)}
  # average across all plates per metabolite
  qc_avg <- colMeans(D1 %>% assay() %>% t() %>% data.frame(), na.rm = TRUE)
  # average per plate divided by average across all plates
  qc_val <- dplyr::bind_cols(D1 %>% assay() %>% t() %>% data.frame(),
                    D1 %>% colData() %>% data.frame() %>% select(!!plate_col)) %>%
            dplyr::group_by(!!as.name(plate_col)) %>%
            dplyr::summarise_at(vars(-starts_with(!!plate_col)), my_mean) %>%
            ungroup() %>% data.frame()
  row_names <- qc_val %>% pull(!!plate_col)
  qc_val <- apply((qc_val %>% select(-!!plate_col)), 1, FUN=function(x) x/ qc_avg) %>%
    t() %>% data.frame() %>% summarise_all(replace_nan)
  rownames(qc_val) <- row_names
  ## Correct the data using qc values
  X <- D %>% assay() %>% t()
  plates <- D %>% colData() %>% data.frame() %>% pull(!!plate_col)
  for(i in c(1:nrow(X))){
    X[i, ] <- as.numeric(as.matrix(X[i, ] / qc_val[rownames(qc_val)%in%plates[i], ]))
  }
  assay(D) <- t(X)
  ## add status information
  funargs <- mti_funargs()
  D %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("Corrected assay data based on NIST samples!")
    )
  ## return
  D
}
