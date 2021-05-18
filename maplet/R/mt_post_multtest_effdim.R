#' Multiple Testing Correction using Effective Number of Independent Tests
#'
#' Adjust output of statistical test using the method described in Li J, Ji L (2005) Adjusting multiple testing in
#' multilocus analyses using the eigenvalues of a correlation matrix. Heredity 95:221-227
#' \href{https://www.ncbi.nlm.nih.gov/pubmed/16077740}{https://www.ncbi.nlm.nih.gov/pubmed/16077740}.
#'
#' @param D \code{SummarizedExperiment} input.
#' @param stat_name Name of the statistical comparison to adjust.
#' @param p_col Name of p-value column in statistical table. Default: p.value.
#'
#' @return $results[[stat_name]]$output: p.adj column added to statistical table.
#'
#' @examples
#' \dontrun{# correct the statistical comparison "comp" using effective dimension method
#' ... %>%
#'  mt_post_multtest_effdim(stat_name="comp") %>% ...}
#'
#' @author KC
#'
#' @export
mt_post_multtest_effdim <- function(D,
                                    stat_name,
                                    p_col = p.value) {
  # validate argument
  if(missing(stat_name)){
    stop("stat_name must be provided")
  }

  # check no missing values in assay
  if(any(is.na(assay(D)))){
    stop("NA values in data")
  }

  # get entry
  stat_id <- metadata(D)$results %>%
    purrr::map_lgl(~"stats" %in% .x$fun && .x$output$name == stat_name) %>%
    which()
  if(length(stat_id) == 0)
    stop("stat element with name ", stat_name, " does not exist")
  if(length(stat_id)  > 1)
    stop("there are multiple stat elements with name ", stat_name)

  # calculate effective number of independent tests
  effdim <- calc_effdim(t(assay(D)))

  p_col <- dplyr::enquo(p_col)

  # perform correction
  metadata(D)$results[[stat_id]]$output$table %<>%
    dplyr::mutate(p.adj=(!!p_col * effdim) %>%
                    pmin(.,1))

  ## add status information & plot
  funargs <- mti_funargs()
  D %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("Multiple testing correction of '%s' using effective dimension method", stat_name),
      output = list(effdim = effdim)
    )

  ## RETURN
  D

}

#' Calculate Effective Number of Independent Tests
#'
#' Calculate the effective number of independent tests using the method described in Li J, Ji L (2005) Adjusting multiple
#' testing in multilocus analyses using the eigenvalues of a correlation matrix. Heredity 95:221-227
#'
#' @description
#' NA values are not allowed.
#'
#' @param X Dataframe with samples as rows and features as columns.
#'
#' @return effective number of idependent tests
#'
#' @noRd
calc_effdim <- function(X){

  correlationmatrix <- as.data.frame(abs(cor(X)))

  evals<-eigen(t(correlationmatrix),symmetric=T)$values

  var<-var(evals)
  M<-length(evals)
  L<-(M-1)

  newevals<-evals
  for(i in 1:length(newevals)) {
    if(newevals[i] < 0) {
      newevals[i] <- 0
    }
  }

  IntLinewevals<-newevals

  for(i in 1:length(IntLinewevals)) {
    if(IntLinewevals[i] >= 1 ) {
      IntLinewevals[i] <- 1
    }
    if(IntLinewevals[i] < 1 ) {
      IntLinewevals[i] <- 0
    }
  }

  NonIntLinewevals <- newevals-floor(newevals)

  effdim <- sum(NonIntLinewevals+IntLinewevals)

  effdim

}
