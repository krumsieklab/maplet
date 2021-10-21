#' Confounding correction using stepwise AIC
#'
#' Each feature is corrected by lm with the confounding variables determined by stepwise aic and residuals are returned
#' in lieu of uncorrected features expressions.
#'
#' @param D \code{SummarizedExperiment} input.
#' @param cols_to_correct column names from colData to correct for.
#' @param cols_to_exclude optional, column names from colData to exclude for correction.
#' @param n_cores Number of cores to use in parallelization. Default: 1.
#'
#' @return assay: Corrected data.
#' @return $results$output: Returns data.frame with the feature, its covars and the fit pval and rsq values.
#'
#' @examples
#'  \dontrun{#... %>% mt_pre_confounding_correction_stepaic(cols_to_correct = c(1, 4, 5), n_cores = 10)
#'  }
#'
#'
#' @author AS, RB
#'
#' @export

mt_pre_confounding_correction_stepaic <- function(D, cols_to_correct, cols_to_exclude=NULL, n_cores = 1) {
  
  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  
  # coldata for the SE
  Y <- D %>% colData() %>% data.frame()
  # all column names in colData
  all_cols <- Y %>% names()
  
  # identify the columns with single level to be excluded
  # exclude the covariates with constant values for all patients
  noNA <- function(x){x <- x[!is.na(x) & x!="NA" & x!=""];x}
  for (icol in (cols_to_correct)) {
    # just one factor represented?
    if (Y %>% .[[icol]] %>% unique() %>% noNA() %>% length() <= 1) {
      cols_to_exclude <- c(cols_to_exclude, icol)
    } 
  }   
  # remove the col number of covariates to exclude
  if(length(cols_to_exclude)>0){
    cols_to_correct <- setdiff(cols_to_correct, cols_to_exclude)
    cols_to_exclude <- toString(cols_to_exclude)
    warning(sprintf("Covariates in column numbers %s with constant values were removed from correction!", cols_to_exclude))
  }
  if(length(cols_to_correct)<1) stop("No covariates left to correct for!")
  
  # which samples have atleast one NA covariate information ?
  samples_to_exclude <- apply(Y [, cols_to_correct], 1, FUN=function(x)
    length(which(is.na(x)))>0) %>% which()
  # subset the input based on missing covariate info
  if((length(samples_to_exclude) > 0)) {
    if(length(samples_to_exclude) < nrow(Y)){
      D <- D [, -samples_to_exclude]  
    } else stop("No samples with any covariate info!")
  } 
  
  # loop over analytes
  X <- t(assay(D))
  # bind covariate columns with features for modelling
  model_data <- cbind.data.frame(X, colData(D)[, cols_to_correct])
  outlist <- parallel::mclapply(1:ncol(X), FUN=function(i) {
    # all covariates and one feature
    form <- paste(colnames(X)[i], "~", paste(cols_to_correct, collapse=" + "))
    # compute coefficients
    mod <- stats::lm(stats::as.formula(form), data = model_data)
    # stepwise selection of covariates that have an effect on feature
    stepmod <- MASS::stepAIC(mod, direction = "backward", trace=0, k = log(nrow(model_data)))
    # selected covariates
    selected_covars <- as.character(stats::formula(stepmod)[3])
    # if no covariates are selected...
    if(selected_covars==1){
      met_vec <- X[, i] # ...return the original feature value
      selected_covars <- fit_pval <- fit_rsq <- NA
    } else { # otherwise
      # fetch the names of the covariates
      selected_covars <- unlist(strsplit(selected_covars, split="+", fixed=TRUE))
      # collapse separated by semicolon
      selected_covars <- paste(gsub(" ", "", selected_covars), collapse=";")
      # get the significance of the model
      fit_pval <- signif(stats::pf(summary(stepmod)$fstatistic[1], summary(stepmod)$fstatistic[2], summary(stepmod)$fstatistic[3], lower.tail = FALSE), 3)
      # rsq of the model
      fit_rsq <- signif(summary(stepmod)$r.squared, 3)
      # return the residual --> corrected feature values
      met_vec <- stepmod$residuals + stepmod$coefficients["(Intercept)"]
    }
    # result vector
    res <- list(met_vec, list(colnames(X)[i], selected_covars, fit_pval, fit_rsq))
    return(res)
  }, mc.cores= n_cores)
  
  ###### results gathering : bind outputs into assay for summarized experiment
  covar_adjusted <- do.call(rbind, lapply(outlist, function(x) unlist(x[[1]])))
  covars_log <- do.call(rbind, lapply(outlist, function(x) unlist(x[2])))
  colnames(covars_log) <-  c("feature", "covariates", "model.rsq", "model.pvalue")
  colnames(covar_adjusted) <- colnames(assay(D))
  rownames(covar_adjusted)<- rownames(assay(D))
  assay(D) <- covar_adjusted
  
  # add status information
  funargs <- maplet:::mti_funargs()
  D %<>% 
    maplet:::mti_generate_result(
      funargs = funargs,
      logtxt = sprintf('Adjusted for covariate effects with stepAIC'),
      output = covars_log
    )
  D
}
