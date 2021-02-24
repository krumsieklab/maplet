#' Confounding correction
#'
#' Each feature is corrected by lm with the confounding variables and residuals are returned in lieu of uncorrected
#' features expressions.
#'
#' @param D \code{SummarizedExperiment} object.
#' @param formula Formula including column names from sample annotation, i.e, ~batch+age.
#' @param strata Strata classes for stratified confounding correction, should be a column from colData(D). Default: NULL.
#' @param scale_data Should data be scaled after correction? Default: F.
#'
#' @return assay: Corrected data.
#' @return $results$output: The lm.fit pvalue per features to show confounding affect on them.
#'
#' @examples
#'  \dontrun{# not run
#'  #... %>% mt_pre_confounding_correction( formula = ~batch + age, strata = "RUN_DAY" )
#'  }
#'
#' @author MB
#'
#' @export
mt_pre_confounding_correction <- function(D, formula, strata = NULL, scale_data = F) {

  # validate and extract arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  X = assay(D)

  # covariates of samples
  dfc = colData(D)
  frm = stats::as.formula(paste0(c("y",as.character(formula)), collapse = ""))
  if(!is.null(strata)) strata = dfc[,strata]
  else strata = rep(1, nrow(dfc))

  # residual function which supports NA
  fresid <- function(fit){
    resids = fit$residuals
    naa = fit$na.action
    if(is.null(naa)) return(resids)

    resh <- rep(NA, length(resids)+length(naa))
    names(resh) = seq(length(resh))
    names(resh)[naa] = names(naa)
    names(resh)[-naa] = names(resids)
    resh[-naa] = resids
    resh
  }

  # function that returns fit p value
  fpv<- function(fit){
    f <- summary(fit)$fstatistic
    p <- stats::pf(f[1],f[2],f[3],lower.tail=F)
    log(unname(p))
  }

  # run correction
  lm.fits <- apply(X %>% t,2, function(y){
    pocs = lapply(unique(sort(strata)) %>% {names(.)=. ; .},
                  function(i) stats::lm(frm, data.frame(y=y, dfc)[strata==i,]))

    # indices to go back original indices
    rh = unlist(lapply(pocs, fresid))[order(order(strata))]
    # model pvalues
    p.lms = sapply(pocs, fpv)
    list(rh=rh, p.lms=p.lms)
  }) #residuals(%>% scale_data %>% t

  # # overall effect of covariates to each variable as lm pvalue
  # p.lms  <- sapply(lm.fits, function(fit){
  #   f <- summary(fit)$fstatistic
  #   p <- stats::pf(f[1],f[2],f[3],lower.tail=F)
  #   log(unname(p))
  # })

 Xh = do.call(rbind, lapply(lm.fits, `[[`,"rh"))
 p.lms = sapply(lm.fits, `[[`,"p.lms")

 rm(lm.fits)

 if(!identical(dim(X), dim(Xh)))
   stop("missing values in confounding variables not supported")

  # make sure colnames and rownames are kept
  colnames(Xh) <- colnames(X)
  rownames(Xh) <- rownames(X)

  # replace original assay with corrected verison
  assay(D) = Xh

  # add status information
  funargs <- mti_funargs()
  metadata(D)$results %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("covariates adjustment: '%s'", as.character(formula)),
      output = list(
        pvals   = p.lms,
        formula = formula
      )
    )

  # return
  D

}

# maplet
#
# covariate adjustment of features
#
# last update: 2019-03-29
# authors: MB
#
#! returns corrected data scaled
#
# consider missing values
# has to be colnames

