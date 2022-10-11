#' Fast linear models using MatrixEQTL package
#'
#' Substantially faster than regular lm as implemented in mt_stats_univ_lm.
#' Only supports standard linear models of the form: outcome ~ feature + [covariates].
#'
#' @param D \code{SummarizedExperiment} input.
#' @param formula Right-hand side of formula to be put into glm function.
#' @param stat_name Name under which this comparison will be stored, must be unique to all other statistical results.
#' @param samp_filter Term which samples to filter to first (e.g. used if the data contains >2 groups but the user wants
#'    to run a two-group comparison).
#'
#' @return $results$output: Statistics object.
#'
#' @import MatrixEQTL
#'
#' @examples
#' \donttest{# run lm with no confounders, "Group" as outcome
#' # filter to groups "Li_2" and "Li_5"
#' # name the comparison "Li's"
#' ... %>%
#'  mt_stats_univ_lm_matrixeqtl(
#'    formula      = ~ Group,
#'    samp_filter = (Group %in% c("Li_2","Li_5")),
#'    stat_name         = "Li's"
#'  ) %>% ...
#'  }
#'
#' @author JK
#'
#' @export
mt_stats_univ_lm_matrixeqtl <- function(D, formula, stat_name, samp_filter) {

  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))

  # make sure name does not exist yet
  if (stat_name %in% unlist(maplet::mtm_res_get_entries(D, "stats") %>% purrr::map("output") %>% purrr::map("name"))) stop(sprintf("stat element with name '%s' already exists",stat_name))

  # merge data with sample info
  Ds <- D %>% mti_format_se_samplewise() # NOTE: No explosion of dataset size, no gather() - 8/17/20, JK

  ## FILTER SAMPLES
  if(!missing(samp_filter)) {

    filter_q <- dplyr::enquo(samp_filter)
    num_samp <- ncol(D)
    samples.used <- mti_filter_samples(Ds, filter_q, num_samp)
    Ds <- Ds[samples.used,]

  } else {
    samples.used = rep(T, ncol(D))
  }

  # validate formula (can't contain random or interaction effects)
  strform <- formula %>% format() 
  if (grepl("\\|", strform)) stop(sprintf("Formula seems to contain a random effect. Linear mixed models not supported for MatrixEQTL. %s", strform))
  # this could be implemented later:
  if (grepl("\\:", strform) || grepl("\\*", strform)) stop(sprintf("Formula seems to contain a interaction effect. Currently not supported for MatrixEQTL. %s", strform))

  # dissect formula
  vars <- attr(stats::terms(formula, keep.order = T),"factors") %>% rownames()
  outvar <- vars[1]
  covars <- vars[-1]

  # check that outcome is either binary or numerical
  outvec <- Ds[[outvar]]
  cl <- outvec %>% class()
  if (("character" %in% cl) || ("factor" %in% cl)) {
    if ((outvec %>% as.factor() %>% levels() %>% length()) != 2) {
      stop("If outcome is a factor, it must have exactly two levels")
    }
    # now convert to 0/1
    outvec %<>% as.factor() %>% as.numeric()
    # save groups
    outgroups <- Ds[[outvar]] %>% as.factor() %>% levels()
  } else {
    outgroups <- NULL
  }

  # build SliceData objects
  Y <- MatrixEQTL::SlicedData$new()
  Y$CreateFromMatrix(assay(D[,samples.used])) # leave in features X samples form
  X <- MatrixEQTL::SlicedData$new()
  X$CreateFromMatrix(matrix(outvec, nrow = 1))
  CO <- MatrixEQTL::SlicedData$new()
  CO$CreateFromMatrix(Ds[,covars] %>% t())

  # run MatrixEQTL
  # from http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/features.html:
  # modelLINEAR: expression = alpha + sum_k beta_k*covariate_k + gamma * genotype_additive
  # we want to model feature ~ outcome + covariates
  me <- MatrixEQTL::Matrix_eQTL_main(
    snps = X,
    gene = Y,
    cvrt = CO,
    pvOutputThreshold = 1, # output all p-values
    verbose = F
  )

  # restructure result table so it fits the standard form
  table <- me$all$eqtls %>%
    dplyr::select(-snps) %>%
    dplyr::rename(var=gene) %>%
    dplyr::rename(p.value=pvalue)
  # rearrange back to original order
  o <- match(rownames(D),table$var)
  stopifnot(!any(is.na(o))) # sanity check
  table <- table[o,]

  ## add status information & results
  funargs <- mti_funargs()
  D %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("MatrixEQTL lm, %s", strform),
      output = list(
        table = table,
        formula = paste0("~", strform),
        name    = stat_name,
        groups = outgroups,
        samples.used = samples.used,
        outcome = outvar
      )
    )

  ## return
  D

}
