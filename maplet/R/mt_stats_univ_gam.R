#' Generalized additive models using mgcv package
#'
#' Only supports standard linear models of the form: outcome ~ feature + [covariates].
#'
#' @param D \code{SummarizedExperiment} input.
#' @param outcome Outcome to be used for response.
#' @param outcome_type Outcome type - numeric/binary/ordinal. Default: numeric
#' @param rval For ordinal outcomes. Default: NULL.
#' @param conf_formula Confounders to be corrected for.
#' @param int_w_analyte Name of covariate that interacts with metabolite. Default: NULL.
#' @param stat_name Name under which this comparison will be stored, must be unique to all other statistical results.
#' @param samp_filter Term which samples to filter to first (e.g. used if the data contains >2 groups but the user wants
#'    to run a two-group comparison).
#'
#' @return $results$output: Statistics object.
#'
#' @import mgcv
#'
#' @examples
#' \donttest{# run lm with no confounders, "Group" as outcome
#' # filter to groups "Li_2" and "Li_5"
#' # name the comparison "Li's"
#' ... %>%
#'  mt_stats_univ_gam(
#'    conf_formula      = ~ Group,
#'    samp_filter = (Group %in% c("Li_2","Li_5")),
#'    stat_name         = "Li's"
#'  ) %>% ...
#'  }
#'
#' @author RB
#'
#' @export
mt_stats_univ_gam <- function(D,
                              outcome,
                              outcome_type='numeric',
                              rval=NULL,
                              conf_formula,
                              int_w_analyte = NULL,
                              stat_name,
                              samp_filter) {

  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))

  # make sure name does not exist yet
  if (stat_name %in% unlist(maplet::mtm_res_get_entries(D, "stats") %>% purrr::map("output") %>% purrr::map("name"))) stop(sprintf("stat element with name '%s' already exists",stat_name))

  # merge data with sample info
  Ds <- D %>% maplet:::mti_format_se_samplewise() # NOTE: No explosion of dataset size, no gather() - 8/17/20, JK

  ## FILTER SAMPLES
  if(!missing(samp_filter)) {

    filter_q <- dplyr::enquo(samp_filter)
    num_samp <- ncol(D)
    samples.used <- maplet:::mti_filter_samples(Ds, filter_q, num_samp)
    Ds <- Ds[samples.used,]

  } else {
    samples.used = rep(T, ncol(D))
  }

  # metabolites in data
  mets <- D %>% assay () %>% rownames ()

  # loop over metabolites
  univ_stats <- lapply(mets, function(x) {
    # formula for this metabolite
    # with interaction term?
    if (is.null(int_w_analyte) == F) {
      this_formula <-
        as.formula(glue('{outcome} ~ {x} + {conf_formula} + {int_w_analyte}*{x}'))
      # with covariates?
    } else {
      this_formula <-
        as.formula(glue('{outcome} ~ {x} + {conf_formula}'))
      # only metabolite?
    }
    # univariate analysis of numeric outcomes
    if (outcome_type == 'numeric') {
      # turn this outcome variable into numeric
      Ds %<>% mutate (!!sym(outcome) := as.numeric(as.matrix(!!sym(outcome))))
      # gaussian spline regression
      this_fit <- mgcv::gam(this_formula, family="gaussian", data=Ds)
    } else if (outcome_type == 'twofactor') {
      # turn the outcome variable into factor
      Ds %<>% mutate(!!sym(outcome) := as.factor(as.matrix(!!sym(outcome))))
      # binomial spline regression
      this_fit <- mgcv::gam(this_formula, family="binomial", data=Ds)
    } else if (outcome_type == 'ordinal') {

      # turn the outcome variable into factor
      Ds %<>% mutate(!!sym(outcome) := as.numeric(as.matrix(!!sym(outcome))))
      # ordinal spline regression
      this_fit <- mgcv::gam(as.formula(this_formula), family=ocat(R=rval), data=Ds)
    }

    # formating output accordingly
    this_res <- this_fit %>% anova() %>% .$p.table
    # empty matrix
    tmp <- tmp_names <- NULL
    # loop over rows
    for(i in c(2:nrow(this_res))){
      tmp <- c(tmp, as.matrix(this_res[i, ]))
      tmp_names <- c(tmp_names, paste0(rownames(this_res)[i], sep = '_', c('estimate', 'std_error', 'statistic', 'p_value')))
    }
    # remove specific metabolite name as that will be an issue when binding all results
    tmp_names <- sub('*..+\\:', 'analyte:', tmp_names)
    names(tmp) <- tmp_names
    # remove specific metabolite name as that will be an issue when binding all results
    names(tmp)[1:4] <- c('estimate', 'std_error','statistic', 'p_value')
    this_res <- tmp %>% data.matrix() %>% t() %>% data.frame()
    # results summary
    tmp <- this_fit %>% anova()
    if(!is.null(tmp$s.table)){
      tmp <- tmp$s.table %>% data.frame() %>% .[1, ]
      names(tmp) <- c('edf', 'ref_df', 'f_value', 'sp_value')
      this_res <- bind_cols(this_res, tmp)
    }

    # format results
    this_res %<>% mutate(
      analyte = x,
      outcome = outcome,
      covariates = conf_formula)
    # order output columns
    this_res %<>% select(analyte, outcome, everything())
    return(this_res)
  })%>% # create data from of results
    do.call(rbind, .) %>% data.frame()

  ## rename stat table columns to be maplet compatible
  univ_stats %<>% dplyr::rename(dplyr::all_of(c(var = "analyte", p.value = "p_value")))

  ## add status information & results
  funargs <- mti_funargs()
  D %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("GAM, %s", conf_formula),
      output = list(
        table = univ_stats,
        formula = sprintf('%s ~ feat + %s + interaction of %s', outcome, conf_formula, int_w_analyte),
        name    = stat_name,
        samples.used = samples.used,
        outcome = outcome
      )
    )

  ## return
  D

}
