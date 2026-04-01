#' Univariate Concordance Test with LOD-Censored Missing Values (rox)
#'
#' Computes a rank-order concordance test for each feature using the rox package,
#' which handles missing values due to limit of detection (LOD) without imputation.
#' Missing values are treated as left-censored observations: they are ranked below all
#' observed values while preserving their partial ordering information.
#'
#' @description
#' \enumerate{
#'   \item Treats the first term of the formula as outcome (same convention as \code{mt_stats_univ_lm}).
#'   \item If only an outcome is specified (no confounders), uses \code{rox::rox()} for a simple
#'         concordance test.
#'   \item If confounders are present, uses \code{rox::rox_mv()} with concordance regression
#'         (\code{concreg}) or Cox approximation (\code{coxph}) to obtain confounder-adjusted results.
#'   \item The \code{estimate} column contains the Somers' D statistic (2*concordance - 1), ranging
#'         from -1 to 1, providing a signed effect size compatible with downstream pathway direction plots.
#' }
#'
#' @param D \code{SummarizedExperiment} input.
#' @param formula Right-hand side formula specifying outcome and optional confounders (e.g. \code{~ Group}
#'   or \code{~ Group + Age + Sex}). First term is the outcome, remaining terms are confounders.
#'   Left-hand side must be empty.
#' @param stat_name Name under which this comparison will be stored, must be unique to all other statistical results.
#' @param samp_filter Term which samples to filter to first (e.g. used if the data contains >2 groups but the user wants to
#'    run a two-group comparison).
#' @param switching Logical. If TRUE (default), rox will fall back to complete-case analysis when no
#'    evidence of an LOD effect is found. See \code{rox::rox()} for details.
#' @param mv_func For multivariable models (when confounders are present), the fitting method to use.
#'    One of \code{"concreg"} (concordance regression, default) or \code{"coxph"} (Cox partial likelihood,
#'    faster for large datasets). Only relevant when confounders are in the formula.
#' @param n_cores Number of cores to use for mclapply. More than one core will not work on Windows platforms. Default: 1.
#'
#' @return $result$output: Statistics object with table containing columns: var, term, estimate,
#'   concordance, se, statistic, p.value, w, formula.
#'
#' @examples
#' \donttest{# run rox with no confounders, "Group" as outcome
#' # filter to groups "AD" and "Control"
#' ... %>%
#'  mt_stats_univ_rox(
#'    formula   = ~ Group,
#'    samp_filter = (Group %in% c("AD","Control")),
#'    stat_name = "AD_vs_Control_rox"
#'  ) %>% ...
#'
#' # with confounders
#' ... %>%
#'  mt_stats_univ_rox(
#'    formula   = ~ Group + Age + Sex,
#'    stat_name = "AD_vs_Control_rox_adj"
#'  ) %>% ...
#' }
#'
#' @author JK (wrapper)
#'
#' @export
mt_stats_univ_rox <- function(D,
                              formula,
                              stat_name,
                              samp_filter,
                              switching = TRUE,
                              mv_func = "concreg",
                              n_cores = 1) {

  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  if (!requireNamespace("rox", quietly = TRUE))
    stop("Package 'rox' is required. Install with: devtools::install_github('krumsieklab/rox', subdir='rox')")
  if (length(formula) == 3)
    stop("Left-hand side of formula must be empty")
  if (!(mv_func %in% c("concreg", "coxph")))
    stop("mv_func must be one of 'concreg' or 'coxph'")

  # make sure stat_name does not exist yet
  if (stat_name %in% unlist(mtm_res_get_entries(D, "stats") %>%
                            purrr::map("output") %>% purrr::map("name")))
    stop(sprintf("stat element with name '%s' already exists", stat_name))

  # merge data with sample info
  Ds <- D %>% mti_format_se_samplewise()

  # apply sample filter
  if (!missing(samp_filter)) {
    filter_q <- dplyr::enquo(samp_filter)
    num_samp <- ncol(D)
    samples.used <- mti_filter_samples(Ds, filter_q, num_samp)
    Ds <- Ds[samples.used, ] %>% droplevels()
  } else {
    samples.used <- rep(TRUE, ncol(D))
  }

  # extract outcome variable (first term) and confounders (remaining terms)
  all_terms <- attr(stats::terms(formula, keep.order = TRUE), "term.labels")
  if (length(all_terms) == 0) stop("Formula must contain at least one term (the outcome variable)")
  outvar <- all_terms[1]
  confounders <- if (length(all_terms) > 1) all_terms[-1] else NULL
  has_confounders <- !is.null(confounders)

  # validate that all formula variables exist in data
  all_formula_vars <- c(outvar, confounders)
  missing_vars <- all_formula_vars[!(all_formula_vars %in% colnames(Ds))]
  if (length(missing_vars) > 0)
    stop(sprintf("Column(s) %s do not exist in data", paste(missing_vars, collapse = ", ")))

  # convert outcome: ensure it is numeric for rox
  outcome_vec <- Ds[[outvar]]
  if (is.character(outcome_vec)) outcome_vec <- as.factor(outcome_vec)
  if (is.factor(outcome_vec)) {
    outgroups <- levels(outcome_vec)
    outcome_vec <- as.numeric(outcome_vec)
  } else {
    outgroups <- NULL
  }

  # build confounder data.frame for rox_mv if needed
  if (has_confounders) {
    conf_df <- Ds[, confounders, drop = FALSE]
  }

  # per-feature rox test
  do_rox <- function(m) {
    xh <- Ds[[m]]
    y <- outcome_vec

    # rox requires y to have no NAs; xh NAs are the whole point
    valid <- !is.na(y)
    if (has_confounders) {
      valid <- valid & stats::complete.cases(conf_df)
    }
    xh <- xh[valid]
    y  <- y[valid]

    if (all(is.na(xh)) || (!any(is.na(xh)) && stats::var(xh, na.rm = TRUE) == 0)) {
      return(NULL)
    }

    result <- tryCatch({
      if (has_confounders) {
        ys_df <- data.frame(y = y, conf_df[valid, , drop = FALSE])
        fit <- rox::rox_mv(xh, ys_df, switching = switching, func = mv_func)
        y_row <- fit$stats_model["y", , drop = TRUE]
        # concreg's p from 2*pnorm(-z) is wrong for negative z; compute two-sided
        z_val <- unname(y_row["z"])
        list(
          estimate    = unname(y_row["B"]),
          concordance = unname(fit$stats_overall["d"]),
          se          = unname(y_row["se"]),
          statistic   = z_val,
          p.value     = 2 * pnorm(-abs(z_val)),
          w           = fit$w,
          fit         = fit
        )
      } else {
        fit <- rox::rox(xh, y, switching = switching)
        list(
          estimate    = unname(2 * fit$stats["d"] - 1),
          concordance = unname(fit$stats["d"]),
          se          = unname(fit$stats["se"]),
          statistic   = unname(fit$stats["z"]),
          p.value     = unname(fit$stats["p"]),
          w           = fit$w,
          fit         = fit
        )
      }
    }, error = function(e) {
      mti_logwarning(glue::glue("rox failed for feature {m}: {e$message}"))
      NULL
    })

    result
  }

  # run tests for all features
  rox_results <- parallel::mclapply(rownames(D), do_rox, mc.cores = n_cores) %>%
    stats::setNames(rownames(D))

  # assemble results table
  tab <- purrr::imap_dfr(rox_results, function(res, varname) {
    if (is.null(res)) {
      tibble::tibble(
        var         = varname,
        term        = outvar,
        estimate    = NA_real_,
        concordance = NA_real_,
        se          = NA_real_,
        statistic   = NA_real_,
        p.value     = NA_real_,
        w           = NA_real_,
        formula     = as.character(formula)[2]
      )
    } else {
      tibble::tibble(
        var         = varname,
        term        = outvar,
        estimate    = res$estimate,
        concordance = res$concordance,
        se          = res$se,
        statistic   = res$statistic,
        p.value     = res$p.value,
        w           = res$w,
        formula     = as.character(formula)[2]
      )
    }
  })

  # strip heavy objects from stored results to keep SE size manageable
  lstobj <- purrr::map(rox_results, function(res) {
    if (is.null(res)) return(NULL)
    fit <- res$fit
    fit$dfit$call <- NULL
    if ("dfit_overall" %in% names(fit)) fit$dfit_overall$call <- NULL
    fit
  })

  # mark samples with NA outcome as not used
  samples.used[is.na(Ds[[outvar]])] <- FALSE

  # construct formula string for logging
  formula_str <- as.character(formula)[2]
  if (has_confounders) {
    logtxt <- sprintf("univariate rox (confounder-adjusted via %s), %s", mv_func, formula_str)
  } else {
    logtxt <- sprintf("univariate rox, %s", formula_str)
  }

  # add status information & results
  funargs <- mti_funargs()
  D %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt  = logtxt,
      output  = list(
        table        = tab,
        formula      = formula_str,
        name         = stat_name,
        lstobj       = lstobj,
        groups       = outgroups,
        samples.used = samples.used,
        outcome      = outvar
      )
    )

  # return
  D
}
