#' Compute fold-change
#'
#' Add feature fold-changes to statistical results table.
#'
#' @param D \code{SummarizedExperiment} input.
#' @param stat_name Name of the statistical comparison.
#' @param fun Function with respect to which fold changes are computed (e.g. difference of the log means or median). Default:
#'    function(x){mean(x,na.rm=T)}.
#'
#' @return $results[[stat_name]]$output: fc column added to statistical table
#'
#' @examples
#' \dontrun{# add fold-changes to the result table of the statistical comparison called "comparison1", after correcting for variable "age"
#' ... %>%
#'  mt_post_fold_change(stat_name="comparison1") %>% ...}
#'
#'
#' @author JZ
#'
#' @export
mt_post_fold_change <- function(D,
                       stat_name,
                       fun = function(x){mean(x,na.rm=T)}){

  ## FOLDCHANGE FUNCTION (CONSIDER PREVIOUS LOG)
  if ((length(maplet::mtm_res_get_entries(D, c("pre","trans","log"))) != 1) &&
      (length(maplet::mtm_res_get_entries(D, c("load","flag","logged"))) != 1))
    stop("fold-changes can only be calculated for log-scale data")

  ## stat
  if(missing(stat_name))
    stop("stat_name must be given")
  ## FIND ENTRY
  stat_id <- metadata(D)$results %>%
    purrr::map_lgl(~"stats" %in% .x$fun && .x$output$name == stat_name) %>%
    which()
  if(length(stat_id) == 0)
    stop("stat element with name ", stat_name, " does not exist")
  if(length(stat_id)  > 1)
    stop("there are multiple stat elements with name ", stat_name)

  ## stop if the results are coming from an ANOVA, determine from first model
  if (class( metadata(D)$results[[ stat_id ]] %>%
             .$output %>%
             ## get first model in list (doesn't matter which)
             .$lstobj %>%
             .[[1]])[[1]] == "anova") {
    stop("Cannot add fold change for ANOVA analysis.")
  }

  ## GET OUTCOME (works for both lm and lmer)
  formula <- metadata(D)$results[[ stat_id ]]$output$formula
  terms <- all.vars(as.formula(formula))
  outcome <- terms[1]

  # find any random effects... they cannot be corrected for here
  # parse out by finding any variable inside (x | y)
  re_terms <- stringr::str_match(formula, "\\((.*?)\\|(.*?)\\)")[-1] %>% trimws()
  terms <- setdiff(terms, re_terms)


  ## CORRECT FOR CONFOUNDER
  if (length(terms)>1) {
    # construct new formula
    conf_form <- as.formula(sprintf("~%s", paste0(terms[2:length(terms)], collapse = "+")))
    maplet:::mti_logstatus(glue::glue("correcting for {as.character(conf_form)}"))
    D1 <- maplet:::mti_correctConfounder(D, conf_form) # store in different object to not change original SummarizedExperiment
    # get data-frame from corrected one
  } else {
    # no confounding correction
    D1 <- D
  }


  ## EXTRACT DATA
  d_fc <- cbind(colData(D)[,outcome,drop=F],t(assay(D1))) %>% as.data.frame() # only required reference to D1, since this uses the assay
  keep <- metadata(D)$results[[ stat_id ]]$output$samples.used
  d_fc <- d_fc[keep,,drop=F]
  d_fc[[outcome]] <- as.factor(d_fc[[outcome]])
  d_fc[[outcome]] <- droplevels(d_fc[[outcome]])

  ## CHECK TYPE OF OUTCOME
  if(!(class(d_fc[[ outcome ]]) %in% c("factor", "character")))
    stop(glue::glue("Fold-changes are only meaningful for factor/character, but {outcome} is a {class(d_fc[[outcome]])}"))
  if(length(levels(as.factor((d_fc[[ outcome ]])))) != 2)
    stop(glue::glue("Fold-changes are only meaningful for 2 groups, but {outcome} has {length(levels(as.factor((d_fc[[ outcome ]]))))}"))

  ## GET LEVELS IN ORDER OF FACTOR LEVELS
  outcome_q <- rlang::sym(outcome)
  if( "factor" %in% class(d_fc[[outcome]])){
    levels <- levels(d_fc[[outcome]])
  }else{
    levels <- unique(d_fc[[outcome]])
  }
  levels_1  <- rlang::sym(levels[1])
  levels_2  <- rlang::sym(levels[2])

  ## FOLDCHANGE FUNCTION (CONSIDER PREVIOUS LOG)
  d_fc <- d_fc %>%
    tidyr::gather(var, value, one_of(rownames(D))) %>%
    dplyr::group_by(var,!!outcome_q) %>%
    dplyr::summarise(value = fun(value)) %>%
    dplyr::ungroup() %>%
    tidyr::spread(!!outcome_q, value) %>%
    dplyr::mutate(fc = !!levels_1 - !!levels_2) %>%
    dplyr::select(var, fc)

  ## ADD TO RESULTS
  metadata(D)$results[[stat_id]]$output$table %<>%
    dplyr::left_join(d_fc, by = "var")

  ## make sure fold change has the same sign as statistic
  ## this is a debug solution, it should come out properly from the code above... but this fixes the bug for now
  metadata(D)$results[[stat_id]]$output$table$fc <-
    metadata(D)$results[[stat_id]]$output$table$fc %>% abs() * sign(metadata(D)$results[[stat_id]]$output$table$statistic)


  ## add status information & plot
  funargs <- maplet:::mti_funargs()
  D %<>% 
    maplet:::mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("Calculated foldchanges for %s", stat_name)
    )
  ## RETURN
  D
}

