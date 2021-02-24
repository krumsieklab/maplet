#' Univariate GLMs
#'
#' Computes univariate GLM for each feature.
#'
#' @description
#' \enumerate{
#'   \item Will treat the first term of the formula as outcome.
#'   \item If outcome has >2 factors, will perform ANOVA.
#'   \item If random effect term is contained, will used lmer function (can also be used for paired analysis, e.g. before/after).
#' }
#'
#' @param D \code{SummarizedExperiment} input.
#' @param formula Left-hand side of formula to be put into glm function.
#' @param stat_name Name under which this comparison will be stored, must be unique to all other statistical results.
#' @param samp_filter Term which samples to filter to first (e.g. used if the data contains >2 groups but the user wants to
#'    run a two-group comparison).
#' @param n_cores Number of cores to use for mclapply. More than one core will not work on Windows platforms. Default: 1.
#'
#' @return $result$output: Statistics object.
#'
#' @examples
#' \donttest{# run lm with no confounders, "Group" as outcome
#' # filter to groups "Li_2" and "Li_5"
#' # name the comparison "Li's"
#' ... %>%
#'  mt_stats_univ_lm(
#'    formula      = ~ Group,
#'    samp_filter = (Group %in% c("Li_2","Li_5")),
#'    stat_name         = "Li's"
#'  ) %>% ...
#'  }
#'
#' @author JK, JZ
#'
#' @export
mt_stats_univ_lm <- function(D,
                             formula,
                             stat_name,
                             samp_filter,
                             n_cores = 1) {

  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))

  # make sure name does not exist yet
  if (stat_name %in% unlist(maplet::mtm_res_get_entries(D, "stats") %>% purrr::map("output") %>% purrr::map("stat_name"))) stop(sprintf("stat element with name '%s' already exists",stat_name))

  # merge data with sample info
  Ds <- D %>% maplet:::mti_format_se_samplewise() # NOTE: No explosion of dataset size, no gather() - 6/2/20, JK

  ## FILTER SAMPLES
  if(!missing(samp_filter)) {

    filter_q <- dplyr::enquo(samp_filter)
    num_samp <- ncol(D)
    samples.used <- maplet:::mti_filter_samples(Ds, filter_q, num_samp)
    Ds <- Ds[samples.used,] %>% droplevels()

  } else {
    samples.used = rep(T, ncol(D))
  }

  ## save outcome variable
  outvar       <- attr(stats::terms(formula, keep.order = T),"term.labels")[1]
  outvar_label <- outvar
  is_interaction <- FALSE
  if(stringr::str_detect(outvar, ":")){
    outvar <- strsplit(outvar, ":") %>% unlist()
    is_interaction <- TRUE
    maplet:::mti_logmsg("calculating interaction term")
  }
  if(any(!(outvar %in% colnames(Ds))))
    stop(sprintf("column %s do not exist in data", stringr::str_c(outvar[ !(outvar %in% colnames(Ds))], collapse = ", ")))

  # handle factor outcomes
  do_anova    <- FALSE
  outvar_term <- outvar
  for(o in seq_along(outvar)){
    v = outvar[o]
    if (is.character(Ds[[v]]))
      Ds[[ v ]] = as.factor(Ds[[ v ]])
    if (is.factor(Ds[[v]])) {
      # check that there are exactly two levels
      if (length(levels( Ds[[ v]] ))!=2){
        # stop(sprintf("factor outcomes must have exactly two levels, '%s' has %d", o, length(levels(v))))
        do_anova <- TRUE
        # sanity check: if there are more than 30 levels: crash
        if (length(levels( Ds[[ v]] ))>30) {
          stop(sprintf("More than 30 factors in outcome '%s'... forgot to convert string to numeric?",v))
        }
      }
      # remember the level name of the second (will be deleted later on)
      outvar_term[o] <- stringr::str_c(outvar[o], levels(Ds[[ v ]])[2])
    }
  }
  if(is_interaction)
    outvar_term <- stringr::str_c(outvar_term, collapse = ":")


  ## choose lm functions
  has_random_eff <- FALSE
  if(stringr::str_detect(as.character(formula[2]), "\\|")){
    maplet:::mti_logmsg("random effect detected; using lmer")
    f_lm     <- lmerTest::lmer
    f_tidy   <- broom.mixed::tidy # this calls summary()
    has_random_eff <- TRUE
  }else{
    f_lm   <- lm
    f_tidy <- broom::tidy # this calls summary()
  }
  f_tidy_tidy <- function(m, ...){
    if(is.null(m))
      return(tibble::tibble(term = outvar_term))
    f_tidy(m) %>%
      dplyr::mutate(formula = as.character(attr(m,'terms'))) %>%
      dplyr::select(term, formula, dplyr::everything())
  }


  ## validate formula
  if (!is.null(formula.tools::lhs(formula))) stop("Left-hand side of formula must be empty")
  ## will crash sort of meaningfully if variables don't exist
  if(has_random_eff){
    mm <- lme4::lFormula(stats::update.formula(formula, stringr::str_c(rownames(D)[[1]], "~.")), Ds)
  }else{
    mm <- stats::model.matrix(formula,Ds)
  }

  ## check if anova is necessary for interaction effects
  ## (happens if 2 factors with 2 levels each where individual variables are not included)
  if(is_interaction){
    interactions <- purrr::map(outvar, function(.x){
      if(is.factor(Ds[[.x]]))
        stringr::str_c(.x, unique(Ds[[.x]]))
      else
        .x
    })%>%
      expand.grid() %>%
      apply(MARGIN = 1, stringr::str_c, collapse = ":")
    if(length(dplyr::intersect(colnames(mm), interactions)) > 1){
      do_anova <- TRUE
      stop("for interactions effects between factors, individual terms should be included")
    }
  }


  ## function to create models
  ## does checks to minimise errors
  do_lm <- function(m){
    ## run glm with updated formula
    form <- stats::update(formula, sprintf("%s~.",m))
    ## check for constant confounders
    trms <- all.vars(formula) %>%
      purrr::discard(~stringr::str_detect(.x, ":")) %>%
      c(m)
    clss <- purrr::map_chr(trms, ~class(Ds[[.x]])) %>%
      stats::setNames(trms)
    ## subset to complete data
    d <- Ds %>%
      dplyr::select(dplyr::one_of(trms), !!rlang::sym(m)) %>%
      dplyr::filter(stats::complete.cases(.))
    ## check for invariant vonfounders
    conf_invar_num <- clss %>%
      purrr::keep(~.x %in% c("integer", "numeric")) %>%
      purrr::imap(~stats::var(d[[.y]])) %>%
      purrr::keep(~.x == 0)
    conf_invar_fct <- clss %>%
      purrr::discard(~.x %in% c("integer", "numeric")) %>%
      purrr::imap(~length(unique(d[[.y]]))) %>%
      purrr::keep(~.x == 1)
    conf_invar <- c(conf_invar_num, conf_invar_fct)
    if(length(conf_invar) > 0){
      ## terminate if feature or outcome are invariant
      if(m %in% names(conf_invar)){
        maplet:::mti_logwarning(glue::glue("feature {m} invariant "))
        return(NULL)
      }
      if(any(outvar %in% names(conf_invar))){
        maplet:::mti_logwarning(glue::glue("outcome {outvar} invariant for feature {m}"))
        return(NULL)
      }
      for(c in names(conf_invar)){
        maplet:::mti_logwarning(glue::glue("confounder {c} invariant for feature {m}, removing from formula"))
        form <- stats::update.formula(form, stringr::str_c(". ~ . -", c))
      }
    }

    ## CALCLUATE ACTUAL MODEL
    mod <- f_lm(
      data    = Ds,
      formula = form
    )
    terms <- terms(mod)
    ## DO ANOVA IF MULTIPLE FACTOR LEVELS
    if(do_anova)
      mod <- stats::anova(mod)
    ## attach linear model terms as attribute to the variable
    ## (only way to make this compatible for both lm and anova, because they are very different data structures)
    attr(mod, 'terms') <- terms
    ## RETURN
    mod
  }

  ## run tests for all features
  models <- parallel::mclapply(rownames(D), do_lm, mc.cores = n_cores) %>%
    stats::setNames(rownames(D))

  # broom it up, subselect to term, rename term
  if (do_anova) {
    tab <- purrr::map_dfr(models, f_tidy_tidy, conf.int = T, .id = "var") %>%
      dplyr::filter(term == outvar) %>%
      dplyr::mutate(term =  outvar_label)
  } else {
    tab <- purrr::map_dfr(models, f_tidy_tidy, conf.int = T, .id = "var") %>%
      dplyr::filter(term == outvar_term) %>%
      dplyr::mutate(term =  outvar_label)
  }

  ## tidy up a bit more
  tab <- tab %>%
    dplyr::select(-dplyr::matches("^(effect|group)$"))

  ## construct output groups variable
  if (is.factor(Ds[[outvar]])) {
    outgroups <- levels(Ds[[outvar]])
  } else {
    outgroups <- NULL
  }

  ## remove environments from all models (blows up when saving to file)
  for (i in 1:length(models)) {
    attr(models[[i]], "terms") <- NULL
    if (!("lmerModLmerTest" %in% class(models[[i]]))) {
      # regular LM
      attr(models[[i]], "terms") <- NULL
      attr(models[[i]]$model,"terms") <- NULL
      models[[i]]$terms <- NULL
    } else {
      # LMER
      attr(models[[i]]@frame, "terms") <- NULL
      attr(models[[i]]@frame, "formula") <- NULL
    }
  }
  # make sure that NAs in the outcome are set to FALSE in the list of used samples
  samples.used[is.na(Ds[, outvar])] <- F

  ## add status information & results
  funargs <- maplet:::mti_funargs()
  metadata(D)$results %<>%
    maplet:::mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("univariate lm, %s", as.character(formula)),
      output = list(
        table   = tab,
        formula = dput(formula) %>% as.character(),
        name    = stat_name,
        lstobj  = models,
        groups = outgroups,
        samples.used = samples.used,
        outcome = outvar
      )
    )

  ## return
  D
}
