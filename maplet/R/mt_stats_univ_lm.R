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
#' @param all_terms Whether to return all terms in formula in the statistical results table. Default: FALSE.
#' @param n_cores Number of cores to use for mclapply. More than one core will not work on Windows platforms. Default: 1.
#' @param return_warnings Whether to add a Boolean column indicating whether a feature model
#'    returned a warning. Currently only checks for a failure to converge warnings. See
#'    \code{\link[lme4]{isSingular()}}. Default: FALSE.
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
                             all_terms = FALSE,
                             n_cores = 1,
                             return_warnings = FALSE) {

  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))

  # make sure name does not exist yet
  if (stat_name %in% unlist(maplet::mtm_res_get_entries(D, "stats") %>% purrr::map("output") %>% purrr::map("name"))) stop(sprintf("stat element with name '%s' already exists",stat_name))

  # merge data with sample info
  Ds <- D %>% mti_format_se_samplewise() # NOTE: No explosion of dataset size, no gather() - 6/2/20, JK

  # apply sample filter
  if(!missing(samp_filter)) {

    filter_q <- dplyr::enquo(samp_filter)
    num_samp <- ncol(D)
    samples.used <- mti_filter_samples(Ds, filter_q, num_samp)
    Ds <- Ds[samples.used,] %>% droplevels()

  } else {
    samples.used = rep(T, ncol(D))
  }

  ## SAVE OUTCOME VARIABLE ##
  # Extract the fist attribute of the formula to use as outcome variable
  # Ff outcome variable is interaction term (contains ':'), set is_interaction flag to TRUE.
  outvar       <- attr(stats::terms(formula, keep.order = T),"term.labels")[1]
  outvar_label <- outvar
  is_interaction <- FALSE
  if(stringr::str_detect(outvar, ":")){
    outvar <- strsplit(outvar, ":") %>% unlist()
    is_interaction <- TRUE
    mti_logmsg("calculating interaction term")
  }
  if(any(!(outvar %in% colnames(Ds))))
    stop(sprintf("column %s do not exist in data", stringr::str_c(outvar[ !(outvar %in% colnames(Ds))], collapse = ", ")))

  ## HANDLE FACTOR OUTCOMES ##
  # If outcome variable is a factor, set do_anova flag to TRUE if > 2 levels. Crash if > 30 levels.
  # Extract outvar_term (outcome variable + second factor level). Required for filtering in lines [ADD LINE NUMBERS].
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


  ## CHOOSE LM & TIDY FUNCTIONS ##
  # If random effect term is present in the formula (i.e. a term containing '|'), set has_random_eff
  # flag to TRUE. Select the versions of lm and tidy functions to use depending on whether or not
  # a random effect term is present.
  has_random_eff <- FALSE
  if(stringr::str_detect(as.character(formula[2]), "\\|")){
    mti_logmsg("random effect detected; using lmer")
    f_lm     <- lmerTest::lmer
    f_tidy   <- broom.mixed::tidy # this calls summary()
    has_random_eff <- TRUE
  }else{
    f_lm   <- lm
    f_tidy <- broom::tidy # this calls summary()
  }
  # Define wrapper function for extracting information for the models and converting to a table.
  # Uses selected tidy method from above, called at line [INSERT LINE NUMBER]
  f_tidy_tidy <- function(m, ...){
    if(is.null(m))
      return(tibble::tibble(term = outvar_term))
    if(return_warnings){
      f_tidy(m) %>%
        dplyr::mutate(formula = as.character(attr(m,'terms'))) %>%
        dplyr::select(term, formula, dplyr::everything()) %>%
        dplyr::mutate(singular_warning = attr(m, 'singular_warning_generated'))
    }else{
      f_tidy(m) %>%
        dplyr::mutate(formula = as.character(attr(m,'terms'))) %>%
        dplyr::select(term, formula, dplyr::everything())
    }
  }


  ## VALIDATE FORMULA ##
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


  ## FUNCTION TO CREATE MDOELS ##
  # For a feature m:
  #   add m to LHS of formula
  #   check whether m, outcome, and / or confounders are invariant
  #     if m or outcome invariant, return NULL for model
  #     if other confounders invariant, remove them from formula
  #   run the linear model
  do_lm <- function(m){
    # add m to LHS of formula
    form <- stats::update(formula, sprintf("%s~.",m))
    # extra formula terms and attach classes
    trms <- all.vars(formula) %>%
      purrr::discard(~stringr::str_detect(.x, ":")) %>%
      c(m)
    clss <- purrr::map_chr(trms, ~class(Ds[[.x]])) %>%
      stats::setNames(trms)

    ## subset to complete data
    d <- Ds %>%
      dplyr::select(dplyr::one_of(trms), !!rlang::sym(m)) %>%
      dplyr::filter(stats::complete.cases(.))

    ## check for invariant confounders
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
      # return NULL for model if feature or outcome is invariant
      if(m %in% names(conf_invar)){
        mti_logwarning(glue::glue("feature {m} invariant "))
        return(NULL)
      }
      if(any(outvar %in% names(conf_invar))){
        mti_logwarning(glue::glue("outcome {outvar} invariant for feature {m}"))
        return(NULL)
      }
      # remove invariant confounders from formula
      for(c in names(conf_invar)){
        mti_logwarning(glue::glue("confounder {c} invariant for feature {m}, removing from formula"))
        form <- stats::update.formula(form, stringr::str_c(". ~ . -", c))
      }
    }

    # calculate linear model
    mod <- f_lm(
      data    = Ds,
      formula = form
    )
    terms <- terms(mod)

    # perform ANOVA if multiple factor levels
    if(do_anova)
      mod <- stats::anova(mod)
    # attach linear model terms as attribute to the variable
    # (only way to make this compatible for both lm and anova, because they are very different data structures)
    attr(mod, 'terms') <- terms
    if(return_warnings == TRUE) attr(mod, 'singular_warning_generated') <- lme4::isSingular(mod)

    # return
    mod
  } # END OF do_lm() FUNCTION DEFINITION

  ## RUN TESTS FOR ALL FEATURES ##
  models <- parallel::mclapply(rownames(D), do_lm, mc.cores = n_cores) %>%
    stats::setNames(rownames(D))

  ## CONVERT LIST OF MODELS TO STATS TABLE ##
  if(all_terms){
    # get statistical table for all terms
    tab <- purrr::map_dfr(models, f_tidy_tidy, conf.int = T, .id = "var")
    form_terms <- tab$term %>% unique() %>% .[.!= outvar_term] %>% c(outvar, .)
    exclude_cols <- c("var", "term", "formula", "effect", "group")
    statres_cols <- dplyr::setdiff(colnames(tab), exclude_cols)
    ordered_cols <-  as.vector(outer(statres_cols, form_terms, paste, sep="_"))

    # format and order table
    term_key <- outvar
    names(term_key) <- outvar_term
    tab %<>% dplyr::select(-dplyr::any_of(c("effect", "group"))) %>%
      dplyr::mutate(term=dplyr::recode(term, !!!term_key)) %>%
      tidyr::pivot_wider(names_from = term, values_from = dplyr::all_of(statres_cols)) %>%
      dplyr::mutate(term = outvar) %>%
      dplyr::select(var, term, formula, dplyr::all_of(ordered_cols)) %>%
      dplyr::select_if(~!all(is.na(.)))

    # rename the outcome columns
    old_names <- tab %>% dplyr::select(ends_with(outvar)) %>% colnames()
    new_names <- gsub(paste0("_", outvar), "", old_names)
    tab %<>% dplyr::rename_with(~ new_names[which(old_names == .x)], .cols = dplyr::all_of(old_names))

  }else{
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
    # skip model if NULL
    # if we don't do this, a very strange corruption of the models list can occur.
    if (is.null(models[[i]])) {
      next  # Skip to the next iteration if NULL
    }
    # set everything to NULL that is not needed
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
  funargs <- mti_funargs()
  D %<>%
    mti_generate_result(
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
