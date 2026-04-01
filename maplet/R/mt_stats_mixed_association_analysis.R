#' Mixed Association Analysis
#' 
#' Performs association analysis using regular GLM/LM or mixed models depending on whether
#' random effects are specified. Supports numeric, binary, and ordinal outcomes.
#' 
#' @param D \code{SummarizedExperiment} input.
#' @param stat_name Name of the statistical result.
#' @param outcome Outcome to be used for response.
#' @param outcome_type Outcome type - one of numeric/binary/ordinal.
#' @param conf_formula Confounders to be corrected for. Default: NULL.
#' @param random_effects Random effects specification (e.g., "(1|subject_id)"). If NULL, uses regular GLM/LM. Default: NULL.
#' @param log_progress If TRUE, logs progress for each metabolite being processed. Default: FALSE.
#'
#' @return $results$output: List with 'name' (stat_name), 'table' (dataframe with columns: 'estimate', 'std.error', 'df', 'statistic', 'p.value', 'p.adj'), 'outcome', and 'samples.used'.
#'
#' @examples
#' \dontrun{# ADD EXAMPLE}
#'
#'
#' @author RB
#' 
#' @export
mt_stats_mixed_association_analysis <- function(
  D,
  stat_name,
  outcome,
  outcome_type,
  conf_formula=NULL,
  random_effects=NULL,
  log_progress=FALSE
){
  
  # clean conf_formula if provided (remove leading '+' with optional whitespace)
  if(!is.null(conf_formula)){
    conf_formula <- trimws(gsub("^\\s*\\+\\s*", "", conf_formula))
  }
  
  # merge data with sample info
  Ds <- D %>% mti_format_se_samplewise()
  
  # convert outcome variable once based on type
  if(outcome_type=='binary'){
    Ds %<>% mutate(!!sym(outcome) := as.factor(as.matrix(!!sym(outcome))))
  } else if(outcome_type=='numeric'){
    Ds %<>% mutate(!!sym(outcome) := as.numeric(as.matrix(!!sym(outcome))))
  } else if(outcome_type=='ordinal'){
    Ds %<>% mutate(!!sym(outcome) := as.ordered(as.factor(as.matrix(!!sym(outcome)))))
  } else {
    stop(sprintf("Invalid outcome_type '%s'. Must be one of: 'numeric', 'binary', 'ordinal'.", outcome_type))
  }
  
  # metabolites in data
  mets <- D %>% assay() %>% rownames()
  
  # loop over metabolites - fit models and collect results
  results_list <- lapply(seq_along(mets), function(i){
    x <- mets[i]
    
    if(log_progress){
      mti_logstatus(sprintf("Processing metabolite %d of %d: %s", i, length(mets), x))
    }
    
    # wrap metabolite name in backticks to handle special characters
    x_safe <- paste0('`', x, '`')
    
    # formula for this metabolite
    if(is.null(conf_formula)==F){
      base_formula <- glue::glue('{outcome} ~ {x_safe} + {conf_formula}')
    } else{
      base_formula <- glue::glue('{outcome} ~ {x_safe}')
    }
    
    # add random effects if specified
    if(!is.null(random_effects)){
      this_formula <- as.formula(glue::glue('{base_formula} + {random_effects}'))
    } else {
      this_formula <- as.formula(base_formula)
    }
    
    formula_str <- paste(deparse(this_formula), collapse = " ")
    
    if(outcome_type=='binary'){
      if(!is.null(random_effects)){
        this_fit <- lme4::glmer(this_formula, data =Ds, family = binomial,
                                nAGQ = 2)
      } else {
        this_fit <- glm(this_formula, data =Ds, family = binomial)
      }
      this_res <- this_fit %>% summary() %>% coefficients() %>% data.frame()
      names(this_res) <- c('estimate', 'std.error', 'statistic', 'p.value')
      this_res$df <- NA
    } else if(outcome_type=='numeric'){
      if(!is.null(random_effects)){
        this_fit <- lme4::lmer(this_formula, data =Ds)
        this_res <- this_fit %>% lmerTest::as_lmerModLmerTest() %>% 
          summary() %>% coefficients() %>% data.frame()
        names(this_res) <- c('estimate', 'std.error', 'df', 'statistic', 'p.value')
      } else {
        this_fit <- lm(this_formula, data =Ds)
        this_res <- this_fit %>% summary() %>% coefficients() %>% data.frame()
        names(this_res) <- c('estimate', 'std.error', 'statistic', 'p.value')
        this_res$df <- NA
      }
    } else if(outcome_type=='ordinal'){
      if(!is.null(random_effects)){
        this_fit <- ordinal::clmm(this_formula, data =Ds)
      } else {
        this_fit <- ordinal::clm(this_formula, data =Ds)
      }
      this_res <- this_fit %>% summary() %>% coefficients() %>% data.frame()
      names(this_res) <- c('estimate', 'std.error', 'statistic', 'p.value')
      this_res$df <- NA
    }
    
    # extract the metabolite coefficient by name for ordinal (thresholds come first),
    # by row index for lm/glm (row 1 = intercept, row 2 = metabolite)
    if(outcome_type=='ordinal'){
      if(x %in% rownames(this_res)){
        this_res <- this_res[x, , drop=FALSE]
      } else if(x_safe %in% rownames(this_res)){
        this_res <- this_res[x_safe, , drop=FALSE]
      } else {
        stop(sprintf("Metabolite '%s' not found in ordinal model coefficients. Available: %s", 
                     x, paste(rownames(this_res), collapse=", ")))
      }
    } else {
      this_res <- this_res[2, , drop=FALSE]
    }
    
    this_res <- as.data.frame(this_res)
    rownames(this_res) <- NULL
    
    this_res$var <- x
    this_res$term <- x
    this_res$formula <- formula_str
    
    return(list(res = this_res, model = this_fit))
  })
  
  # extract models and results separately
  models <- lapply(results_list, function(x) x$model)
  names(models) <- mets
  
  univ_stats <- lapply(results_list, function(x) x$res) %>%
    bind_rows()
  
  # merge with rowData
  rowdata_df <- D %>% rowData() %>% data.frame()
  
  if("name" %in% colnames(rowdata_df)){
    rowdata_df$join_key <- make.names(rowdata_df$name)
  } else {
    rowdata_df$join_key <- rownames(rowdata_df)
  }
  
  univ_stats <- univ_stats %>% 
    left_join(rowdata_df, by = c("var" = "join_key")) %>%
    select(var, term, formula, estimate, std.error, df, statistic, p.value, everything()) %>%
    dplyr::rename(orig_name = name)
  
  # adjust pvalues
  univ_stats %<>% mutate(p.adj = p.adjust(p.value, method='BH'))
  univ_stats <- univ_stats[order(univ_stats$p.adj), ]
  
  univ_stats <- as.data.frame(univ_stats)
  rownames(univ_stats) <- NULL
  
  ## construct output groups variable (if outcome is factor)
  if (is.factor(Ds[[outcome]])) {
    outgroups <- levels(Ds[[outcome]])
  } else {
    outgroups <- NULL
  }
  
  ## create output list with required structure (matching mt_stats_univ_lm)
  output_list <- list(
    table = univ_stats,
    formula = if(!is.null(conf_formula)) conf_formula else "",
    name = stat_name,
    lstobj = models,
    groups = outgroups,
    samples.used = !is.na(Ds[[outcome]]),
    outcome = outcome
  )
  
  ## construct formula string for logging
  if(!is.null(conf_formula) & !is.null(random_effects)){
    log_formula <- sprintf("%s ~ metabolite + %s + %s", outcome, conf_formula, random_effects)
  } else if(!is.null(conf_formula)){
    log_formula <- sprintf("%s ~ metabolite + %s", outcome, conf_formula)
  } else if(!is.null(random_effects)){
    log_formula <- sprintf("%s ~ metabolite + %s", outcome, random_effects)
  } else {
    log_formula <- sprintf("%s ~ metabolite", outcome)
  }
  
  ## construct model description for logging
  if(outcome_type == 'binary'){
    if(!is.null(random_effects)){
      model_desc <- "binary outcome, logistic mixed model"
    } else {
      model_desc <- "binary outcome, logistic regression"
    }
  } else if(outcome_type == 'numeric'){
    if(!is.null(random_effects)){
      model_desc <- "numeric outcome, linear mixed model"
    } else {
      model_desc <- "numeric outcome, linear regression"
    }
  } else if(outcome_type == 'ordinal'){
    if(!is.null(random_effects)){
      model_desc <- "ordinal outcome, cumulative link mixed model"
    } else {
      model_desc <- "ordinal outcome, cumulative link model"
    }
  }
  
  ## add status information & plot
  funargs <- mti_funargs()
  D %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("mixed association analysis script</br>%s</br>%s", model_desc, log_formula),
      output = output_list
    )
  ## RETURN
  D
  
}
