#' Mixed Association Analysis
#' 
#' ADD DESCRIPTION
#' 
#' @param D \code{SummarizedExperiment} input.
#' @param stat_name Name of the statistical result.
#' @param outcome Outcome to be used for response.
#' @param outcome_type Outcome type - one of numeric/binary/ordinal.
#' @param conf_formula Confounders to be corrected for. Default: NULL.
#'
#' @return $results$output: List of values including 'estimate', 'std_error', 'df', 'statistic', and 'p_value'.
#'
#' @examples
#' \dontrun{# ADD EXAMPLE}
#'
#'
#' @author RB
#' 
#' @export
mt_stats_mixed_association_analysis <- function(D,
                                                stat_name,
                                                outcome,
                                                outcome_type,
                                                conf_formula=NULL
){
  
  
  # merge data with sample info
  Ds <- D %>% maplet:::mti_format_se_samplewise() # NOTE: No explosion of dataset size, no gather() - 6/2/20, JK
  
  # metabolites in data
  mets <- D %>% assay() %>% rownames()
  
  # loop over metabolites
  univ_stats <- lapply(mets, function(x){
    # formula for this metabolite
    if(is.null(conf_formula)==F){
      this_formula <- as.formula(glue('{outcome} ~ {x} + {conf_formula}'))  
    } else{
      this_formula <- as.formula(glue('{outcome} ~ {x}'))  
    }
    if(outcome_type=='binary'){
      # turn the outcome variable into factor
      Ds %<>% mutate(!!sym(outcome) := as.factor(as.matrix(!!sym(outcome))))
      # logistic regression
      this_fit <- lme4::glmer(this_formula, data =Ds, family = binomial,
                              nAGQ = 2)
      # results summary
      this_res <- this_fit %>% summary() %>% coefficients() %>% data.frame()
      names(this_res) <- c('estimate', 'std_error', 'statistic', 'p_value')
    } else if(outcome_type=='numeric'){
      # turn the outcome variable into factor
      Ds %<>% mutate(!!sym(outcome) := as.numeric(as.matrix(!!sym(outcome))))
      # linear regression
      this_fit <- lme4::lmer(this_formula, data =Ds)
      # results summary
      this_res <- this_fit %>% lmerTest::as_lmerModLmerTest() %>% 
        summary() %>% coefficients() %>% data.frame()
      
      names(this_res) <- c('estimate', 'std_error', 'df', 'statistic', 'p_value')
    }
    this_res <- this_res [2, ]
    this_res %<>% mutate(analyte=rownames(this_res),
                         outcome=outcome, 
                         covariates=conf_formula)
    # order output columns
    this_res %<>% select(analyte, outcome, everything())
    return(this_res)
  }) %>% # create data from of results
    do.call(rbind,.) %>% data.frame() %>% 
    # bind rowData
    bind_cols(D %>% rowData() %>% data.frame() %>% select(name, everything()))
  # adjust pvalues
  univ_stats %<>% mutate(adj_p = p.adjust(p_value, method='BH'))
  # order by adjusted p-values
  univ_stats <- univ_stats[order(univ_stats$adj_p), ]
  
  ## add status information & plot
  funargs <- mti_funargs()
  D %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("Performed mixed association analysis for outcome %s.", outcome),
      output = univ_stats
    )
  ## RETURN
  D
  
}
