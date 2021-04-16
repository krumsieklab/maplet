# # maplet
#
# Helper functions.
#
# last update: 2020-09-17
# authors: JK,MB, KC
#

#' Retrieve ML model results by name
#'
#' Returns a list of ML results given a result name
#'
#' @param D SummarizedExperiment object
#' @param name the name of the ML result of interest
#'
#' @return res_list: list containing the output and output2 lists for a given name
#'
#' @noRd
mti_get_ml_res_by_name <- function(D, name){
 # code is identical to mtm_get_stat_by_name, but returns output and output2 lists
  stopifnot("SummarizedExperiment" %in% class(D))

  if(! ("results" %in% names(metadata(D)))){
    stop("no results element found in D")
  }

  stats <- maplet:::mti_res_get_path(D,"ml")

  if(length(stats) == 0){
    stop("no stats element found in D")
  }

  names  <- stats %>% purrr::map_chr(~.x$output$name)
  output <- which(names == name)

  if(length(output) == 0){
    stop("stat element with name ", name, " does not exist")
  }
  if(length(output)  > 1){
    stop("there are multiple stat elements with name ", name)
  }

  res_list <- list()
  res_list$output <- stats[[ output ]]$output
  res_list$output2 <- stats[[ output ]]$output2

  res_list

}


#' Safely add to end of list, even if list does not exist
#'
#' Adds an item to the end of a list, even if the list is empty or NULL
#'
#' @param lst list to be extended
#' @param element element to be added
#' @param oname name of the new element
#'
#' @return extended list
#'
#' @noRd
mti_add_to_list <- function(lst, element, oname) {
  # init if needed
  if (is.null(lst) || length(lst)==0) lst=c()
  # add and return
  lst[[length(lst)+1]] = element
  names(lst)[length(lst)] = oname

  lst
}


#' Concatenate data matrix with sample annotations
#'
#' Returns data frame as samples X variables, merges all sample annotations, and adds sample rownames as "merge.primary" field
#'
#' @param D SummarizedExperiment input
#'
#' @return Concatenated dataframe
#'
#' @noRd
mti_format_se_samplewise <- function(D){
  # coldata and assay cannot have overlapping names
  # (this should be caugh earlier in the pipeline, but here is where it causes trouble)
  inters <- intersect(colnames(colData(D)), rownames(D))
  if (length(inters)>0) {
    stop(sprintf("There are metabolites and colData variables with the same name: %s", paste0(inters, collapse = ", ")))
  }
  # cbind
  cbind(colData(D),
        t(assay(D))) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("merge.primary")
}


#' Helper function to generate metadata(.)$results entry
#'
#' Automatically generates fun and args fields, and sends logtxt through logmsg()
#'
#' @param D \code{SummarizedExperiment} input
#' @param funargs output from mti_funargs() that should be collected by calling function
#' @param logtxt log text describing what the function did
#' @param output output structure of that function (e.g. plot, statistical results)... default is NULL (i.e. no output)
#' @param output2 optional second output... default is NULL
#'
#' @return list ready to be stored in metadata()$results
#'
#' @noRd
mti_generate_result <- function(
  D,
  funargs,
  logtxt="",
  output=NULL,
  output2=NULL
) {
  
  # ensure structure of funargs
  stopifnot("fun" %in% names(funargs))
  stopifnot("args" %in% names(funargs))

  this.uuid = uuid::UUIDgenerate()
  # check if assay updated
  save_assays = mti_get_setting(D, "save_all_assays")
  if(save_assays){
    if(assays(D) %>% length() != 0){
      metadata(D)$assays %<>% mti_assay_ptr(res_assay=assay(D))
      assay_ptr <- metadata(D)$assays$head_ptr
    }else{
      assay_ptr <- NULL
    }
  }else{
    assay_ptr <- NULL
  }
  
  # assemble list
  metadata(D)$results %<>% mti_add_to_list(
    list(
      fun=funargs$fun,
      args=funargs$args,
      logtxt=mti_logmsg(logtxt),
      uuid=this.uuid,
      output=output,
      output2=output2,
      assay_ptr=assay_ptr
    ),
    oname = paste(paste(funargs$fun,collapse = "_"), this.uuid, sep = ".")
  )
  
  D
}


#' Return pointer to current assay
#' 
#' Checks whether assay has been updated. If true, adds new assay to $results$assay and assigns new head_assay and returns new 
#' pointer. Else, returns pointer to current head_assay.
#' 
#' @param assays List containing $head_ptr (pointer to current version of assay) and assay_lst (list of all versions of assay).
#' @param res_assays Current assay data frame from latest function call.
#' 
#' @return assays List with updated $head_ptr and assay_lst (if assay changed).
#' 
#' @noRd
mti_assay_ptr <- function(assays, res_assay){
  
  if(is.null(assays)) assays <- list(head_ptr=NULL, assay_lst=list())
  
  a_uuid = uuid::UUIDgenerate()
  
  # no head, add current assay and update pointer
  if(is.null(assays$head_ptr) && !is.null(res_assay)){
    assay_ptr <- paste0("assay.", a_uuid)
    assays$assay_lst[[assay_ptr]] <- res_assay
    assays$head_ptr <- assay_ptr
  }else{
    # get head_assay and compare to res_assay
    head_assay <- assays$assay_lst[[assays$head_ptr]]
    assay_changed <- !identical(head_assay, res_assay)
    
    # if assay changed, add to list and update
    if(assay_changed){
      assay_ptr <- paste0("assay.", a_uuid)
      assays$assay_lst[[assay_ptr]] <- res_assay
      assays$head_ptr <- assay_ptr
    }
  }
  
  assays
  
}


#' Return list of all arguments supplied to function
#'
#' Removes first argument $D if it exists.
#' Removes "mt" from function name list
#' Returns list of $fun, already exploded (e.g. c("plots","boxplot")), and $args
#'
#' @param ... Not actually used (TODO delete?). Accesses arguments of parent function
#'
#' @return List of arguments, see description.
#'
#' @noRd
mti_funargs <- function(...) {
  call <- evalq(match.call(expand.dots = FALSE), parent.frame(1))
  formals <- evalq(formals(), parent.frame(1))

  for(i in dplyr::setdiff(names(formals), names(call)))
    call[i] <- list( formals[[i]] )

  # assemble results
  raw <-  as.list(match.call(sys.function(sys.parent(1)), call))
  if(typeof(raw[[1]])=="language") raw[[1]] <- as.character(raw[[1]])[[3]]
  res <- list(
    fun = strsplit(gsub(".*::", "", as.character(raw[[1]])), '_')[[1]],
    args = raw[-1]
  )
  # remove "mt"
  res$fun = res$fun[res$fun!="mt"]
  # make sure D does not exist in args
  res$args$D = NULL
  # return
  res

}


#' Turns ... argument into string, key1=value1, key2=value2 etc.
#'
#' @param ... Arbitary list of input arguments
#'
#' @return String key/value list.
#'
#' @noRd
mti_dots_to_str <- function(...) {
  l = eval(substitute(alist(...)))
  paste(sapply(names(l), function(k){sprintf('%s=%s',k,as.character(l[[k]]))}), collapse = ', ')
}


#' Confounder correction for a SummarizedExperiment
#'
#' Used specifically for boxplot function to generate confounder-corrected residuals for plotting.
#'
#' @param D \code{SummarizedExperiment} input
#' @param formula formula for correction
#'
#' @returns SummarizedExperiment with corrected data
#'
#' @noRd
mti_correctConfounder <- function(D, formula){
  d <- D %>% mti_format_se_samplewise() # NOTE: No explosion of dataset size, no gather() - 6/2/20, JK
  d_cor <- rownames(D) %>%
    purrr::map_dfc(function(m){
      f   <- stats::update.formula(formula, stringr::str_c(m, "~."))
      mod <- stats::lm(f, data = d, na.action = na.exclude)
      res <- stats::resid(mod)
      res
    }) %>%
    stats::setNames(rownames(D)) %>%
    as.matrix() %>% t()
  colnames(d_cor) <- colnames(D)
  assay(D)        <- d_cor
  D
}


#' Add left- and right-aligned x axis labels to ggplot
#'
#' ggplot is missing the functionality to add x-axis labels that are left- and right-aligned. This function adds those.
#'
#' @param ggplot object
#' @param left text on the left
#' @param right text on the right
#'
#' @examples
#' \dontrun{# To demonstrate function behavior on empty plot:
#' (ggplot(mapping = aes(x=0,y=0)) + geom_point()) %>% mti_add_leftright_gg("left label","right label")
#' }
#'
#' @return ggplot object
#'
#' @author MB, JK
#' @noRd
mti_add_leftright_gg <- function(gg, left, right) {

  # with backward compatibility
  if (utils::compareVersion(as.character(utils::packageVersion("ggplot2")),"3.3.0")>=0) { # at least 3.3.0
    # ggplot2 version >= 3.3.0
    ggbld <- ggplot2::ggplot_build(gg)
    xticks =  ggbld$layout$panel_params[[1]]$x$get_breaks() # needed for >=3.3.0
    xticks.minor = ggbld$layout$panel_params[[1]]$x$get_breaks_minor() # needed for >=3.3.0
    xlims = ggbld$layout$panel_params[[1]]$x.range

    # add positions of annotation labels
    xticks = c(xlims[1], xticks, xlims[2])
    # get breaks labels
    xtlabs = ggbld$layout$panel_params[[1]]$x$get_labels() # needed for >=3.3.0

  } else {
    # ggplot2 version < 3.3.0
    ggbld <- ggplot2::ggplot_build(gg)
    xticks = ggbld$layout$coord$labels(ggbld$layout$panel_params)[[1]]$x.major_source # needed for <3.3.0
    xticks.minor = ggbld$layout$coord$labels(ggbld$layout$panel_params)[[1]]$x.minor_source # needed for <3.3.0
    xlims = ggbld$layout$panel_params[[1]]$x.range

    # add positions of annotation labels
    xticks = c(xlims[1], xticks, xlims[2])
    # get breaks labels
    xtlabs = ggbld$layout$coord$labels(ggbld$layout$panel_params)[[1]]$x.labels # needed for <3.3.0
  }

  # align with \n
  txt <- c(left, right)
  txt = paste0("\n", txt)
  xtlabs = paste0(xtlabs, "\n")
  xtlabs = c(txt[1], xtlabs, txt[2])

  # return
  gg + ggplot2::scale_x_continuous(breaks = xticks, labels = xtlabs, minor_breaks = xticks.minor)

}


# extracts variables from a list of terms
mti_extract_variables <- function(lst) {
  # filter down only to the variables needed for plotting
  # need to parse x and ... list
  # browser()
  # q <- quos(...)
  vars = c()
  if (length(q) > 0) {
    vars <- lst %>% unlist() %>% lapply(function(x){x %>% all.vars()}) %>% unlist() %>% as.vector()
  }
  vars <- unique(vars)

  # return
  vars
}


#' ADD TITLE
#'
#' ADD DESCRIPTION
#'
#' @param x ADD PARAM DESCRIPTION
#'
#' @return ADD RETURN DESCRIPTION
#'
#' @author whoWroteIt?
#' @noRd
fixorder = function(x){o= unique(as.character(x)); gdata::reorder.factor(x, new.order=o)} # fix order of a factor


#' ADD TITLE
#'
#' ADD DESCRIPTION
#'
#' @param Ds Concatenated dataframe returned by mti_format_se_samplewise
#' @param sample_filter term which samples to filter to first
#'
#' @return logical vector of samples to keep
#'
#' @author JK, JZ, KC
#' @noRd
mti_filter_samples <- function(Ds, filter_q, num_samp){

  Ds <- Ds %>%
    dplyr::mutate(tmpsamplenum = 1:nrow(Ds)) %>%
    dplyr::filter(!!filter_q) %>%
    droplevels()
  # message("filter metabolites: ", metab_filter_q, " [", nrow(stat), " remaining]")
  # did we leave 0 rows?
  if (nrow(Ds)==0) stop("Filtering left 0 rows")
  if (nrow(Ds)==num_samp) maplet:::mti_logwarning('filtering did not filter out any samples')

  # store used samples
  samples.used <- rep(F, num_samp)
  samples.used[Ds$tmpsamplenum] <- T

  # return
  samples.used

}

#' Calculate evaluation measures for classifier predictions
#'
#' Calculates the following evaluation measures: sensitivity (sens), specificity (spec),
#'  accuracy (acc), positive-predictive value (ppv), and F1-measure (F1)
#'
#' @param tpfn a list of the four confusion matrix outcomes: True Positive (TP), False Positive (FP), False Negative (FN), and
#'  True Negative (TN)
#'
#' @return list of evaluation measures
#' @noRd
mti_measures <- function(tpfn) {
  # construct sensitivity, specificity, accuracy, PPV and F1 score; return as list
  x = list(
    sens = tpfn$TP/(tpfn$TP+tpfn$FN),
    spec = tpfn$TN/(tpfn$TN+tpfn$FP),
    acc = (tpfn$TP+tpfn$TN)/(tpfn$TP+tpfn$FN+tpfn$FP+tpfn$TN),
    PPV = tpfn$TP/(tpfn$TP+tpfn$FP)
  )
  x$F1 = 2*(x$sens*x$spec)/(x$sens+x$spec)
  x
}

#' Calculate Confusion Matrix Outcomes
#'
#' Calculates the following four confusion matrix outcomes: True Positive (TP), False Positive (FP), False Negative (FN), and
#'  False Positive (FP)
#'
#' @param trueclass a boolean vector of class labels
#' @param predclass a boolean vector of predicted classes
#'
#' @return list of confusion matrix outcomes
#' @noRd
mti_TPFN <- function(trueclass, predclass) {
  list(
    TP = sum(trueclass & predclass),
    FP = sum(!trueclass & predclass),
    FN = sum(trueclass & !predclass),
    TN = sum(!trueclass & !predclass)
  )
}

#' Get Evaluation Measures List
#'
#' Get a list of evaluation measures given a vector of predicted class probabilities and class labels. The following
#' evaluation measures are included in the list: sensitivity (sens), specificity (spec), accuracy (acc), positive-predictive
#' value (ppv), and F1-measure (F1).
#'
#' @param trueclass a boolean vector of class labels
#' @param predprob a vector of predicted class probabilities
#'
#' @return
#' @noRd
mti_get_measures_list <- function(trueclass, predprob) {
  # initialize arrays
  sensvals = c()
  specvals = c()
  ppvvals = c()
  accvals = c()
  f1vals = c()
  # loop over all predprob values as cutoff
  predprob.sorted = sort(predprob)
  # trick: add -Inf to make the curve work
  predprob.sorted = c(-Inf, predprob.sorted)
  # loop over all values
  for (i in 1:length(predprob.sorted)) { # could also be done via sapply()
    # do the cut
    predclass = predprob > predprob.sorted[i]
    # quality measures
    m = mti_measures(mti_TPFN(trueclass, predclass))
    sensvals[i] = m$sens
    specvals[i] = m$spec
    ppvvals[i] = m$PPV
    accvals[i] = m$acc
    f1vals[i] = m$F1
  }
  # calculate AUC under ROC curve
  # we have to take the negative, because we built the curve the backwards
  AUC = -pracma::trapz(1-specvals,sensvals)

  predprob.sorted[1] <- 0
  # return list
  list(sensvals=sensvals, specvals=specvals, AUC=AUC, ppvvals=ppvvals, thresholds=predprob.sorted, accvals=accvals, f1vals=f1vals)
}

#' Check if Data is Logged
#'
#' Check to see if either mt_pre_trans_log or mt_load_flag_logged was called in the pipeline. If mt_pre_trans_log called multiple
#' times or logging reversed, throw a warning and set is_logged to FALSE.
#'
#' @param D SummarizedExperiment
#'
#' @return is_logged (boolean)
#' @noRd
mti_check_is_logged <- function(D){

  is_logged <- FALSE

  # get names of the functions called
  called_functions <- names(metadata(D)$results)

  # if data is logged more than once, throw warning and return is_logged = FALSE
  num_times_logged <- called_functions %>% startsWith("pre_trans_log") %>% sum()
  if(num_times_logged > 1){
    warning("Data logged multiple times! Cannot determine log status. is_logged will be set to FALSE.")
    return(is_logged)
  }

  # if logging is reversed, throw warning and return is_logged = FALSE
  if(num_times_logged==1){
    log_index <- called_functions %>% startsWith("pre_trans_log") %>% which()
    exp_after_logged <- called_functions[(first_log_index+1):length(called_functions)] %>% startsWith("pre_trans_exp") %>% any()
    if(exp_after_logged){
      warning("Logging of data has been reversed! Cannot determine log stats. is_logged will be set to FALSE.")
      return(is_logged)
    }else{
      is_logged <- TRUE
    }
  }

  if(any(startsWith(called_functions, "load_flag_logged"))){
    # if pre_trans_log called AND load_flag_logged called, return is_logged as False
    if(is_logged){
      is_logged <- FALSE
      return(is_logged)
    }else{
      is_logged <- TRUE
    }
  }

  is_logged

}
