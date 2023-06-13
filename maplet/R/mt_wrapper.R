#' Wrapper function for non-maplet code blacks
#'
#' This function allows the user to apply a block of non-maplet code on a maplet-SE object while
#' keeping the object within a maplet pipeline and ensuring the SE object is not altered in a way
#' that would make it maplet-incompatible. The code block must satisfy the following conditions:
#' \itemize{
#'    \item Must be passed as a function
#'    \item Function must take exactly one argument (corresponding to D)
#'    \item Function must not change dimensions of assay (this must be done with mt_modify_filter)
#'    \item Function must not filter the rows of colData or rowData (this must be done with
#'       mt_modify_filter); removal of columns is allowed
#'    \item Function must not change the colnames or rownames of the assay
#'    \item Function must return a list containing elements corresponding to these accepted names:
#'    \itemize{
#'       \item{df - assay}
#'       \item{cd - colData data frame}
#'       \item{rd - rowData data frame}
#'       \item{st - a named list containing the following objects:}
#'       \itemize{
#'          \item{table - statistical data frame; must contain the columns: var, statistic, and p.value}
#'          \item{name - name of the statistical comparison, must be unique in the pipeline}
#'          \item{samples.used - logical vector of samples used in the statistical comparison}
#'       }
#'       \item{pt - a ggplot object or list of ggplot objects}
#'    }
#' }
#'
#' @return Named list of objects. See description.
#'
#' @examples
#' \dontrun{# remove unnecessary columns in colData
#' ...  %>%
#'  mt_wrapper(code = my_code <- function(D){
#'     cd <- colData(D)
#'     cd <- dplyr::select(-c("GROUP_ID", "Age_at_diag", "C_Response", "Pathology_info"))
#'     list(cd = cd)
#'  }) %>% ...}
#'
#' @author KC
#'
#' @export
mt_wrapper <- function(D, code){

  VALID_RESULT_NAMES <- c("df", "cd", "rd", "pt", "st")

  # verify D is SE
  stopifnot("SummarizedExperiment" %in% class(D))

  # verify code is a function that takes exactly one argument
  if(!("function" %in% class(code))) stop("Argument \'code\' must be a function.")
  arg <- formals(code)
  if(length(arg) > 1) stop("Function must take exactly one argument.")


  # try to run code - crash if error, continue if warning
  tryCatch({
    results <- code(D)
  },
  warning = function(warn) {
    warning(glue::glue("Code block produced the following warning: {warn}"))
  },
  error = function(err) {
    stop(glue::glue("Code block crashed with the following error:\n {err}"))
  })

  # check results list is valid
  res_names <- names(results)
  if(!all(res_names %in% VALID_RESULT_NAMES)) stop("One or more element names in list are invalid. See documentation.")

  # perform checks, then update assay
  if("df" %in% res_names){
    if(dim(results$df) != dim(assay(D))) stop("Dimensions of new assay do not match old assay!")
    if(colnames(results$df) != colnames(assay(D))) stop("Column names of new assay do not match old assay!")
    if(rownames(results$df) != rownames(assay(D))) stop("Column names of new assay do not match old assay!")
    assay(D) <- results$df

  }

  # perform checks, then update colData
  if("cd" %in% res_names){
    if(dim(results$cd)[1] != dim(colData(D))[1]) stop("Dimensions of new colData data frame do not match old colData data frame!")
    colData(D) <- S4Vectors::DataFrame(results$cd)
  }

  # perform checks, then update rowData
  if("rd" %in% res_names){
    if(dim(results$rd)[1] != dim(rowData(D))[1]) stop("Dimensions of new rowData data frame do not match old rowData data frame!")
    rowData(D) <- S4Vectors::DataFrame(results$rd)
  }

  # perform checks, add statistical results to output
  if('st' %in% res_names){
    # check all required list elements present
    st_req_names <- c('table', 'name', 'samples.used')
    req_missing <- setdiff(st_req_names, names(results$st))
    if(length(req_missing)!=0) stop(glue::glue("The following required names are missing from returned stats list: {req_missing}"))

    # check all required columns present
    tab_req_cols <- c('var', 'statistic', 'p.value')
    req_missing <- setdiff(tab_req_cols, colnames(results$st$table))
    if(length(req_missing)!=0) stop(glue::glue("The following required columns are missing from returned statistical table: {glue::glue_collapse(req_missing, sep = ', ')}"))

    # check samples.used is equal to the number of samples and boolean
    if(length(results$st$samples.used) != ncol(D)) stop("Length of returned vector samples.used must be equal to the number of samples.")
    if(!is.logical(results$st$samples.used)) stop("Returned vector samples.used must be logical.")

    # check stat_name is unique
    if (results$st$name %in% unlist(maplet::mtm_res_get_entries(D, "stats") %>% purrr::map("output") %>% purrr::map("name"))) stop(sprintf("stat element with stat_name '%s' already exists",results$st$name))

    output = results$st
  }else{
    output = NULL
  }

  # perform checks, add plots to output2
  if('pt' %in% res_names){
    # check if object or element list can be plotted
    if(("gg" %in% class(results$pt)) | ("gg" %in% class(results$pt[[1]]))){
      output2 = results$pt
    }else{
      stop("Returned plot object must by a ggplot object or a list of ggplot objects.")
    }

  }else{
    output2 = NULL
  }


  funargs <- mti_funargs()
  D %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = paste0(c("The following code block was executed: ", body(code)), collapse ='\n'),
      output = output,
      output2 = output2
    )

  # return
  D

}
