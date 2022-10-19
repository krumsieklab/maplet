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
#'    }
#' }
#'
#' @return
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

  VALID_RESULT_NAMES <- c("df", "cd", "rd")

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

  # ADD: code for storing plots (?)

  # ADD: code for storing statistical results (?)

  funargs <- mti_funargs()
  D %<>%
    mti_generate_result(
      funargs = funargs,
      #logtxt = cat(paste0(c("The following code block was executed: ", body(code)), collapse ='\n'))
      logtxt = paste0(c("The following code block was executed: ", body(elisa_code)), collapse ='\n')
    )

  # return
  D

}
