#' Compare SE objects for equality for each step in a pipeline
#'
#' Compares reference objects for each step in a pipeline to the individual results of each step in a pipeline
#' being run in real-time.
#' Currently, plot objects and other non-stat table forms of output are ignored. Comparisons will only be made between two SE objects
#' with identical function calls.
#'
#' @param file Name of test script file.
#' @param ref_path Path to dir containing step-wise \code{SummarizedExperiment} input objects.
#' @param obj_pre Prefix of maplet pipeline test objects.
#'
#' @return Data frame with columns: step, fun, test, and result.
#'
#' @examples
#' \dontrun{
#'      mtt_compare_stepwise(file = "stats_step_test.R", ref_path="REFERENCE_OBJECTS/", obj_pre = "stats_pipe")
#' }
#' @export
mtt_compare_stepwise <- function(file, ref_path, obj_pre){

  code <- parse(file)
  tokens <- as.list(code)

  step_res_lst <- list()

  # step through each step in the pipeline
  for (ii in seq_along(tokens)) {
    part <- tokens[[ii]]

    # parse and call functions
    if(length(part) > 2){
      if(as.character(part[[2]]) == "D"){
        if(is.call(part[[3]]) && startsWith(as.character(part[[3]][[1]]), "mt_")){

          # parse function name and arguments
          fun_arg_lst <- part[[3]]
          mt_fun <- as.character(fun_arg_lst[1]) # function name
          fun_arg_lst <- fun_arg_lst[-1]
          args_as_char <- as.character(fun_arg_lst)
          if(as.character(part[[1]]) == "<-"){
            # if system file provided, get absolute path
            if(any(startsWith(args_as_char, "system.file"))){
              idx <- which(startsWith(args_as_char, "system.file"))
              fun_arg_lst[idx] <- args_as_char[idx] %>% str2lang() %>% eval()
            }
            step <- 1
          }
          # if system file provided, get absolute path
          if(any(startsWith(args_as_char, "system.file"))){
            idx <- which(startsWith(args_as_char, "system.file"))
            fun_arg_lst[idx] <- args_as_char[idx] %>% str2lang() %>% eval()
          }
          arg_lst <- fun_arg_lst %>% as.list() # argument values
          #names(arg_lst) <- names(part[[3]])[-1] # argument names
          if(as.character(part[[1]]) == "%<>%") arg_lst$D <- D_step

          # call the maplet function
          mti_logstatus(glue::glue("   Running step {step} of test script."))
          D_step <- do.call(mt_fun, args = arg_lst)

          # load reference object
          mti_logstatus(glue::glue("   Loading step {step} reference object."))
          ref_file <- file.path(ref_path, paste0(obj_pre, step, ".rds"))
          D_ref <- readRDS(ref_file)

          # compare if SEs equal up until this step, save result
          mti_logstatus(glue::glue("   Comparing step {step} of test script to reference object."))
          step_test_res <- mtt_equal(D_step, D_ref)
          step_test_res <- cbind(step=step, fun=mt_fun, step_test_res)
          step_res_lst[[step]] <- step_test_res

          step <- step + 1
        } # end of startsWith(.,"mt_") if-block
      } # end of part[[2]] == "D" if-block
    }# end length(part)>2 block
  } # end of for-loop

  # return step_res_list as data frame
  step_res_df <- do.call(rbind, step_res_lst)
  step_res_df


}
