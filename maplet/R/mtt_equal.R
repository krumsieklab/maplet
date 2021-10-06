#' Compare two SE objects for global equality
#'
#' Compares the following attributes of two SE objects:
#' \itemize{
#' \item assay - checks assay data frames are identical
#' \item rowData - checks rowData data frames are identical
#' \item colData - checks colData data frames are identical
#' \item stat_len - checks both objects have the same number of statistical result tables
#' \item stat_# - checks statistical result tables are identical, # is a number identifying the statistical result
#' \item pathways_len - checks both objects have the same number of pathway annotation tables
#' \item pathways_# - checks pathway annotation data frames are identical, # is a number identifying the pathway annotation data frame
#' }
#' Currently, plot objects and other non-stat table forms of output are ignored. Comparisons will only be made between two SE objects
#' with identical function calls.
#'
#' @param D1 First \code{SummarizedExperiment} input. When running tests, assume D1 is the reference SE object.
#' @param D2 Second \code{SummarizedExperiment} input.
#'
#' @return Data frame with columns: test and result.
#'
#' @examples
#' \dontrun{
#'      # compare two maplet pipleine objects
#'      mtt_equal(D1, D2)
#' }
#' @export
mtt_equal <- function(D1, D2){

  # check D1 and D2 are SummarizedExperiment objects
  if("SummarizedExperiment" %in% class(D1)==F) stop("D1 is not a SummarizedExperiment object.")
  if("SummarizedExperiment" %in% class(D2)==F) stop("D2 is not a SummarizedExperiment object.")

  # check that function calls are exactly the same - otherwise can't compare
  d1_calls <- mtm_res_get_uuids(D1) %>% dplyr::select(-c(UUIDs)) %>% set_rownames(.$call_order)
  d2_calls <- mtm_res_get_uuids(D2) %>% dplyr::select(-c(UUIDs)) %>% set_rownames(.$call_order)
  if(!identical(d1_calls, d2_calls)) stop("D1 and D2 must have the same function calls in the same order.")

  test_list <- list()

  if(length(D1) != 0){
    # compare the three main data frames - rowData, colData, and assay
    mti_logstatus("   Comparing the assay, rowData, and colData data frames.")
    test_list[["assay"]] <- ifelse(identical(assay(D1), assay(D2)), "Pass", "Fail")
    test_list[["rowData"]] <- ifelse(identical(as.data.frame(rowData(D1)), as.data.frame(rowData(D2))), "Pass", "Fail")
    test_list[["colData"]] <- ifelse(identical(as.data.frame(colData(D1)), as.data.frame(colData(D2))), "Pass", "Fail")
  }

  # compare any stats tables
  d1_stat_tables <- mtm_res_get_entries(D1, name_list = c("stats")) %>% purrr::map("output") %>% purrr::map("table")
  d2_stat_tables <- mtm_res_get_entries(D2, name_list = c("stats")) %>% purrr::map("output") %>% purrr::map("table")
  if(length(d1_stat_tables) == length(d2_stat_tables)){
    test_list[["stats_len"]] <- "Pass"
    if(length(d1_stat_tables)!=0){
      mti_logstatus("   Comparing the statistical results tables.")
      for(i in 1:length(d1_stat_tables)){
        d1_stat_df <- d1_stat_tables[[i]]
        d2_stat_df <- d2_stat_tables[[i]]
        test_list[[paste0("stats_", i)]] <- ifelse(identical(d1_stat_df, d2_stat_df), "Pass", "Fail")
      }
    }
  }else{
    test_list[["stats_len"]] <- "Fail"
  }

  # compare any pathway annotation data frames
  d1_pw_anno_dfs <- metadata(D1)$pathways
  d2_pw_anno_dfs <- metadata(D2)$pathways
  if(length(d1_pw_anno_dfs) == length(d2_pw_anno_dfs)){
    test_list[["pathways_len"]] <- "Pass"
    if(length(d1_pw_anno_dfs)!=0){
      mti_logstatus("   Comparing the pathway annotation data frames.")
      for(i in 1:length(d1_pw_anno_dfs)){
        d1_pw_df <- d1_pw_anno_dfs[[i]]
        d2_pw_df <- d2_pw_anno_dfs[[i]]
        test_list[[paste0("pathways_", i)]] <- ifelse(identical(d1_pw_df, d2_pw_df), "Pass", "Fail")
      }
    }
  }else{
    test_list[["pathways_len"]] <- "Fail"
  }

  # compare assay lists, if present
  if(!is.null(metadata(D1)$assays) && !is.null(metadata(D2)$assays)){
    mti_logstatus("   Comparing the assay version lists.")
    # DO SOMETHING!

  }

  test_res_df <- data.frame(test=names(test_list), result=unlist(test_list)) %>% set_rownames(1:length(test_list))

  test_res_df

}
