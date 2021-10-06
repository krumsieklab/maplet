#' Export statistical results from pipeline into an Excel file
#'
#' Writes out the statistics data frame with p-values, adjusted p-values, fold changes, test statistics etc. into an Excel sheet.
#'
#' @param D \code{SummarizedExperiment} input.
#' @param file Output filename to write to.
#' @param stat_list Names of one or more statistical comparison to be written out. If NULL, will export all.
#' @param sort_by_p Automatically sort by p-values? Default: F.
#' @param feat_col OPTIONAL. Name of rowdata column to include in the output. Default: NULL.
#'
#' @return Does not change the \code{SummarizedExperiment} object.
#'
#' @examples
#' \dontrun{# Write out all results
#' ... %>% mt_write_stats(file="results.xlsx") %>%
#' # Write out specific result]
#' ... %>% mt_write_stats(file="results.xlsx", stat_list="comp1") %>%}
#'
#' @author JK, RB
#'
#' @export
mt_write_stats <- function(D, file, stat_list=NULL, sort_by_p=F, feat_col=NULL) {

  # verify input arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  stopifnot(is.character(file))

  # get all stats entries
  S <- D %>% mtm_res_get_entries("stats")
  allcomps <- S %>% purrr::map("output") %>% purrr::map("name") %>% unlist()

  # restrict to one or output all?
  if (!is.null(stat_list)) {
    # find the ones to filter to
    m <- match(stat_list, allcomps)
    # throw error if any are not found
    if (any(is.na(m))) {
      stop(sprintf("Cannot find the following stats entries to export: %s", paste0(stat_list[is.na(m)], collapse = ", ")))
    }
    # subset
    S <- S[m]
  }

  # export all
  wb <- openxlsx::createWorkbook()
  for (i in 1:length(S)) {
    # add to workbook
    df <- S[[i]]$output$table
    # if feat_col is provided and it matches a column in row data add it to the output
    # else just silently skip this step
    if(is.null(feat_col)==F && is.na(match(feat_col, names(rowData(D))))==F){
        df <- cbind.data.frame(df, feat_col=unlist(data.frame(rowData(D))[feat_col]))
    }
    # sort?
    if (sort_by_p) {
      df %<>% arrange(p.value)
    }
    # add direction?
    if ("groups" %in% names(S[[i]]$output) && length(S[[i]]$output$groups)==2) {
      # generate vector of directions for indexing
      inds <- as.numeric(df$estimate>0)+1
      # translate into names
      df$dir_highin <- S[[i]]$output$groups[inds]
    }
    # output
    name <- S[[i]]$output$name
    ws=openxlsx::addWorksheet(wb,sheetName=name)
    openxlsx::writeData(wb=wb, sheet=name, x=df)
  }
  openxlsx::saveWorkbook(wb, file, overwrite=T)

  # add status information
  funargs <- mti_funargs()
  D %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("Exported sheets '%s' to Excel file '%s'",
                       S %>% purrr::map("output") %>% purrr::map("name") %>% unlist() %>% paste0(collapse = ', '), file)
    )

  # pass SummarizedExperiment back, so pipeline can keep running
  D
}
