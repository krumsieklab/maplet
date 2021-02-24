#' Write out all pathway annotations, possible redundant pathways, and feature annotations
#'
#' Creates an Excel file that contains detailed feature-to-pathway mapping information. Writes out 4 different forms of the mapping
#' information. Also includes information about which pathways are redundant, i.e. contain the same set of features.
#'
#' @param D \code{SummarizedExperiment} input.
#' @param pw_col Name of the pathway annotation column.
#' @param file Output filename to write to.
#'
#' @return Does not change the \code{SummarizedExperiment} object.
#'
#' @examples
#' \dontrun{... %>%
#' mt_anno_pathways_graphite(in_col = "KEGG",
#'                           out_col = "kegg_db",
#'                           pwdb_species = "hsapiens",
#'                           pwdb_name = "kegg") %>%
#' mt_anno_pathways_remove_redundant(feat_col = "KEGG", pw_col = "kegg_db") %>%
#' mt_write_pathways(pw_col='kegg_db', file='pwannos.xlsx')}
#'
#' @author JK
#'
#' @export
mt_write_pathways <- function(D, pw_col, file) {

  wb <- openxlsx::createWorkbook()
  # helper function
  add_sheet <- function(df, sheetname) {
    openxlsx::addWorksheet(wb, sheetname)
    openxlsx::writeData(wb, sheetname, df, rowNames = F, colNames=T)
  }

  ### pathways with one annotation per row ("long")
  # just generate, don't export yet
  pws.long <- D %>%
    metadata() %>%
    .$pathways %>%
    .[[pw_col]] %>%
    tidyr::separate_rows(all_IDs, all_pathway_names, sep='\\|')

  pws.wide <- D %>%
    metadata() %>%
    .$pathways %>%
    .[[pw_col]]


  ### feature annotations
  mets.long <- D %>%
    # build list of pathway IDs seperated by |
    rowData() %>%
    tibble::as.tibble() %>%
    dplyr::select(dplyr::one_of(c("name",pw_col))) %>%
    dplyr::rename(pw=dplyr::one_of(pw_col)) %>%
    dplyr::mutate(pw_str= sapply(pw, paste, collapse='|')) %>%
    dplyr::select(name, pw_str) %>%
    # explode into one row per annotation
    tidyr::separate_rows(pw_str, sep='\\|') %>%
    # add original pathway names
    dplyr::left_join(pws.wide, by = c("pw_str" = "ID")) %>%
    dplyr::select(name, all_IDs, all_pathway_names) %>%
    # explode redundant ones
    tidyr::separate_rows(all_IDs, all_pathway_names, sep='\\|') %>%
    # clean up
    dplyr::rename_( .dots=stats::setNames(list('all_IDs'), pw_col)) %>%
    dplyr::rename(pathway_name = all_pathway_names)
  # write out sorted by features and by pathway_name
  mets.long %>% dplyr::arrange(name) %>% add_sheet(sheetname = 'fullmap_by_feat')
  mets.long %>% dplyr::arrange(pathway_name) %>% add_sheet(sheetname = 'fullmap_by_pw')

  ### generate wide version of features
  mets.wide <- mets.long %>%
    dplyr::group_by_(pw_col, 'pathway_name') %>%
    dplyr::summarise(features = stringr::str_c(name, collapse = "|")) %>%
    dplyr::select(dplyr::one_of(pw_col,'features'))

  ### export long pathway list
  pws.long %>%
    dplyr::mutate(is.redundant=  ifelse(duplicated(pws.long$pathway_name), 'yes', 'no')) %>%
    dplyr::left_join(mets.wide, by=c("ID"=pw_col)) %>%
    add_sheet(sheetname = 'pathways_long')

  ### export pathway annotations with redundant ones in list form ("wide")
  D %>%
    metadata() %>%
    .$pathways %>%
    .[[pw_col]] %>%
    dplyr::left_join(mets.wide, by=c("ID"=pw_col)) %>%
    add_sheet(sheetname = 'pathways_compact')


  # write out Excel file
  openxlsx::saveWorkbook(wb, file=file, overwrite=TRUE)


  # add status information
  funargs <- mti_funargs()
  metadata(D)$results %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("Pathway annotations in '%s' exported to Excel file '%s'", pw_col, file)
    )

  # pass SummarizedExperiment back, so pipeline can keep running
  D

}
