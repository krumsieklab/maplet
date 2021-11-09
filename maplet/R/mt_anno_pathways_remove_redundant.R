#' Remove redundant pathway annotations
#'
#' Remove identical pathways from an existing \code{SummarizedExperiment} data structure that has a column of pathway annotations.
#'
#' @param D \code{SummarizedExperiment} input.
#' @param feat_col Column containing feature IDs.
#' @param pw_col Column containing pathway IDs.
#'
#' @return rowData: Redundant pathway annotation from SE pw_col column filtered.
#'
#' @examples
#' # first annotate features using smp_db and then remove redundant pathways
#' \dontrun{... %>%
#'   mt_anno_pathways_HMDB(in_col = "HMDb_ID",
#'                          out_col = "smp_db",
#'                          pwdb_name = "SMP",
#'                          db_dir = system.file("extdata", "precalc/hmdb/", package = "maplet")) %>%
#'   mt_anno_pathways_remove_redundant(feat_col = "HMDb_ID", pw_col = "smp_db") %>%
#' ...}
#'
#' @author PG
#'
#' @export
mt_anno_pathways_remove_redundant <- function(D, feat_col, pw_col ) {

  # check arguments
  stopifnot("SummarizedExperiment" %in% class(D))

  if(!feat_col %in% names(rowData(D)))
    stop(glue::glue("feat_col is invalid. please enter an existing column name."))

  if(!pw_col %in% names(rowData(D)))
    stop(glue::glue("pw_col is invalid. please enter an existing column name."))

  row_data <-
    D %>%
    rowData() %>%
    as.data.frame() %>%
    dplyr::select(feat_col = !!feat_col,
                  pw_col = !!pw_col)

  mti_logstatus(glue::glue("creating grouping indices for the pathway IDs in {pw_col} "))
  row_data_indexed <-
    row_data %>%
    dplyr::filter(!is.na(feat_col),
           pw_col != "NULL") %>%
    tidyr::unnest_longer(pw_col) %>%
    dplyr::group_by(pw_col) %>%
    dplyr::arrange(feat_col) %>%
    dplyr::mutate(feat_cols = stringr::str_c(feat_col, collapse = "|")) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(feat_col,
              pw_col,
              pw_idx = dplyr::group_indices(., feat_cols))

  pw_id_replacement <- row_data_indexed %>%
    dplyr::group_by(pw_idx) %>%
    dplyr::arrange(pw_col) %>%
    dplyr::mutate(tmp_ID = dplyr::first(pw_col),
           all_IDs = stringr::str_c(pw_col %>% unique(), collapse = "|")) %>%
    dplyr::ungroup() %>%
    dplyr::select(feat_col,
                  ID = tmp_ID,
                  all_IDs) %>%
    dplyr::distinct() %>%
    dplyr::group_by(feat_col) %>%
    dplyr::summarise(ID = stringr::str_c(ID, collapse = ", ")) %>%
    dplyr::mutate(ID = stringr::str_split(ID, ", "))
  
  # match the nested pathways to our lieblings IDs
  match_idx <-match(row_data$feat_col, pw_id_replacement$feat_col)
  pw_id_replacement <- pw_id_replacement$ID[match_idx]

  pwdb_summary_replacement <-
    dplyr::inner_join(metadata(D)$pathways[[pw_col]],
               dplyr::select(row_data_indexed, -feat_col) %>%
                 dplyr::distinct(),
               by = c("ID" = "pw_col")) %>%
    dplyr::group_by(pw_idx) %>%
    dplyr::arrange(ID) %>%
    dplyr::mutate(tmp_ID = dplyr::first(ID),
           all_IDs = stringr::str_c(ID, collapse = "|"),
           tmp_name = dplyr::first(pathway_name),
           all_pathway_names = stringr::str_c(pathway_name, collapse = "|")) %>%
    dplyr::filter(ID == tmp_ID) %>%
    dplyr::ungroup() %>%
    dplyr::select(-c(tmp_ID, tmp_name, pw_idx))

  # replace rowData and metadata of D
  rowData(D)[[pw_col]] <- pw_id_replacement
  metadata(D)$pathways[[pw_col]] <- pwdb_summary_replacement

  funargs <- mti_funargs()
  D %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf('removed redundant pathway annotations using the %s column', pw_col)
    )

  D
}
