#' Add pathway information
#'
#' Modified from mt_anno_pathways_hmdb. Adds custom pathways to the already existing \code{SummarizedExperiment} data taken from
#' KEGGREST.
#'
#' @param D \code{SummarizedExperiment} input.
#' @param in_col Column to use for pathway fetching. The selected column must contain protein Uniprot identifiers.
#' @param out_col A new column name for rowData to output pathway information to.
#' @param raw_db_outfile OPTIONAL. Name of file to export the pathway database to. Must be a string containing the path name with a
#'    .xlsx extension.
#'
#' @return rowData: New pathway annotation column for proteins.
#' @return $results$pathways: A dataframe of pathway information.
#'
#' @examples
#' # annotate proteins using kegg
#' \dontrun{... %>%
#'   mt_anno_pathways_uniprot(in_col = "Uniprot_ID", out_col = "kegg_db") %>%
#' ...}
#'
#' @author RB
#'
#' @importFrom data.table :=
#'
#' @export
mt_anno_pathways_uniprot <- function(D, in_col, out_col, raw_db_outfile) {

  # check arguments
  stopifnot("SummarizedExperiment" %in% class(D))

  if(!in_col %in% names(rowData(D)))
    stop(glue::glue("in_col is invalid. please enter an existing column name."))

  # load look-up table of HMDB to KEGG identifiers
  load(system.file("extdata", "precalc/uniprot/ProteinMapping.rds", package = "maplet"))
  pwdb <- path_uniprot_map
  # create a dataframe that enables the mapping between pathway
  # names and IDs. Included also are num_total, num_measured,
  # num_pw_total, and num_pw_measured (see below for further details)

  # num_total - overall number of proteins in entire database (this will
  # be a redundant, repeating number, identical in every row… but it’s the
  # easiest way to store it right now)
  num_total <- pwdb$uniprot_id %>% unique() %>% length()

  # num_measured - overall number of measured features found in the DB
  num_measured =
    pwdb$uniprot_id %>%
    # for this intersect, similar to the one in pwdb_summary (below),
    # it is assumed that the HMDB IDs in the dataset D is non-redundant.
    # else this number may not be accurate.
    intersect(rowData(D)[[in_col]]) %>%
    length()

  mti_logstatus(glue::glue("summarizing KEGG pathway information"))
  pwdb_summary <- pwdb
  # using methods from data.table reduces runtime by almost 10x as compared
  # to dplyr
  pwdb_summary <-
    data.table::setDT(pwdb_summary)[, `:=`(

      num_total = num_total,
      num_measured = num_measured,

      # num_pw_total - the total number of proteins in that pathway
      num_pw_total = data.table::uniqueN(uniprot_id),
      # num_pw_measured - the number of measured proteins in that pathway
      # (for the M type of analysis on the actual measured background).
      num_pw_measured =
        sum(uniprot_id %in% rowData(D)[[in_col]], na.rm = TRUE)
    ),
    by = path_id] %>%
    # some IDs might have more than two names, however, these will be discarded
    # for now
    unique(by = c("path_id")) %>%
    subset(!is.na(path_id),
           select = -c(uniprot_id))

  mti_logstatus(glue::glue("nesting KEGG pathway IDs"))
  # nest all the pathway IDs given our lieblings input IDs (Uniprot IDs for now)
  pwdb_reduced <-
    pwdb %>%
    dplyr::group_by(uniprot_id) %>%
    dplyr::filter(!is.na(path_id)) %>%
    dplyr::distinct(uniprot_id, path_id) %>%
    tidyr::nest(path_id, .key = path_id) %>%
    dplyr::mutate(path_id=
             path_id %>%
             unlist(recursive = FALSE) %>%
             as.list())
  # match the nested pathways to our lieblings IDs
  match_idx <-
    match(rowData(D)[[in_col]],
          pwdb_reduced$uniprot_id)

  pw_col <- pwdb_reduced$path_id[match_idx]

  # add the pathway IDs into D
  rowData(D)[[out_col]] <- pw_col

  # add pathway map to the metadata of D
  pwdb_summary %<>% dplyr::rename(ID=path_id) # must be ID else the mt_anno_pathways_remove_redundant
  metadata(D)$pathways[[out_col]] <- pwdb_summary

  if(!missing(raw_db_outfile)) {
    openxlsx::write.xlsx(pwdb, raw_db_outfile)
  }


  funargs <- mti_funargs()
  metadata(D)$results %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = 'added pathway annotations using the KEGG pathway database'
    )

  D
}
