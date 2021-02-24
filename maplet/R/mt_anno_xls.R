#' Load annotations from Excel file
#'
#' @description
#' Loads annotations and merges them into current SummarizedExperiment.
#' Performs "left-joins", i.e. leaves the original SE unchanged and just adds information where it can be mapped.
#' Can load annotations for both features (rowData) and samples (colData).
#'
#' @description
#' If annotation fields are already existing, this function will fill up any NAs with the values from the new file.
#' Existing values are not overwritten.
#'
#' @param D \code{SummarizedExperiment} input.
#' @param file Name of input Excel file.
#' @param sheet Name or number of sheet.
#' @param anno_type Either "samples" (colData) or "features" (rowData).
#' @param anno_id_col Column in annotation file that contains ID information for mapping.
#' @param data_id_col Column in existing data that contains ID information for mapping. Default: Equal to anno_id_col.
#' @param no_map_err Throw error (T) or warning (F) if something does not map. Default: F.
#' @param replace_names_col OPTIONAL. Name of column from new annotation file to use to overwrite colnames (i.e. sample names) of SE. Default: none.
#'
#' @return rowData or colData: New annotation columns added.
#'
#' @examples
#' \dontrun{
#'   # Load data, two sheets with sample annotations, and one sheet with feature annotations from the same file
#'   D <-
#'   # load raw data
#'   mt_load_xls(file=file, sheet="data", samples_in_rows=T, id_col="SAMPLE_NAME") %>%
#'   # sample annotations from metabolomics run
#'   mt_anno_xls(file=file, sheet="sampleinfo", anno_type="samples", anno_id_col = "SAMPLE_NAME") %>%
#'   # sample annotations from clinical table
#'   mt_anno_xls(file=file, sheet="clinicals", anno_type="samples", anno_id_col="SAMPLE_NAME") %>%
#'   # feature annotations`
#'   mt_anno_xls(file=file, sheet="metinfo", anno_type="features", anno_id_col="BIOCHEMICAL", data_id_col = "name") %>%
#'   ...}
#'
#' @author JK
#'
#' @export
mt_anno_xls <-function(D,
                       file,
                       sheet,
                       anno_type,
                       anno_id_col,
                       data_id_col = anno_id_col,
                       no_map_err = F,
                       replace_names_col = NA) {

  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  if (!(anno_type %in% c("samples","features"))) stop("anno_type must be either 'samples' or 'features'")

  # load excel sheet
  df <- as.data.frame(readxl::read_excel(path=file,sheet=sheet,col_names=T))
  # ensure that sample ID column exists
  if (!(anno_id_col %in% colnames(df))) stop(glue::glue("sample ID column '{anno_id_col}' does not exist in '{basename(file)}, sheet '{sheet}'"))
  if (any(is.na(df[[anno_id_col]]))) stop(glue::glue("sample ID column '{anno_id_col}' contains empty cells, '{basename(file)}, sheet '{sheet}'"))
  rownames(df) <- make.names(df[[anno_id_col]], unique=T)

  # slightly different behavior for samples or features
  if (anno_type=="samples") {
    # ensure the data_id_col column exists
    if (!(data_id_col %in% colnames(colData(D)))) stop(glue::glue("ID column '{data_id_col}' does not exist in current sample annotations of SE"))
    # check that all samples are found in the colnames of the existing dataset
    m <- match(df[[anno_id_col]], colData(D)[[data_id_col]])
    if (any(is.na(m))) {
      msg <- sprintf("The following sample IDs could not be found in the existing data matrix: %s",paste0(df[[data_id_col]][is.na(m)],collapse=","))
      if (no_map_err) stop(msg)
      # else warning(msg)
    }
    # make everything a string
    df[[anno_id_col]] %<>% as.character()
    colData(D)[[data_id_col]] %<>% as.character()
    # merge data frames
    newdf <- coalesce_join(data.frame(colData(D), check.names=F), df, by = setNames(anno_id_col ,data_id_col), join = dplyr::left_join)
    newdf[[anno_id_col]] <- newdf[[data_id_col]] # make sure anno column name also exists (if different from data column name)
    stopifnot(all.equal(newdf[[data_id_col]],colData(D)[[data_id_col]])) # to make sure nothing was mixed up
    rownames(newdf) <- colnames(D)
    colData(D) <- DataFrame(newdf)

    # replace_names_col?
    if (!is.na(replace_names_col)) {
      # copy field, set colnames
      cn = newdf[[replace_names_col]]
      colnames(D) = cn
    }


  } else if (anno_type=="features") {
    # ensure the data_id_col column exists
    if (!(data_id_col %in% colnames(rowData(D)))) stop(glue::glue("ID column '{data_id_col}' does not exist in current feature annotations of SE"))
    # check that all features are found in the $name column of the existing dataset
    m <- match(df[[anno_id_col]], rowData(D)[[data_id_col]])
    if (any(is.na(m))) {
      msg <- sprintf("The following feature IDs could not be found in the existing data matrix: %s",paste0(df[[data_id_col]][is.na(m)],collapse=","))
      if (no_map_err) stop(msg)
      # else warning(msg)
    }
    # make everything a string
    df[[anno_id_col]] %<>% as.character()
    rowData(D)[[data_id_col]] %<>% as.character()
    # merge data frames
    newdf <- coalesce_join(data.frame(rowData(D)), df, by = setNames(anno_id_col,data_id_col), join=dplyr::left_join)
    newdf[[anno_id_col]] <- newdf[[data_id_col]] # make sure anno column name also exists (if different from data column name)
    stopifnot(all.equal(newdf[[data_id_col]],rowData(D)[[data_id_col]])) # to make sure nothing was mixed up
    rownames(newdf) <- rownames(D)
    rowData(D) <- DataFrame(newdf)

    # replace_names_col?
    if (!is.na(replace_names_col)) {
      stop("replace_names_col cannot be used for feature annotations")
    }

  } else
    stop("bug")

  # colData and rowData cannot have overlapping names
  inters <- intersect(colnames(colData(D)), rownames(D))
  if (length(inters)>0) {
    stop(sprintf("There are features (rowData) and sample (colData) variables with the same name: %s", paste0(inters, collapse = ", ")))
  }


  # add status information
  funargs <- mti_funargs()
  metadata(D)$results %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = glue::glue("loaded {anno_type} annotations from Excel file '{basename(file)}, sheet '{sheet}'")
    )

  # return
  D


}


# from https://alistaire.rbind.io/blog/coalescing-joins/
# edited to be able to handle completely disjunct column sets (as in normal left_join). - JK 6/21/20
# edited to make type safe, if two merged columns are incompatible (one numeric, one factor/string) - JK 6/21/20
coalesce_join <- function(x,
                          y,
                          by = NULL,
                          suffix = c(".x", ".y"),
                          join = dplyr::full_join, ...) {
  joined <- join(x, y, by = by, suffix = suffix, ...)
  # names of desired output
  cols <- union(names(x), names(y))
  to_coalesce <- names(joined)[!names(joined) %in% cols]

  if (length(to_coalesce) > 0) {
    suffix_used <- suffix[ifelse(endsWith(to_coalesce, suffix[1]), 1, 2)]
    # remove suffixes and deduplicate
    to_coalesce <- unique(substr(
      to_coalesce,
      1,
      nchar(to_coalesce) - nchar(suffix_used)
    ))
    # coalesce in a type-safe way (if one of them is a character or factor, convert both to character)
    coalesced <- lapply(to_coalesce, function(var) {
      # make sure types are ok
      v1 <- joined[[paste0(var, suffix[1])]]
      v2 <- joined[[paste0(var, suffix[2])]]
      if (is.character(v1) || is.character(v2) || is.factor(v1) || is.factor(v2) ) {
        v1 %<>% as.character()
        v2 %<>% as.character()
      }
      # combine
      dplyr::coalesce(v1, v2)
    }) %>% as.data.frame() %>% dplyr::as_tibble()
    # original code, not type safe
    # coalesced <- purrr::map_dfc(to_coalesce, ~dplyr::coalesce(
    #   joined[[paste0(.x, suffix[1])]],
    #   joined[[paste0(.x, suffix[2])]]
    # ))

    names(coalesced) <- to_coalesce
    cbind(joined, coalesced)[cols] %>% dplyr::as_tibble()
  } else {
    # return unchanged
    joined
  }
}

