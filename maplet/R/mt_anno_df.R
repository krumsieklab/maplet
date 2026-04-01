#' Load annotations from a data frame
#'
#' @description
#' Loads annotations from a data.frame (or tibble) and merges them into the
#' current SummarizedExperiment. Performs "left-joins": leaves the original SE
#' unchanged and only adds information where it can be mapped. Supports both
#' features (rowData) and samples (colData).
#'
#' If annotation fields are already present, this function fills any NAs with
#' values from the new data frame. Existing non-NA values are not overwritten.
#'
#' @param D \code{SummarizedExperiment} input.
#' @param anno_df A data frame (or tibble) containing annotations to merge.
#' @param anno_type Either "samples" (colData) or "features" (rowData).
#' @param anno_id_col Column in \code{anno_df} that contains ID information for mapping.
#' @param data_id_col Column in existing data that contains ID information for mapping. Default: equal to \code{anno_id_col}.
#' @param no_map_err Throw error (TRUE) or warning (FALSE) if something does not map. Default: FALSE.
#' @param replace_names_col OPTIONAL. Name of column from \code{anno_df} (post-prefix) to use to overwrite colnames (i.e., sample names) of SE (samples only). Default: none.
#' @param prefix OPTIONAL. If non-empty, all annotation columns except \code{anno_id_col} are prefixed before mapping.
#'
#' @return The input \code{SummarizedExperiment} with new annotation columns added to rowData or colData.
#'
#' @examples
#' \dontrun{
#'   D %>%
#'     mt_anno_df(anno_df = sample_info, anno_type = "samples", anno_id_col = "SAMPLE_NAME") %>%
#'     mt_anno_df(anno_df = metinfo,     anno_type = "features", anno_id_col = "BIOCHEMICAL", data_id_col = "name")
#' }
#'
#' @author JK (adapted)
#' @import tidyverse
#' @export
mt_anno_df <- function(D,
                       anno_df,
                       anno_type,
                       anno_id_col,
                       data_id_col = anno_id_col,
                       no_map_err = FALSE,
                       replace_names_col = NA,
                       prefix = "") {

  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  if (!(anno_type %in% c("samples","features"))) stop("anno_type must be either 'samples' or 'features'")
  atype_str <- ifelse(anno_type=="samples", "sample", "feature")

  # coerce to plain data.frame
  df <- as.data.frame(anno_df, stringsAsFactors = FALSE)

  # ensure that annotation ID column exists and contains no empty cells
  if (!(anno_id_col %in% colnames(df))) stop(glue::glue("{atype_str} ID column '{anno_id_col}' does not exist in provided data frame"))
  if (any(is.na(df[[anno_id_col]]) | df[[anno_id_col]] == "")) stop(glue::glue("{atype_str} ID column '{anno_id_col}' contains empty/NA values"))
  rownames(df) <- make.names(df[[anno_id_col]], unique = TRUE)

  # optionally prefix all annotation columns except the mapping column
  use_prefix <- !is.null(prefix) && !is.na(prefix) && nzchar(prefix)
  if (use_prefix) {
    cols_to_prefix <- setdiff(colnames(df), anno_id_col)
    colnames(df)[match(cols_to_prefix, colnames(df))] <- paste0(prefix, cols_to_prefix)
  }

  # effective column name for replace_names_col after optional prefixing
  replace_names_col_eff <- replace_names_col
  if (!is.na(replace_names_col_eff) && use_prefix) {
    replace_names_col_eff <- paste0(prefix, replace_names_col_eff)
  }

  # slightly different behavior for samples or features
  if (anno_type == "samples") {
    # ensure the data_id_col column exists
    if (!(data_id_col %in% colnames(colData(D)))) stop(glue::glue("ID column '{data_id_col}' does not exist in current sample annotations of SE"))
    # ensure no columns in annotation df already exist in colData
    anno_col_names <- df %>% dplyr::select(-dplyr::one_of(anno_id_col)) %>% colnames()
    data_col_names <- colData(D) %>% as.data.frame() %>% dplyr::select(-dplyr::one_of(data_id_col)) %>% colnames()
    if (any(anno_col_names %in% data_col_names)) stop("Annotation column names already exist in colData.")
    # check that all samples are found in the existing dataset
    m <- match(df[[anno_id_col]], colData(D)[[data_id_col]])
    if (any(is.na(m))) {
      msg <- sprintf("The following sample IDs could not be found in the existing data matrix: %s",
                     paste0(df[[anno_id_col]][is.na(m)], collapse = ","))
      if (no_map_err) stop(msg)
    }
    # check no annotations are duplicated
    if (any(duplicated(df[!is.na(m), anno_id_col]))) stop(glue::glue("The ID column '{anno_id_col}' of the annotation data frame contains duplicated values."))
    # make everything a string
    df[[anno_id_col]] %<>% as.character()
    colData(D)[[data_id_col]] %<>% as.character()
    # merge data frames (left join) and coalesce
    newdf <- coalesce_join(data.frame(colData(D), check.names = FALSE),
                           df,
                           by   = setNames(anno_id_col, data_id_col),
                           join = dplyr::left_join)
    newdf[[anno_id_col]] <- newdf[[data_id_col]]
    stopifnot(isTRUE(all.equal(newdf[[data_id_col]], colData(D)[[data_id_col]])))
    rownames(newdf) <- colnames(D)
    colData(D) <- S4Vectors::DataFrame(newdf)

    # replace sample names?
    if (!is.na(replace_names_col_eff)) {
      if (!(replace_names_col_eff %in% colnames(newdf))) {
        stop(glue::glue("replace_names_col '{replace_names_col}' (effective: '{replace_names_col_eff}') not found in annotations."))
      }
      cn <- as.vector(newdf[[replace_names_col_eff]])
      colnames(D) <- cn
    }

  } else if (anno_type == "features") {
    # ensure the data_id_col column exists
    if (!(data_id_col %in% colnames(rowData(D)))) stop(glue::glue("ID column '{data_id_col}' does not exist in current feature annotations of SE"))
    # ensure no columns in annotation df already exist in rowData
    anno_col_names <- df %>% dplyr::select(-dplyr::one_of(anno_id_col)) %>% colnames()
    data_col_names <- rowData(D) %>% as.data.frame() %>% dplyr::select(-dplyr::one_of(data_id_col)) %>% colnames()
    if (any(anno_col_names %in% data_col_names)) stop("Annotation column names already exist in rowData.")
    # check that all features are found in the existing dataset
    m <- match(df[[anno_id_col]], rowData(D)[[data_id_col]])
    if (any(is.na(m))) {
      msg <- sprintf("The following feature IDs could not be found in the existing data matrix: %s",
                     paste0(df[[anno_id_col]][is.na(m)], collapse = ","))
      if (no_map_err) stop(msg)
    }
    # check no annotations are duplicated
    if (any(duplicated(df[!is.na(m), anno_id_col]))) stop(glue::glue("The ID column '{anno_id_col}' of the annotation data frame contains duplicated values."))
    # make everything a string
    df[[anno_id_col]] %<>% as.character()
    rowData(D)[[data_id_col]] %<>% as.character()
    # merge data frames (left join) and coalesce
    newdf <- coalesce_join(data.frame(rowData(D), check.names = FALSE),
                           df,
                           by   = setNames(anno_id_col, data_id_col),
                           join = dplyr::left_join)
    newdf[[anno_id_col]] <- newdf[[data_id_col]]
    stopifnot(isTRUE(all.equal(newdf[[data_id_col]], rowData(D)[[data_id_col]])))
    rownames(newdf) <- rownames(D)
    rowData(D) <- S4Vectors::DataFrame(newdf)

    # replace_names_col is not supported for features
    if (!is.na(replace_names_col_eff)) {
      stop("replace_names_col cannot be used for feature annotations")
    }

  } else {
    stop("bug")
  }

  # colData and rowData cannot have overlapping variable names
  inters <- intersect(colnames(colData(D)), colnames(rowData(D)))
  if (length(inters) > 0) {
    stop(sprintf("There are features (rowData) and sample (colData) variables with the same name: %s",
                 paste0(inters, collapse = ", ")))
  }

  # add status information
  funargs <- mti_funargs()
  D %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = glue::glue("loaded {anno_type} annotations from data frame")
    )

  # return
  D
}
