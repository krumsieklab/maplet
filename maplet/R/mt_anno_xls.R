#' Load annotations from Excel file
#'
#' @description
#' Loads annotations and merges them into current SummarizedExperiment.
#' Performs "left-joins": leaves the original SE unchanged and adds info where it can be mapped.
#' Works for both features (rowData) and samples (colData).
#'
#' If annotation fields already exist, this function fills NAs from the new file.
#' Existing non-NA values are not overwritten.
#'
#' @param D SummarizedExperiment input.
#' @param file Name of input Excel file.
#' @param sheet Name or number of sheet.
#' @param anno_type Either "samples" (colData) or "features" (rowData).
#' @param anno_id_col Column in annotation file that contains ID information for mapping.
#' @param data_id_col Column in existing data that contains ID information for mapping. Default: equal to anno_id_col.
#' @param no_map_err Throw error (TRUE) or warning (FALSE) if something does not map. Default: FALSE.
#' @param replace_names_col OPTIONAL. Column from new annotation (post-prefix) to overwrite sample names (colnames). Only valid for samples.
#' @param prefix OPTIONAL. If non-empty, all annotation columns except `anno_id_col` are prefixed before mapping.
#'
#' @return SummarizedExperiment with new annotation columns added.
#'
#' @examples
#' \dontrun{
#'   D <-
#'   mt_load_xls(file=file, sheet="data", samples_in_rows=T, id_col="SAMPLE_NAME") %>%
#'   mt_anno_xls(file=file, sheet="sampleinfo", anno_type="samples", anno_id_col = "SAMPLE_NAME") %>%
#'   mt_anno_xls(file=file, sheet="clinicals", anno_type="samples", anno_id_col="SAMPLE_NAME") %>%
#'   mt_anno_xls(file=file, sheet="metinfo", anno_type="features", anno_id_col="BIOCHEMICAL", data_id_col = "name") %>%
#'   ...}
#'
#' @author JK
#'
#' @import tidyverse
#'
#' @export
mt_anno_xls <- function(D,
                        file,
                        sheet,
                        anno_type,
                        anno_id_col,
                        data_id_col = anno_id_col,
                        no_map_err = FALSE,
                        replace_names_col = NA,
                        prefix = "") {

  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  if (!(anno_type %in% c("samples", "features")))
    stop("anno_type must be either 'samples' or 'features'")
  atype_str <- ifelse(anno_type == "samples", "sample", "feature")

  # load excel sheet
  df <- as.data.frame(readxl::read_excel(path = file, sheet = sheet, col_names = TRUE))

  # ensure that annotation ID column exists and contains no empty cells
  if (!(anno_id_col %in% colnames(df)))
    stop(glue::glue("{atype_str} ID column '{anno_id_col}' does not exist in '{basename(file)}', sheet '{sheet}'"))
  if (any(is.na(df[[anno_id_col]])))
    stop(glue::glue("{atype_str} ID column '{anno_id_col}' contains empty cells, '{basename(file)}', sheet '{sheet}'"))

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

  if (anno_type == "samples") {
    # ensure the data_id_col column exists
    if (!(data_id_col %in% colnames(colData(D))))
      stop(glue::glue("ID column '{data_id_col}' does not exist in current sample annotations of SE"))

    # check that all samples are found
    m <- match(df[[anno_id_col]], colData(D)[[data_id_col]])
    if (any(is.na(m))) {
      msg <- sprintf("The following sample IDs could not be found in the existing data: %s",
                     paste0(df[[anno_id_col]][is.na(m)], collapse = ","))
      if (no_map_err) stop(msg)
    }

    # check no annotations are duplicated for mapped rows
    if (any(duplicated(df[!is.na(m), anno_id_col])))
      stop(glue::glue("The ID column '{anno_id_col}' of the annotation data frame contains duplicated values."))

    # make everything a string (for safe joins)
    df[[anno_id_col]] %<>% as.character()
    colData(D)[[data_id_col]] %<>% as.character()

    # merge with NA-coalescing; allow existing columns (do not stop on name collisions)
    newdf <- coalesce_join(
      data.frame(colData(D), check.names = FALSE),
      df,
      by   = setNames(anno_id_col, data_id_col),
      join = dplyr::left_join
    )

    # ensure mapping column also exists under anno_id_col name
    newdf[[anno_id_col]] <- newdf[[data_id_col]]
    stopifnot(all.equal(newdf[[data_id_col]], colData(D)[[data_id_col]]))

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
    if (!(data_id_col %in% colnames(rowData(D))))
      stop(glue::glue("ID column '{data_id_col}' does not exist in current feature annotations of SE"))

    # check that all features are found
    m <- match(df[[anno_id_col]], rowData(D)[[data_id_col]])
    if (any(is.na(m))) {
      msg <- sprintf("The following feature IDs could not be found in the existing data: %s",
                     paste0(df[[anno_id_col]][is.na(m)], collapse = ","))
      if (no_map_err) stop(msg)
    }

    # check no annotations are duplicated for mapped rows
    if (any(duplicated(df[!is.na(m), anno_id_col])))
      stop(glue::glue("The ID column '{anno_id_col}' of the annotation data frame contains duplicated values."))

    # make everything a string (for safe joins)
    df[[anno_id_col]] %<>% as.character()
    rowData(D)[[data_id_col]] %<>% as.character()

    # merge with NA-coalescing; allow existing columns (do not stop on name collisions)
    newdf <- coalesce_join(
      data.frame(rowData(D), check.names = FALSE),
      df,
      by   = setNames(anno_id_col, data_id_col),
      join = dplyr::left_join
    )

    # ensure mapping column also exists under anno_id_col name
    newdf[[anno_id_col]] <- newdf[[data_id_col]]
    stopifnot(all.equal(newdf[[data_id_col]], rowData(D)[[data_id_col]]))

    rownames(newdf) <- rownames(D)
    rowData(D) <- S4Vectors::DataFrame(newdf)

    # replace_names_col is not supported for features
    if (!is.na(replace_names_col_eff)) {
      stop("replace_names_col cannot be used for feature annotations")
    }

  } else {
    stop("bug")
  }

  # colData and rowData must not reuse the same variable names
  inters <- intersect(colnames(colData(D)), colnames(rowData(D)))
  if (length(inters) > 0) {
    stop(sprintf(
      "There are feature (rowData) and sample (colData) variables with the same name: %s",
      paste0(inters, collapse = ", ")
    ))
  }

  # add status information
  funargs <- mti_funargs()
  D %<>% mti_generate_result(
    funargs = funargs,
    logtxt  = glue::glue("loaded {anno_type} annotations from Excel file '{basename(file)}', sheet '{sheet}'")
  )

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

    names(coalesced) <- to_coalesce
    cbind(joined, coalesced)[cols] %>% dplyr::as_tibble()
  } else {
    # return unchanged
    joined
  }
}
