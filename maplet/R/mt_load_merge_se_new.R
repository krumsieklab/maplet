#' Merges Two maplet SE objects
#'
#' Merges one maplet SE object provided as argument (D1) with another provided as input SE object
#' (D2). Combines assay, rowData, and colData data frames. Creates new metadata result list or
#' can continue from one of the existing lists from either SE input (D1 or D2).
#' This function can handle cases of unique, identical or overlapping samples and features. It also
#' allows the user to ensure the uniqueness of the sample and feature IDs from one or both of the
#' SE objects by adding suffixes. In cases of overlapping samples and features, the merge function
#' will attempt to combine assay data frames unless any of the sample-feature pairs contains
#' non-missing values for both D1 and D2 assays.
#' We recommend the user ensure that the colData and rowData annotation column names are
#' identical between D1 and D2. If this is not the case, annotation columns from each SE object
#' will be suffixed with d1 or d2, depending on the SE object of origin, and combined into one
#' large annotation data frame. This may create issues of annotation redundancy and the user must
#' ensure this will not create downstream issues. In the case of overlapping samples or features,
#' annotation columns will be suffixed with d1 and d2 and combined into one data frame as described
#' above to avoid the potential of overwritting annotaiton values.
#'
#' @param D1 First \code{SummarizedExperiment} object to combine.
#' @param D2 Second \code{SummarizedExperiment} object to combine.
#' @param samp_id_col1 From the first SE, name of the column in colData containing sample IDs.
#' @param samp_id_col2 From the second SE, name of the column in colData containing sample IDs.
#' @param feat_id_col1 From the first SE, name of the column in rowData containing feature IDs.
#' @param feat_id_col2 From the second SE, name of the column in rowData containing feature IDs.
#' @param samp_suffix1 Suffix to add to sample names in D1 in order to ensure uniqueness.
#' @param samp_suffix2 Suffix to add to sample names in D2 in order to ensure uniqueness.
#' @param feat_suffix1 Suffix to add to feature names in D1 in order to ensure uniqueness.
#' @param feat_suffix2 Suffix to add to feature names in D2 in order to ensure uniqueness.
#' @param anno_cols_identical Flag indicating whether to force annotation columns for both rowData
#'    and colData to be identical. If TRUE, function will crash if annotation column names are not
#'    identical. If FALSE, function will keep both common and unique columns between the two SE
#'    objects and warn the user to ensure that this will not create any downstream problems. In
#'    cases where one or more sample / features are identical, this flag will be ignored and all
#'    columns will be combined and kept separately. Default: TRUE.
#' @param results_from Which SE object to take result list. Must be one of: "none", "first", or "second".
#'    Default: "none" (new pipeline).
#'
#' @return Merges assay, rowData, and colData of two SE objects.
#' @return $result$output: metadata of first, second or both SE objects.
#'
#' @examples
#' \dontrun{D1 <-
#'   # merge D2 into D1
#'   mt_load_merge_se(id_col1="Sample_ID", D2=D2, id_col2="Sample_ID") %>%
#'   ...}
#'
#' @author KC
#'
#' @export
mt_load_merge_se_new <- function(D1,
                                 D2,
                                 samp_id_col1,
                                 samp_id_col2,
                                 feat_id_col1,
                                 feat_id_col2,
                                 samp_suffix1,
                                 samp_suffix2,
                                 feat_suffix1,
                                 feat_suffix2,
                                 anno_cols_identical = TRUE,
                                 results_from = c('none', 'first', 'second')
                                 ){

  # NTS: There are inconsistencies in the way we enforce the presence of rownames /
  #      colnames in the assay objects. We need to standardize this

  # -- Get arguments and extract required data frames --
  results_from = match.arg(results_from)

  # make data frame because tibble creates issues with rownames
  assay1 <- assay(D1) %>% tibble::rownames_to_column("Feature_ID") %>% as.data.frame()
  rownames(assay1) <- assay1$Feature_ID
  assay1 %<>% dplyr::select(-Feature_ID)
  assay2 <- assay(D2) %>% tibble::rownames_to_column("Feature_ID") %>% as.data.frame()
  rownames(assay2) <- assay2$Feature_ID
  assay2 %<>% dplyr::select(-Feature_ID)

  rd1 <- rowData(D1) %>% as.data.frame()
  rd2 <- rowData(D2) %>% as.data.frame()

  cd1 <- colData(D1) %>% as.data.frame()
  cd2 <- colData(D2) %>% as.data.frame()

  # -- Validate Summarized Experiments and all required arguments --
  stopifnot("SummarizedExperiment" %in% class(D1))
  stopifnot("SummarizedExperiment" %in% class(D2))
  if(any(missing(samp_id_col1), missing(samp_id_col2), missing(feat_id_col1), missing(feat_id_col2)))
    stop("One or more of the required ID column arguments is missing.")
  # fix this later
  stopifnot(samp_id_col1 %in% colnames(cd1))
  stopifnot(samp_id_col2 %in% colnames(cd2))
  stopifnot(feat_id_col1 %in% colnames(rd1))
  stopifnot(feat_id_col2 %in% colnames(rd2))

  # check that rownames match feature IDs
  feat1 <- sort(rd1[,feat_id_col1])
  feat2 <- sort(rd2[,feat_id_col2])

  stopifnot(identical(feat1, sort(rownames(assay1))))
  stopifnot(identical(feat2, sort(rownames(assay2))))

  # check that colnames match sample IDs
  samp1 <- sort(cd1[,samp_id_col1])
  samp2 <- sort(cd2[,samp_id_col2])

  stopifnot(identical(samp1, sort(colnames(assay1))))
  stopifnot(identical(samp2, sort(colnames(assay2))))

  # throw warning if zeros in assays (zeros should be NAs, if not this can cause issues with
  #  overlapping sample-feature pairs)
  if(any(assay1==0, na.rm = T)) warning("Zeros found in assay1. Zeros should be converted to NAs.")
  if(any(assay2==0, na.rm = T)) warning("Zeros found in assay2. Zeros should be converted to NAs.")

  # ---- merge SE objects ----
  # add prefixes sample and feature names (if provided)
  if(!missing(samp_suffix1)){
    colnames(assay1) <- paste(colnames(assay1), samp_suffix1, sep="_")
    cd1 %<>% dplyr::mutate(orig_samp_id = !!as.name(samp_id_col1)) %>%
      dplyr::mutate(!!as.name(samp_id_col1) := stringr::str_c(!!as.name(samp_id_col1), samp_suffix1, sep = "_")) %>%
      tibble::remove_rownames()
    samp1 <- colnames(assay1)
  }
  if(!missing(feat_suffix1)){
    rownames(assay1) <- paste(rownames(assay1), feat_suffix1, sep="_")
    rd1 %<>% dplyr::mutate(orig_feat_name = !!as.name(feat_id_col1)) %>%
      dplyr::mutate(!!as.name(feat_id_col1) := stringr::str_c(!!as.name(feat_id_col1), feat_suffix1, sep = "_")) %>%
      tibble::remove_rownames()
    feat1 <- rownames(assay1)
  }
  if(!missing(samp_suffix2)){
    colnames(assay2) <- paste(colnames(assay2), samp_suffix2, sep="_")
    cd2 %<>% dplyr::mutate(orig_samp_id = !!as.name(samp_id_col2)) %>%
      dplyr::mutate(!!as.name(samp_id_col2) := stringr::str_c(!!as.name(samp_id_col2), samp_suffix2, sep = "_")) %>%
      tibble::remove_rownames()
    samp2 <- colnames(assay2)
  }
  if(!missing(feat_suffix2)){
    rownames(assay2) <- paste(rownames(assay2), feat_suffix2, sep="_")
    rd2 %<>% dplyr::mutate(orig_feat_name = !!as.name(feat_id_col2)) %>%
      dplyr::mutate(!!as.name(feat_id_col2) := stringr::str_c(!!as.name(feat_id_col2), feat_suffix2, sep = "_")) %>%
      tibble::remove_rownames()
    feat2 <- rownames(assay2)
  }

  # check that no common sample-feature pairs have overlapping non-missing values
  common_samples <- intersect(samp1, samp2)
  common_features <- intersect(feat1, feat2)

  if(length(common_samples) > 0 & length(common_features) > 0){
    assay1 %<>% tibble::rownames_to_column("Feature_ID")
    assay2 %<>% tibble::rownames_to_column("Feature_ID")
    common_df1 <- assay1[assay1$Feature_ID %in% common_features, c("Feature_ID", common_samples)]
    common_df2 <- assay2[assay2$Feature_ID %in% common_features, c("Feature_ID", common_samples)]

    tmp1 <- tidyr::pivot_longer(common_df1, common_samples)
    tmp2 <- tidyr::pivot_longer(common_df2, common_samples)
    colnames(tmp1)[3] <- "first"
    colnames(tmp2)[3] <- "second"

    # check at least one value in overlapping pairs is NA, crash if any are not
    tmp <- plyr::join(tmp1, tmp2, by = c("Feature_ID", "name"))
    tmp %<>% dplyr::mutate(test = (is.na(first) | is.na(second)))
    overlap_pairs <- sum(tmp$test==F)
    rm(list= c("tmp", "tmp1", "tmp2"))
    if(overlap_pairs != 0){
      stop(glue::glue("There are {overlap_pairs} overlapping sample-feature pairs both with non-missing values.
                  Check your data or add prefixes to make sample-feature pairs unique."))
    }

    assay1 %<>% tibble::column_to_rownames("Feature_ID")
    assay2 %<>% tibble::column_to_rownames("Feature_ID")

  }

  # add new unique samples from first SE to second SE assay (columns)
  unique_samp1 <- samp1[(samp1 %in% samp2) == FALSE]
  if(length(unique_samp1) != 0) assay2[,unique_samp1] <- NA
  # add new unique features from first SE to second SE assay (rows)
  unique_feat1 <- feat1[(feat1 %in% feat2) == FALSE]
  if(length(unique_feat1) != 0){
    tmp2 <- data.frame(matrix(ncol = ncol(assay2), nrow = length(unique_feat1))) %>% setNames(colnames(assay2))
    rownames(tmp2) <- unique_feat1
    assay2 <- rbind(assay2, tmp2)
  }

  # add new unique samples from second SE to first SE assay (columns)
  unique_samp2 <- samp2[(samp2 %in% samp1) == FALSE]
  if(length(unique_samp2) != 0) assay1[,unique_samp2] <- NA
  # add new unique features from second SE to first SE assay (rows)
  unique_feat2 <- feat2[(feat2 %in% feat1) == FALSE]
  if(length(unique_feat2) != 0){
    tmp1 <- data.frame(matrix(ncol = ncol(assay1), nrow = length(unique_feat2))) %>% setNames(colnames(assay1))
    rownames(tmp1) <- unique_feat2
    assay1 <- rbind(assay1, tmp1)
  }

  # combine extended data frames
  x <- assay1 %>% tibble::rownames_to_column("Feature_ID")
  y <- assay2 %>% tibble::rownames_to_column("Feature_ID")
  joined <- merge_coalesce_join(x, y, by = "Feature_ID")

  joined %<>% tibble::column_to_rownames("Feature_ID")

  # currently no good solution for combining annotations for duplicate samples
  #  for now, if any duplicates found, keep cd1 and cd2 annotation columns separate
  if(length(common_samples) != 0){

    # add d1 and d2 suffixes to end of each cd1 and cd2 anno columns, respectively (exluding sample ID columns)
    cd1 %<>% dplyr::rename_at(dplyr::vars(-!!as.name(samp_id_col1)),function(x) paste0(x,"_d1"))
    cd2 %<>% dplyr::rename_at(dplyr::vars(-!!as.name(samp_id_col2)), function(x) paste0(x,"_d2"))

    unique_cd1_cols <- cd1 %>% dplyr::select(-!!as.name(samp_id_col1)) %>% colnames()
    unique_cd2_cols <- cd2 %>% dplyr::select(-!!as.name(samp_id_col2)) %>% colnames()

    # subset common samples from cd1 and cd2
    common_cd1 <- cd1 %>% dplyr::filter(!!as.name(samp_id_col1) %in% common_samples)
    common_cd2 <- cd2 %>% dplyr::filter(!!as.name(samp_id_col2) %in% common_samples)

    # merge common samples to combine all columns
    cd <- merge(common_cd1, common_cd2, by.x=samp_id_col1, by.y=samp_id_col2)

    # check if unique samples from cd1
    if(length(unique_samp1) != 0){
      #  if TRUE, add empty cd2 anno columns
      unique_cd1 <- cd1 %>% dplyr::filter(!!as.name(samp_id_col1) %in% unique_samp1)
      unique_cd1[,unique_cd2_cols] <- NA

      #  ensure column order matches common sample df, and rbind unique cd1 sample df
      unique_cd1 <- unique_cd1[,match(colnames(cd), colnames(unique_cd1))]
      cd <- rbind(cd, unique_cd1)
    }
    # check if unique samples from cd2
    if(length(unique_samp2) != 0){
      #  if TRUE, add empty cd1 anno columns
      unique_cd2 <- cd2 %>% dplyr::filter(!!as.name(samp_id_col2) %in% unique_samp2)
      unique_cd2[,unique_cd1_cols] <- NA

      #  ensure column order matches common sample df, and rbind unique cd2 sample df
      unique_cd2 <- unique_cd2[,match(colnames(cd), colnames(unique_cd2))]
      cd <- rbind(cd, unique_cd2)
    }

  }else{
    # check colData columns identical
    unique_cd1_cols <- colnames(cd1)[colnames(cd1) %in% colnames(cd2)==F]
    unique_cd2_cols <- colnames(cd2)[colnames(cd2) %in% colnames(cd1)==F]

    if(length(unique_cd1_cols) != 0 | length(unique_cd2_cols) != 0){

      if(anno_cols_identical) stop("The flag anno_cols_identical is set to TRUE. Annotation columns are not identical.")

      warning("Sample annotation (colData) column names are not identical. This may create downstream
            problems.")

      # add unique cd2 columns to cd1
      if(length(unique_cd2_cols) != 0) cd1[,unique_cd2_cols] <- NA

      # add unique cd1 columns to cd2
      if(length(unique_cd1_cols) != 0) cd2[,unique_cd1_cols] <- NA

    }

    # combine colData data frames
    cd <- rbind(cd1, cd2)
  }


  # currently no good solution for combining annotations for duplicate features
  #  for now, if any duplicates found, keep rd1 and rd2 annotation columns separate
  if(length(common_features) != 0){

    # add d1 and d2 suffixes to end of each cd1 and cd2 anno columns, respectively (exluding sample ID columns)
    rd1 %<>% dplyr::rename_at(dplyr::vars(-!!as.name(feat_id_col1)),function(x) paste0(x,"_d1"))
    rd2 %<>% dplyr::rename_at(dplyr::vars(-!!as.name(feat_id_col2)), function(x) paste0(x,"_d2"))

    unique_rd1_cols <- rd1 %>% dplyr::select(-!!as.name(feat_id_col1)) %>% colnames()
    unique_rd2_cols <- rd2 %>% dplyr::select(-!!as.name(feat_id_col2)) %>% colnames()

    # subset common samples from cd1 and cd2
    common_rd1 <- rd1 %>% dplyr::filter(!!as.name(feat_id_col1) %in% common_features)
    common_rd2 <- rd2 %>% dplyr::filter(!!as.name(feat_id_col2) %in% common_features)

    # merge common samples to combine all columns
    rd <- merge(common_rd1, common_rd2, by.x=feat_id_col1, by.y=feat_id_col2)

    # check if unique samples from rd1
    if(length(unique_feat1) != 0){
      #  if TRUE, add empty rd2 anno columns
      unique_rd1 <- rd1 %>% dplyr::filter(!!as.name(feat_id_col1) %in% unique_feat1)
      unique_rd1[,unique_rd2_cols] <- NA

      #  ensure column order matches common sample df, and rbind unique rd1 sample df
      unique_rd1 <- unique_rd1[,match(colnames(rd), colnames(unique_rd1))]
      rd <- rbind(rd, unique_rd1)
    }
    # check if unique samples from rd2
    if(length(unique_feat2) != 0){
      #  if TRUE, add empty rd1 anno columns
      unique_rd2 <- rd2 %>% dplyr::filter(!!as.name(feat_id_col2) %in% unique_feat2)
      unique_rd2[,unique_rd1_cols] <- NA

      #  ensure column order matches common sample df, and rbind unique rd2 sample df
      unique_rd2 <- unique_rd2[,match(colnames(rd), colnames(unique_rd2))]
      rd <- rbind(rd, unique_rd2)
    }

  }else{
    # check rowData columns identical
    unique_rd1_cols <- colnames(rd1)[colnames(rd1) %in% colnames(rd2)==F]
    unique_rd2_cols <- colnames(rd2)[colnames(rd2) %in% colnames(rd1)==F]

    if(unique_rd1_cols != 0 | unique_rd2_cols != 0){

      if(anno_cols_identical) stop("The flag anno_cols_identical is set to TRUE. Annotation columns are not identical.")

      warning("Feature annotation (rowData) column names are not identical. This may create downstream
            problems. Columns will be combined.")

      # add unique cd2 columns to cd1
      if(length(unique_rd2_cols) != 0) rd1[,unique_rd2_cols] <- NA

      # add unique cd1 columns to cd2
      if(length(unique_rd1_cols) != 0) rd2[,unique_rd1_cols] <- NA
    }

    # combine rowData data frames
    rd <- rbind(rd1, rd2)
  }

  samp_id_col <- samp_id_col1
  feat_id_col <- feat_id_col1

  # order assay columns and rows
  joined
  # rearrange cd rows to match assay column order
  rownames(cd) <- cd[,samp_id_col]
  View(cd[match(colnames(joined), cd[,samp_id_col]),])

  # rearrange rd rows to match assay row order
  rownames(rd) <- rd[,feat_id_col]
  View(rd[match(rownames(joined), rd[,feat_id_col]),])

  if(results_from == 'none'){
    results <- list()
    output <- list(d1_results = metadata(D1)$results, d2_results = metadata(D2)$results)
  }else if(results_from == 'first'){
    results <- metadata(D1)$results
    output <- list(d2_results = metadata(D2)$results)
  }else if(results_from == 'second'){
    results <- metadata(D2)$results
    output <- list(d1_results = metadata(D2)$results)
  }

  D <- SummarizedExperiment::SummarizedExperiment(assay = joined,
                                                  rowData = rd,
                                                  colData = cd,
                                                  metadata = list(sessionInfo=utils::sessionInfo(),
                                                                  results=results))

  # add status information
  funargs <- mti_funargs()
  D %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = "Merged two SE objects.",
      output = output
    )

  # return
  D

}


merge_coalesce_join <- function(x, y,
                          by = NULL, suffix = c(".x", ".y"),
                          join = dplyr::full_join, ...) {
  joined <- join(x, y, by = by, suffix = suffix, ...)
  # names of desired output
  cols <- union(names(x), names(y))

  to_coalesce <- names(joined)[!names(joined) %in% cols]
  suffix_used <- suffix[ifelse(endsWith(to_coalesce, suffix[1]), 1, 2)]
  # remove suffixes and deduplicate
  to_coalesce <- unique(substr(
    to_coalesce,
    1,
    nchar(to_coalesce) - nchar(suffix_used)
  ))

  coalesced <- purrr::map_dfc(to_coalesce, ~dplyr::coalesce(
    joined[[paste0(.x, suffix[1])]],
    joined[[paste0(.x, suffix[2])]]
  ))
  names(coalesced) <- to_coalesce

  dplyr::bind_cols(joined, coalesced)[cols]
}
