#' Merges Two maplet SE objects
#'
#' Merges one maplet SE object provided as argument (D2) into another provided as input SE object (D1). Combines assay, rowData, 
#' and colData data frames and adds metadata of argument SE object (D2) into input SE object (D1). Will only work if sample IDs for
#' each SE are identical.
#'
#' @param D1 First \code{SummarizedExperiment} input to combine into; the one in the pipeline.
#' @param id_col1 Name of the column in the colData of the first SE containing sample IDs.
#' @param D2 Second \code{SummarizedExperiment} input to add.
#' @param id_col2 Name of the column in the colData of the second SE containing sample IDs.
#' @param crash_on_duplicate Crash function if duplicate features are detected? If FALSE, suffices "_1" and "_2" will be added to
#'    the end of the features in D1 and D2, respectively. Default: T.
#' @param sample_anno_from Which SE object to take sample annotations (colData) data frame from. Must be "first" or "second".
#'    Default: "first".
#'
#' @return Merges second SE object with first SE object.
#' @return $result$output: metadata of second SE object.
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
mt_load_merge_se <- function(D1, id_col1, D2, id_col2, crash_on_duplicate = T, sample_anno_from="first"){
  
  ## -- validate arguments -- ##
  
  # validate SE objects, check all required arguments provided 
  stopifnot("SummarizedExperiment" %in% class(D1))
  stopifnot("SummarizedExperiment" %in% class(D2))
  if(missing(id_col1) | missing(id_col2)) stop("Values must be provided for both id_col1 and id_col2.")
  if(id_col1 %in% colnames(colData(D1)) == F) stop(glue::glue("Value for id_col1 \'{id_col1}\' not found in D1 colData."))
  if(id_col2 %in% colnames(colData(D2)) == F) stop(glue::glue("Value for id_col2 \'{id_col2}\' not found in D2 colData."))
  if(sample_anno_from %in% c("first", "second") == F)stop("Value for sample_anno_from must be \"first\" or \"second\".")
  
  # check sample IDs are identical
  id_tmp1 <- D1[[id_col1]]
  id_tmp2 <- D2[[id_col2]]
  if(!identical(id_tmp1[order(id_tmp1)], id_tmp2[order(id_tmp2)])) stop("Sample IDs are not identical between the two SEs.")
  
  # check sample ID column provided is the same as column names
  col_tmp1 <- colnames(D1)
  col_tmp2 <- colnames(D2)
  if(!identical(col_tmp1[order(col_tmp1)],id_tmp1[order(id_tmp1)])) stop("D1 column names do not match Sample ID column.")
  if(!identical(col_tmp2[order(col_tmp2)],id_tmp2[order(id_tmp2)])) stop("D2 column names do not match Sample ID column.")
  
  # check for identical features
  row_tmp1 <- rownames(D1)
  row_tmp2 <- rownames(D2)
  identical_rows <- intersect(row_tmp1, row_tmp2)
  if(length(identical_rows) > 0){
    if(crash_on_duplicate){
      ir_str <- glue::glue_collapse(glue::glue("{identical_rows}"), "\n")
      stop(glue::glue("The following features are duplicated between the SE objects:\n{ir_str}", 
        "\nTo proceed with duplicates, set crash_on_duplicate=F."))
    }else{
      # recreate SEs differentiating assay and rowData rownames from D1 and D2
      dup_idx1 <- rownames(D1) %>% match(identical_rows,.)
      rownames(D1)[dup_idx1] %<>% paste0("_1")
      
      dup_idx2 <- rownames(D2) %>% match(identical_rows,.)
      rownames(D2)[dup_idx2] %<>% paste0("_2")
      
    }
  }
  
  ## -- merge SE objects -- ##
  
  # combine assay
  df1 <- assay(D1)[, order(colnames(D1))]
  df2 <- assay(D2)[, order(colnames(D2))]
  df <- rbind(df1, df2)
  
  # select colData
  if(sample_anno_from=="first"){
    cd <- colData(D1)
  }else{
    cd <- colData(D2)
  }
  
  # combine rowData
  rd1 <- rowData(D1) %>% as.data.frame() 
  rd2 <- rowData(D2) %>% as.data.frame()
  rd <- plyr::rbind.fill(rd1, rd2)
  rownames(rd) <- c(rownames(rd1), rownames(rd2))
  
  # add block_# column to rowData
  # FIND A BETTER WAY TO DO THIS WITH TIDYVERSE
  rd$merge_block <- "block_"
  rd_idx1 <- match(rownames(rd1), rownames(rd))
  rd_idx2 <- match(rownames(rd2), rownames(rd))
  rd[rd_idx1, "merge_block"] %<>% paste0("1")
  rd[rd_idx2, "merge_block"] %<>% paste0("2")
  
  # create new combined SE object
  D <- SummarizedExperiment(assay    = df,
                            colData  = cd,
                            rowData  = rd,
                            metadata = metadata(D1))
  
  ## -- return SE object -- ##
  # add status information
  funargs <- mti_funargs()
  D %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = glue::glue("Merged two SE objects."),
      output = metadata(D2)$results
    )
  
  # return
  D
  
}
