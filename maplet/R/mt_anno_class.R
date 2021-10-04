#' Change the Class Type of Feature or Sample Annotation Columns
#'
#' Uses the lapply/sapply function to convert selected feature or sample annotation columns given an associated list of desired class 
#' types. Column name / class type pairs can be provided one of two ways: (1) provide a named list where the names are the column names
#' and the values are the class type, (2) provide an excel file and sheet with column names in a column called "Column Name" and class
#' types are in a column called "Class". Can only use one method at a time - will crash if values provided for each of the parameters.
#'
#' @param D \code{SummarizedExperiment} input.
#' @param anno_type Either "samples" (colData) or "features" (rowData).
#' @param col_list Named list of column names, from either colData or rowData, with column names as list name and class type as value.
#'    For example, list(samples="factor", age="numeric").
#' @param file Name of excel file containing two columns: "Column Name" for column names and "Class" for class types.
#' @param sheet Name of excel sheet containing column names ("Column Name") and class types ("Class"). Function will crash if incorrect
#'    column names are used. 
#'
#' @return colData or rowData: Changes the class type of the column.
#'
#' @examples
#'  \dontrun{... %>%
#'  # ensure factor for casecontrol variable
#'  mt_anno_apply(anno_type='samples',
#'                col_name='casecontrol',
#'                fun=as.factor) %>%
#'  ...}
#'
#' @author KC
#'
#' @export
mt_anno_class <- function(D, anno_type, col_list, file, sheet){
  
  accepted_class_types = c("factor", "numeric", "character")
  req_col <- c("Column Name", "Class")
  
  # Validate arguments
  if("SummarizedExperiment" %in% class(D)==F) stop("D is not a SummarizedExperiment object.")
  if(missing(anno_type)) stop("The argument 'anno_type' is required.")
  if(!(anno_type %in% c("samples","features"))) stop("Value of argument 'anno_type' must be either 'samples' or 'features'.")
  if(all(missing(col_list), missing(file), missing(sheet))) stop("Values must be provided for parameters: file & sheet OR col_list.")
  if((!missing(file) & missing(sheet)) | (missing(file) & !missing(sheet))) stop("Parameters file and sheet must be supplied together.")
  if(!missing(col_list) & !missing(file)) stop("Values can only be supplied to parameters file & sheet XOR col_list. Only one method can be used at a time.")
  
  # get column names and class types from file
  if(!missing(file) & !missing(sheet)){
    # load excel sheet
    df <- as.data.frame(readxl::read_excel(path=file,sheet=sheet,col_names=T))
    
    # check for missing required columns
    miss_col <- req_col[req_col %in% colnames(df)==F]
    if(length(miss_col) != 0) stop(paste0("Required column names not found:\n", paste0(miss_col, collapse = "\n")))
    
    col_list <- df[,"Class"] %>% as.list()
    names(col_list) <- df[,"Column Name"]
    
  }
  
  # get annotation data frame
  anno_df = if(anno_type == "features"){rowData(D)%>%as.data.frame()}else{colData(D)%>%as.data.frame()}
  
  # check all column names present in df
  miss_col <- names(col_list)[names(col_list) %in% colnames(anno_df)==F]
  if(length(miss_col)!=0) stop(paste0("The following column names are not present in the annotation data frame:\n",paste0(miss_col, collapse = "\n")))
  # check values are accepted class types
  class_types <- col_list %>% unique() %>% unlist()
  if(any(class_types %in% accepted_class_types==F)) stop(paste0("Unrecognized class type(s). Only the following are allowed:\n", paste0(accepted_class_types, "\n")))
  
  # subset anno_df and columns by class types
  anno_df <- anno_df[,names(col_list)]
  fact_cols <- names(col_list[col_list=="factor"])
  num_cols <- names(col_list[col_list=="numeric"])
  char_cols <- names(col_list[col_list=="character"])
  
  # apply desired class type to columns
  anno_df[,fact_cols] = if(length(fact_cols) > 1){lapply(anno_df[,fact_cols], as.factor)}else if(length(fact_cols)==1){sapply(anno_df[,fact_cols], as.factor)}
  anno_df[,num_cols] = if(length(num_cols) > 1){lapply(anno_df[,num_cols], as.numeric)}else if(length(num_cols)==1){sapply(anno_df[,num_cols], as.numeric)}
  anno_df[,char_cols] = if(length(char_cols) > 1){lapply(anno_df[,char_cols], as.character)}else if(length(char_cols)==1){sapply(anno_df[,char_cols], as.character)}
  
  # write back
  if (anno_type=="features") {
    rowData(D)[names(col_list)] <- anno_df
  } else {
    colData(D)[names(col_list)] <- anno_df
  }
  
  ## add status information
  funargs <- mti_funargs()
  D %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = paste(glue::glue("Modified classes of the following {anno_type} annotations columns (column:class):"),
                     paste(paste(names(col_list), unlist(col_list), sep=":"), collapse = ", "))
    )
  ## return
  D
  
}


