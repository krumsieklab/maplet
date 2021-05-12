#' Extract results from entries in metadata by a list of name pieces
#'
#' Extracts metadata()$results entries in a given namespace of arbitrary depth. See examples below.
#'
#' @param D \code{SummarizedExperiment} input
#' @param name_list list of strings that represent pieces of function names
#'
#' @returns List of results from the entries of interest.
#'
#' @examples
#' \dontrun{
#' # get all result entries for plots
#' D %>% mtm_res_get_entries("plots")
#'
#' # get all result entries for stats
#' D %>% mtm_res_get_entries("stats")
#'
#' # check if quotient normalization has been performed
#' q <- D%>% mtm_res_get_entries(c("pre","norm","quot"))
#' if (length(q)==0) stop("No quotient normalization performed.")
#' }
#'
#' @author JK
#' @export
mtm_res_get_entries <- function(D,name_list) {
  # [sorry for the code, but it works]s
  labels <- lapply(metadata(D)$results,'[[', 'fun')
  m <- rep(T,length(labels))
  # exclude non-matches
  for (i in 1:length(name_list)) {
    m <- m & sapply(labels, function(label){
      if (length(label)<i) F
      else {
        if (label[i]!=name_list[i])F
        else T
      }
    })
  }
  # return
  metadata(D)$results[m]
}

#' Return stats output by name
#'
#' Finds a named statistical result from a mt_stats... function and returns the $output$table dataframe.
#'
#' @param D SummarizedExperiment
#' @param name Name of statistical comparison
#' @param fullstruct optional, output entire $output structure, not just the $table inside.
#'
#' @return $output$table dataframe
#' 
#' @export
mtm_get_stat_by_name <- function(D, name, fullstruct=F){
  stopifnot("SummarizedExperiment" %in% class(D))

  if(! ("results" %in% names(metadata(D))))
    stop("no results element found in D")

  stats <- mtm_res_get_entries(D, "stats")

  if(length(stats) == 0)
    stop("no stats element found in D")

  names  <- stats %>% purrr::map_chr(~.x$output$name)
  output <- which(names == name)

  if(length(output) == 0)
    stop("stat element with name ", name, " does not exist")
  if(length(output)  > 1)
    stop("there are multiple stat elements with name ", name)

  if(!fullstruct) {
    output <- stats[[ output ]]$output$table
    if( ! any(c("tibble", "data.frame") %in% class(output)) )
      stop("output of stat result ", stat, " is not a table")

    if( ! ("var" %in% colnames(output)) )
      stop("output of stat result ", name, " does not have 'var' column")
  } else {
    output <- stats[[ output ]]$output
  }
  output
}

#' Extract all plot objects from pipeline
#'
#' Returns just the plots, either retaining the results structure (a nested list) or unlisted.
#'
#' @param D \code{SummarizedExperiment} input.
#' @param unlist Unlist all plots into one long list. Default: T.
#'
#' @return A list (nested or unlisted) containing all plots.
#'
#' @export
mtm_res_get_plots <- function(D,unlist=T){
  l=sapply(mtm_res_get_entries(D,"plots"),'[[','output',simplify=F)
  if(unlist) l <- unlist(l,recursive=F)
  l
}


#' Plot all plots from a pipeline
#'
#' Opens a device (default: PDF), plots all plots, closes device.
#' Works either on list of plots, or on SE
#'
#' @param input List of plots or SummarizedExperiment
#' @param dev Device to plot to (default: PDF)
#' @param ... Further paramaters to be passed to dev() function.
#' 
#' @return None
#'
#' @export
mtm_plot_all_tofile <- function(input, dev=pdf, ...) {
  if ("SummarizedExperiment" %in% class(input)) {
    plots <- mtm_res_get_plots(input)
  } else {
    plots <- input
  }
  dev(...)
  sapply(plots,plot)
  dev.off()
}


#' Extract all UUIDs for each function call
#' 
#' Extracts all of the UUIDs for each function call and returns a data frame with each function call and corresponding UUID.
#' 
#' @param D \code{SummarizedExperiment} input.
#' 
#' @return A data frame with three columns: (1) call_order, (2) function_calls, (3) UUIDs.
#' 
#' @export
mtm_res_get_uuids <- function(D){
  
  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  
  # extract function names and uuids
  fun_calls <- metadata(D)$results %>% names()
  calls <- sapply(fun_calls, function(x){x %>% strsplit("\\.") %>% unlist() %>% extract2(1)})
  uuids <- sapply(fun_calls, function(x){x %>% strsplit("\\.") %>% unlist() %>% extract2(2)})
  
  # combine into data frame
  fun_uuid_df <- data.frame(call_order = 1:length(calls), function_calls = calls, UUIDs = uuids)
  
  fun_uuid_df
}


#' Extract all assay pointers for each function call
#' 
#' Extracts all assay pointers for each function call and returns a data frame with each function call and corresponding 
#' assay pointer. The assay pointer is an identifier used to determine which assay version was used with a particular function
#' call. Will be NULL for all values if the global setting save_all_assays is FALSE.
#'  
#' @param D \code{SummarizedExperiment} input.
#'  
#' @return A data frame with two columns: (1) call_order, (2) function_calls, (3) assay_ptrs.
#'   
#' @export 
mtm_res_get_ptrs <- function(){
  
  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  
  # extract function names
  res <- metadata(D)$results
  fun_calls <- res %>% names()
  calls <- sapply(fun_calls, function(x){x %>% strsplit("\\.") %>% unlist() %>% extract2(1)})
  
  # extract assay pointers
  tmp_assay_ptrs <- sapply(1:length(res), function(i){res[[i]]$assay_ptr})
  # fix NULL values and unlist
  assay_ptrs <- sapply(tmp_assay_ptrs, function(x){ifelse(is.null(x), NA, x)})
  
  # combine into data frame
  assay_ptr_df <- data.frame(call_order = 1:length(calls), function_calls = calls, assay_ptrs=assay_ptrs)
  
  assay_ptr_df
  
}

#' Return a specific assay version given a UUID or tag name
#' 
#' Returns an assay version associated with a specific function call identified by the UUID or a tag name.
#' 
#' @param D \code{SummarizedExperiment} input.
#' @param uuid A UUID associated with a specific function call.
#' @param tag_name Name of a tag given to a section of a maplet pipeline.
#' 
#' @return An assay data frame.
#' 
#' @export
mtm_get_assay_by_id <- function(D, uuid, tag_name){
  
  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  
  if(missing(uuid) & missing(tag_name)) stop("Must provide a value for either uuid or tag name.")
  if(!missing(uuid) & !missing(tag_name)) stop("Only one of these arguments can be passed at a time: uuid, tag_name.")
  
  if(!missing(tag_name)){
    # check tag name exists
    tag_name_list <- maplet::mtm_res_get_entries(D, c("reporting", "tag")) %>% purrr::map("output") %>% unlist()
    tag_name_idx <- match(tag_name, tag_name_list)
    if(is.na(tag_name_idx)) stop(glue::glue("The tag \'{tag_name}\' was not found."))
    
    # get function call
    fun_call <- names(tag_name_list[tag_name_idx])
    
  }else{
    # check uuid exists
    uuid_df <- mtm_res_get_uuids(D)
    if(uuid %in% uuid_df$UUIDs == F) stop(glue::glue("Value for uuid \'{uuid}\' not found in metadata $results list."))
    
    # get function call
    fun_call <- rownames(uuid_df[uuid_df$UUIDs==uuid,])
    
  }
  
  #extract assay version
  assay_ptr <- metadata(D)$results[[fun_call]]$assay_ptr
  assay <- metadata(D)$assays$assay_lst[[assay_ptr]]
  
  assay
  
}

#' Return a specific assay version given an assay pointer
#' 
#' Returns an assay version identified by a specific assay pointer.
#' 
#' @param D \code{SummarizedExperiment} input.
#' @param ptr A string used to identify assay versions.
#' 
#' @return An assay data frame.
#' 
#' @export
mtm_get_assay_by_ptr <- function(D, ptr){
  
  # check assay pointer exists
  all_assay_ptrs <- metadata(D)$assays$assay_lst %>% names()
  if(ptr %in% all_assay_ptrs == F) stop(glue::glue("Value for ptr {ptr} not found in metadata $assays list."))
  
  # extract assay version
  assay <- metadata(D)$assays$assay_lst[[ptr]]
  
  assay
  
}
