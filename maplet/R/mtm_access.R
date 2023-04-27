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
mtm_res_get_ptrs <- function(D){

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

#' Extract plots from entries in metadata by a 'name argument'
#'
#' Extracts plots from metadata()$results entries if they contain a 'naming argument' (see
#' description of plot_name argument) that is equal to the user provided name.
#'
#' @param D \code{SummarizedExperiment} input.
#' @param plot_name Any name associated with a plot as a value of a function argument. 'Naming
#'    arguments' include stat_name, stat_list, ml_name, stat1, stat2, and title.
#'
#' @returns List of plots for the associated with the user-provided name.
#'
#' @examples
#' \dontrun{
#' # make a volcano plot
#' D %>% mt_plots_volcano(stat_name = "Age super pw",
#'                        x = statistic,
#'                        feat_filter = p.adj < 0.05,
#'                        colour = p.adj < 0.05)
#'
#' # get plot with 'name' "Age super pw"
#' D %>% mtm_res_get_plots_by_name("Age super pw")
#'}
#'
#' @author KC
#' @export
mtm_res_get_plots_by_name <- function(D, plot_name){

  # possible arguments containing plot names
  POSS_NAMES <- c("title", "stat_name", "stat1", "stat2", "ml_name", "stat_list")

  # run checks and extract list of all plots
  stopifnot("SummarizedExperiment" %in% class(D))
  if(missing(plot_name)) stop("A value must be provided for argument 'plot_name'.")
  if(! ("results" %in% names(metadata(D)))) stop("no results element found in D")

  all_plots <- mtm_res_get_entries(D, "plots")
  if(length(all_plots) == 0) stop("No plots element found in D.")

  # extract arguments from those plot results that have 'possible name' arguments
  arg_list <- all_plots %>% purrr::map("args") %>% purrr::map(names) %>%
    sapply(function(x){length(intersect(POSS_NAMES, x)) > 0}) %>%
    unlist() %>% all_plots[.] %>% purrr::map('args')

  # find plot result that matches user provided name
  res_name <- sapply(names(arg_list), function(x){
    y <- unlist(arg_list[[x]])
    pn <- y[which(names(y) %in% POSS_NAMES)]
    if(length(pn)==0) return(NULL)
    if(pn == plot_name) return(x)
    NULL
  }) %>% unlist() %>% unname()

  # extract plot object (if present) and return
  if(is.null(q)) return(list())

  if(length(res_name)==1){
    plot_res <- all_plots[[res_name]]
    plot_type <- plot_res$fun
    if("net" %in% plot_type == F) return(plot_res$output)
    list(output=plot_res$output, output2=plot_res$output2)
  }

  if(length(res_name) > 1){
    plot_list <- lapply(all_plots[res_name], function(x){
      plot_type <- x$fun
      if("net" %in% plot_type==F) return(x$output)
      return(list(output=x$output, output2=x$output2))
    })
    plot_list
  }

}

#' Extract plots from entries in metadata by a type
#'
#' Extracts plots from metadata()$results entries of functions that contain the user-provided
#' substring (via the plot_type argument) in their name. Note the substring can contain any portion
#' of the function name except the prefix 'mt_'. See examples.
#'
#' @param D \code{SummarizedExperiment} input.
#' @param plot_type A substring contained within the name of the function(s) of interest.
#'
#' @returns List of plots for the plot types of interest.
#'
#' @examples
#' \dontrun{
#' # make a volcano and sample boxplots plot
#' D %>% mt_plots_volcano(stat_name = "Age super pw",
#'                         'x = statistic,
#'                         'feat_filter = p.adj < 0.05,
#'                         'colour = p.adj < 0.05) %>%
#'   mt_plots_sample_boxplot(color=Diagnosis, title = "Original", plot_logged = T) %>%
#'   mt_pre_batch_median(batch_col = "BOX.NUMBER") %>%
#'   mt_plots_sample_boxplot(color=Diagnosis, title = "After batch correction", plot_logged = T)
#'
#' # get all volcano plots
#' D %>% mtm_res_get_plots_by_type("volcano")
#'
#' # get all sample boxplots
#' D %>% mtm_res_get_plots_by_type("sample_boxplot")
#'
#' # the following also returns all sample boxplots
#' D %>% mtm_res_get_plots_by_type("mple_boxpl")
#'
#' # this doesn't work
#' D %>%  mtm_res_get_plots_by_type("mt_plots_sample_boxplot")
#'}
#'
#' @author KC
#' @export
mtm_res_get_plots_by_type <- function(D, plot_type){

  # run checks and extract list of all plots
  stopifnot("SummarizedExperiment" %in% class(D))
  if(missing(plot_type)) stop("A value must be provided for argument 'plot_type'.")
  if(! ("results" %in% names(metadata(D))))
    stop("No results element found in D.")

  all_plots <- maplet:::mtm_res_get_entries(D, "plots")
  if(length(all_plots) == 0)
    stop("No plots element found in D.")

  # get result names of plots that contain value of plot_type as a substring
  res_names <- all_plots %>% purrr::map('fun') %>%
    sapply(function(x){paste0(x,collapse = '_')}) %>%
    sapply(function(x){grepl(plot_type,x)}) %>% q[.] %>% names()

  # extract plot objects (if present) and return
  if(length(res_names) > 1){
    plot_list <- lapply(all_plots[res_names], function(x){
      if(grepl("net",plot_type)==F) return(x$output)
      return(list(output=x$output, output2=x$output2))
    })
    return(plot_list)
  }

  if(length(res_names)==1){
    x <- all_plots[[res_names]]
    if(grepl("net",plot_type)==F) return(x$output)
    return(list(output=x$output, output2=x$output2))
  }

  return(list())

}
