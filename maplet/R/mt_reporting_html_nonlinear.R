#' Generates a fully automated report version of non-linear pipeline
#'
#' Generates markdown-based HTML output for non-linear pipelines from a list of SummarizedExperiment objects.
#'
#' @description
#' Will generate RMD markdown document that can be adapted and then knitted to HTML.
#'
#' @param D_list A list of \code{SummarizedExperiment} objects.
#' @param file Output HTML filename.
#' @param title Title of RMD document. Default: 'Non-Linear RMD output'.
#' @param output_calls Output detailed info on function calls? Default: F.
#' @param keep_tmp Keep temporary files? Can be used to manually edit RMD afterwards. Default: F.
#'
#' @return Does not change the \code{SummarizedExperiment} objects. This does not pass through \code{SummarizedExperiment} objects.
#'
#' @author JK, KC
#'
#' @export
mt_reporting_html_nonlinear <- function(D_list,
                                        file,
                                        title = 'Non-Linear RMD output',
                                        output_calls=F,
                                        keep_tmp=F) {

  # validate argument
  ## D_list
  stopifnot("list" %in% class(D_list))
  check_SE <- sapply(D_list, function(p){"SummarizedExperiment" %in% class(p)})
  if(any(check_SE==FALSE)){
    stop("Argument D_list must be a list of SummarizedExperiment objects.")
  }
  ## file
  if(missing(file)){
    stop("Argument file must be provided.")
  }


  res <- MTResultCollector$new()
  res$addMultiple(D_list)

  # throw error if multiple roots (i.e. disjointed D_list)
  if(length(res$graph_roots()) > 1){
    stop("D_list can not be disjointed! A common root must exist for all D_list!")
  }

  # unique string
  ustr <- uuid::UUIDgenerate()

  # define file names
  rmdfile <- sprintf("tmp_%s.RMD", ustr)
  rdsfile <-  sprintf("tmp_%s.rds", ustr)
  # generate RMD
  res %>% reporting_generateMD_nonLinear(
    file = rmdfile,
    read_from = rdsfile,
    title = title,
    output_calls = output_calls)
  # save temp file that will be input for the RMD
  save(D, file=rdsfile)
  # knit
  rmarkdown::render(rmdfile)
  # rename to correct name
  file.rename(paste0(tools::file_path_sans_ext(rmdfile),'.html'), file)
  # clean up
  if (!keep_tmp) {
    file.remove(rmdfile)
    file.remove(rdsfile)
  }

}


#' Markdown-based report generator for non-linear pipelines
#'
#' Generates a fully automated report version of non-linear pipeline.
#'
#' Will generate RMD markdown document that can be adapted and then knitted to HTML.
#'
#' @param res An object of class MTResultCollector.
#' @param file Name of the file to be generated. Should end in .RMD (default "MT_NL.RMD")
#' @param read_from Name of the R data file (.rds) from which the MTResultCollector data will be read when rendering R Markdown (default: "mt_nl.rds")
#' @param title Title of RMD document (default: "RMD Output")
#' @param first_heading Name of the first heading (default: "Output")
#' @param output_calls Output full information about all parameters of each function call into RMD? (default: F)
#'
#' @return Nothing. \code{SummarizedExperiment} object is not passed through.
#'
#' @author KC
#'
#' @noRd
reporting_generateMD_nonLinear <- function(res,
                                           file = "MT_NL.RMD",
                                           read_from = "mt_nl.rds",
                                           title = "RMD Output",
                                           first_heading = "Output",
                                           output_calls = F){


  # Check arguments
  ## res
  if(!("MTResultCollector" %in% class(res))){
    stop("The argument res should be an object of class MTResultCollector.")
  }

  ## file
  stopifnot(is.character(file))
  if(!endsWith(file, ".RMD")){
    stop("The argument file must end with extension .RMD")
  }

  ## read_from
  stopifnot(is.character(read_from))
  if(!endsWith(read_from, ".rds")){
    stop("The argument read_from must end with extension .rds")
  }

  ## TO DO - REMOVE REPORTING STEPS FROM LIST
  ## if reporting steps present, warn user they will be ignored



  # define helper functions
  out <- function(str){writeLines(str, h)}

  writechunk <- function(code, params='') {
    if (nchar(params)>0) params=paste0(' ', params)
    out(sprintf("```{r%s}", params))
    out(code)
    out("```\n")
  }

  getMDLvl <- function(res, order, lvl=1){
    # create empty level list and add first entry of list
    lvl_list <- vector(mode = "list", length = length(order))
    lvl_list[[1]] <- lvl
    # create empty header list and add first entry of list
    head_list <- vector(mode = "list", length = length(order))
    head_list[[1]] <- FALSE

    for(i in 1:length(order)){
      next_lvl <- lvl_list[[i]]
      node <- order[i]
      next_node <- res$graph_next(node)

      if(length(next_node) == 1){
        lvl_list[[i+1]] <- next_lvl
        head_list[[i+1]] <- FALSE
      }else if(length(next_node > 1)){
        next_lvl <- next_lvl + 1
        nn_idx <- which(order %in% next_node)
        lvl_list[nn_idx] <- next_lvl
        head_list[nn_idx] <- TRUE
      }

    }
    md_lists <- list(levels = lvl_list, headers = head_list)
  }


  # initialize output file
  h <- file(file, open='wt')


  # markdown header and first heading
  out(glue::glue('
                 ---
                 title: {title}
                 output:
                  html_document:
                    toc: true
                    toc_float: true
                    number_sections: true
                 ---

                 '))


  # global chunk options
  writechunk("# default chunk options\nknitr::opts_chunk$set(warning=F,echo=F,results='hide',message=F)", params = "echo=F")

  # chunk that loads data
  writechunk(glue::glue('# load data\nload("{read_from}")'))

  # start of output
  #out(glue::glue('# {first_heading}'))

  # loop over results
  lvl=1 # current level of results = 1

  dfs_order <- igraph::dfs(res$graph, res$graph_roots())$order %>% igraph::as_ids()
  md_lsts <- getMDLvl(res, dfs_order)
  lvl_lst <- md_lsts["levels"] %>% unlist() %>% unname()
  head_lst <- md_lsts["headers"] %>% unlist() %>% unname()

  # make sure result list in correct order
  res$lst <- res$lst[dfs_order]

  for(i in 1:length(dfs_order)){

    lvl <- lvl_lst[i]
    make_head <- head_lst[i]

    if(make_head == TRUE){
      # add BRANCH header
      out(glue::glue('{strrep("#",lvl-1)} {"BRANCH"}'))
    }

    if(res$lst[[i]]$fun[1] != "reporting"){ # this check won't be necessary after remove reporting

      ## header
      out(glue::glue('{strrep("#",lvl)} {res$lst[[i]]$fun %>% paste(collapse="_")}'))

      ## detailed arguments?
      if (output_calls) {
        L <- res$lst[[i]]$args
        out("*Function arguments:*<br/>")
        out(names(L) %>% lapply(function(x){sprintf("%s=%s",x, toString(L[[x]]))}) %>% paste0(collapse = "<br/>"))
        out("")
      }

      ## log text
      out(glue::glue('*Log text:*<br/>{res$lst[[i]]$logtxt}\n\n'))

      ## plot?
      if (res$lst[[i]]$fun[1]=="plots") {
        writechunk( glue::glue("res$lst[[{i}]]$output"))
      }


      ## statistical result table?
      if (res$lst[[i]]$fun[1]=="stats") {
        # write warning if df too large
        if(nrow(res$lst[[i]]$output$table) > 1000){
          out(glue::glue('WARNING: Large data frame ({nrow(res$lst[[i]]$output$table)} rows). Displaying first
                         1000 rows.'))
        }
        # write out datatable
        writechunk(glue::glue('
                              # extract result table
                              df<-res$lst[[{i}]]$output$table
                              # add feature names
                              rd <- rowData(D)
                              df <- cbind(name=as.data.frame(rd)$name[match(df$var, rownames(rd))], df) %>%
                              dplyr::arrange(p.value)
                              # subset large data frames
                              if(nrow(df) > 1000) df <- df[1:1000, ]
                              # output
                              DT::datatable(df, rownames = FALSE, filter = "top", options = list(pageLength = 20, lengthMenu = c(10*(2^(0:3)), nrow(df)), autoWidth = TRUE, width = 1200, dom = "Bitlrp", buttons = c("copy", "csv", "excel", "pdf", "print")), class = "cell-border stripe", extensions = "Buttons")  %>% DT::formatStyle(columns = c(1:ncol(df)), fontSize = "80%", target= "row", lineHeight="80%")'),
                   params = "results='asis'")

      }

      # empty line as spacer
      out("")

    }


  }

  # clean up
  close(h)


}

