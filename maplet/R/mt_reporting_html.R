#' Generates markdown-based HTML output from \code{SummarizedExperiment}
#'
#' @description
#' Generates a fully automated report version of an entire (linear) pipeline.
#'
#' @description
#' Will generate RMD markdown document that can be adapted and then knitted to HTML.
#'
#' @param D \code{SummarizedExperiment} input.
#' @param file Output HTML filename.
#' @param title Title of RMD document. Default: 'RMD output'.
#' @param output_calls Output detailed info on function calls? Default: F.
#' @param number_sections Number sections and sub-sections? Default: F.
#' @param start_after Tag name or UUID of pipeline step AFTER which to start. Default: NA (i.e. output entire pipeline).
#' @param use_plotly EXPERIMENTAL. Output interactive plotly plots? WARNING: Setting this argument to TRUE can significantly
#'    increase the runtime of the function. If plotly plots are not needed for a report, we recommend keeping this argument set to
#'    FALSE. Default: F.
#' @param keep_tmp Keep temporary files? Can be used to manually edit RMD afterwards. Default: F.
#'
#' @return Does not change the \code{SummarizedExperiment} object.
#'
#' @author JK
#'
#' @export
mt_reporting_html <- function(D,
                              file,
                              title = 'RMD output',
                              output_calls=F,
                              number_sections=F,
                              start_after=NA,
                              use_plotly=F,
                              keep_tmp=F
) {

  # validate argument
  stopifnot("SummarizedExperiment" %in% class(D))

  # unique string
  ustr <- uuid::UUIDgenerate()

  # define file names
  rmdfile <- sprintf("tmp_%s.RMD", ustr)
  rdsfile <-  sprintf("tmp_%s.rds", ustr) # only used of keep_tmp==T

  # generate RMD
  D %>% reporting_generateMD(
    file = rmdfile,
    readfrom = rdsfile,
    title = title,
    output.calls = output_calls,
    number.sections = number_sections,
    start.after=start_after,
    use.plotly = use_plotly,
    keep_tmp = keep_tmp,
    outfile = file)

  # save temp file that will be input for the RMD?
  if (keep_tmp){
    save(D, file=rdsfile)
    # keep temp files, move them to their final destination name
    file.rename(rmdfile, paste0(tools::file_path_sans_ext(file),'.rmd'))
    file.rename(rdsfile, paste0(tools::file_path_sans_ext(file),'.rds'))
    rmdfile <- paste0(tools::file_path_sans_ext(file),'.rmd')
  }

  # knit
  rmarkdown::render(rmdfile, params=list(D=D))

  # clean up
  if (!keep_tmp) {
    # rename to correct name
    file.rename(paste0(tools::file_path_sans_ext(rmdfile),'.html'), file)
    # no temp files left behind
    file.remove(rmdfile)
    # .rds does not need to be deleted because it was never generated
  }

  # return document, in case pipeline is supposed to keep running
  D
}


#' Markdown-based report generator
#'
#' Generates a fully automated report version of an entire (linear) pipeline.
#'
#' Will generate RMD markdown document that can be adapted and then knitted to HTML.
#'
#' @param D \code{SummarizedExperiment} input
#' @param file File to be generated
#' @param readfrom Name of R data file data will be loaded and is supposed to contain SummarizedExperiment "D". Will not actually be loaded in this function, but while knitting the RMD later.
#' @param title Title of RMD document
#' @param firstheading Name of first heading
#' @param use.plotly Output interactive plotly plots? (experimental)
#' @param output.calls Output full information about all parameters of each function call into RMD?
#' @param number.sections Number sections and sub-sections? (default: F)
#' @param start.after UUID of pipeline step AFTER which to start (default: none, i.e. output entire pipeline)
#' @param keep_tmp Keep the intermediate rmarkdown file? Default: F.
#' @param outfile Name of html file to be created. Required if keep_tmp==T.
#'
#' @returns nothing
#'
#' @importFrom magrittr %>% %<>%
#' @import SummarizedExperiment
#'
#' @author JK, KC
#'
#' @noRd
reporting_generateMD <- function(
  D,
  file = 'MT.RMD',
  readfrom = 'mt.rds',
  title = 'RMD output',
  firstheading='Output',
  use.plotly=F,
  output.calls=F,
  number.sections=F,
  start.after=NA,
  keep_tmp=F,
  outfile
) {


  #### helper functions
  out <- function(str){writeLines(str, h)}

  writechunk <- function(code, params='') {
    if (nchar(params)>0) params=paste0(' ', params)
    out(sprintf("```{r%s}", params))
    out(code)
    out("```\n")
  }


  #### initialize output file
  h <- file(file, open='wt')

  #### markdown header and first heading
  out(glue::glue('
---
title: {title}
output:
  html_document:
    toc: true
    toc_float: TRUE
    {if(number.sections){"number_sections: true"}else{""}}
params:
  D: NA
---

    '))

  #### determine where to start from (first pipeline step, or after one)
  start.from = 1 # by default
  if (!is.na(start.after)) {
    # extract all UUIDs
    allids <- D %>% metadata() %>% .$results %>% purrr::map("uuid") %>% unlist()

    # if tag name provided, get uuid
    tag_name_list <- maplet::mtm_res_get_entries(D, c("reporting", "tag")) %>% purrr::map("output") %>% unlist()
    tag_name_idx <- match(start.after, tag_name_list)
    if(!is.na(tag_name_idx)){
      start.after <- tag_name_list[tag_name_idx] %>% names() %>% gsub("reporting_tag.", "", .)
    }

    # find the one to start after
    start.from <- which(allids==start.after)
    # error-check
    if (length(start.from)==0)
      stop(sprintf("Could not find pipeline step to start after: '%s'", start.after))
    if (D %>% metadata() %>% .$results %>% length() == start.from)
      stop(sprintf("Cannot start after pipeline step '%s', because it's the last entry of the pipeline", start.after))
    start.from <- start.from + 1
  }


  #### global chunk options
  writechunk("# default chunk options\nknitr::opts_chunk$set(warning=F,echo=F,results='hide',message=F)", params = "echo=F")

  #### chunk that loads libraries
  if (!use.plotly) {
    writechunk('# load libraries\nlibrary(maplet)\n')
  } else  {
    writechunk('# load libraries\nlibrary(maplet)\nlibrary("plotly")')
  }

  #### chunk that assigns list r from data
  if(keep_tmp){
    rds_file <- paste0(tools::file_path_sans_ext(outfile),".rds")
    writechunk(glue::glue('load("{rds_file}")\nr <- metadata(D)$results'))
  }else{
    writechunk(glue::glue('r <- metadata(D)$results'))
  }

  #### start of output
  out(glue::glue('# {firstheading}'))

  #### loop over results
  lvl=2 # current level of results = 2

  # loop over results
  r <- metadata(D)$results

  for (i in start.from:length(r)) {
    # to ignore?
    if (length(r[[i]]$fun)<2 || r[[i]]$fun[2]!="void") { # ignore void

      # reporting step?
      if (r[[i]]$fun[1]!="reporting") {
        # not reporting, actual pipeline step

        ## header
        out(glue::glue('{strrep("#",lvl)} {r[[i]]$fun %>% paste(collapse="_")}'))

        ## detailed arguments?
        if (output.calls) {
          L <- r[[i]]$args
          out("*Function arguments:*<br/>")
          out(names(L) %>% lapply(function(x){sprintf("%s=%s",x, toString(L[[x]]))}) %>% paste0(collapse = "<br/>"))
          out("")
        }

        ## log text
        out(glue::glue('*Log text:*<br/>{r[[i]]$logtxt}\n\n'))

        ## plot?
        if (r[[i]]$fun[1]=="plots") {

          # special parameters?
          extraparams <- ""
          if (r[[i]]$fun[2]=="stats"&r[[i]]$fun[3]=="pathway"&r[[i]]$fun[4]=="bar") {
            # dynamic height
            if(r[[i]]$output2$nr!=0){
              # set plot height
              height <- (62+(r[[i]]$output2$nr*r[[i]]$output2$npanrow*23.7143)+84)/2304*12 # manually curated using pixel measurements on example
              width <- 5+3*r[[i]]$output2$npancol
            } else{
              # if empty plot, set height to 3
              height <- 3
              width <- 8
            }
            extraparams <- sprintf(",fig.width=%f,fig.height=%f", width, height)
          } else if(length(r[[i]]$fun)>=3){
            if(r[[i]]$fun[2]=="box"&(r[[i]]$fun[3]=="scatter"|r[[i]]$fun[3]=="special")){
              # dynamic height
              if(!is.null(r[[i]]$output2)){
                # set plot height
                height <- (23+(r[[i]]$output2*160.65)+33)/2304*32 # manually curated using pixel measurements on example
                width <- 7
              } else{
                # if empty plot, set height to 3
                height <- 3
                width <- 7
              }
              extraparams <- sprintf(",fig.width=%f,fig.height=%f", width, height)
            }
          }


          # plot
          if (!use.plotly) {
            # use output2 for mt_plots_net
            if(r[[i]]$fun[2]=="net"){
              writechunk(glue::glue("r[[{i}]]$output2"), params = extraparams)
            }else{
              writechunk( glue::glue("r[[{i}]]$output"), params = extraparams)
            }
          } else {
            writechunk( glue::glue("
plotlist = r[[{i}]]$output %>% lapply(ggplotly)
htmltools::tagList(setNames(plotlist, NULL))
                             "), params=paste0('results="show"', extraparams))
          }
          # }

        }

        ## statistical result table?
        if (r[[i]]$fun[1]=="stats") {
          # write warning if df too large
          if(nrow(r[[i]]$output$table) > 1000){
            out(glue::glue('WARNING: Large data frame ({nrow(r[[i]]$output$table)} rows). Displaying first
                           1000 rows.'))
          }
          # write out datatable
          writechunk(glue::glue('
# extract result table
df<-r[[{i}]]$output$table
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


      } else {

        # data or stats reporting step
        if(r[[i]]$fun[2]=="data" || r[[i]]$fun[2]=="stats"){
          ## header
          out(glue::glue('{strrep("#",lvl)} {r[[i]]$fun %>% paste(collapse="_")}'))
          ## detailed arguments?
          if (output.calls) {
            L <- r[[i]]$args
            out("*Function arguments:*<br/>")
            out(names(L) %>% lapply(function(x){sprintf("%s=%s",x, toString(L[[x]]))}) %>% paste0(collapse = "<br/>"))
            out("")
          }
          ## log text
          out(glue::glue('*Log text:*<br/>{r[[i]]$logtxt}\n\n'))
        }


        # special reporting step
        if (r[[i]]$fun[2]=="heading") {
          # add extra heading
          out(glue::glue('{strrep("#",r[[i]]$output$lvl)} {r[[i]]$output$title}'))
          # result level is this heading +1
          lvl = r[[i]]$output$lvl + 1
        }

        # special reporting step
        if (r[[i]]$fun[2]=="text") {
          # enter text
          out(r[[i]]$output$text)
          out("")
        }
      }
    }
  }

  # clean up
  close(h)

}

