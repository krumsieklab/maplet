#' Add a Reference Tag
#'
#' Add a reference tag to a section in the pipeline. For use by mt_reporting_html to indicate where in the pipeline the report should
#' start from.
#'
#' @param D \code{SummarizedExperiment} input.
#' @param tag_name Name of the reference tag.
#'
#' @return $results$output: Add reference tag.
#'
#' @examples
#' # add a tag to a section of the pipeline
#' \dontrun{
#'    # add tag to a spot in the pipeline
#'    mt_reporting_tag(tag_name="preprocessing") %>% ...
#'    D %>% mt_reporting_html(file="example.html", start_after="preprocessing")
#' }
#'
#' @author KC
#'
#' @export
mt_reporting_tag <- function(D, tag_name){

  # valide arguments
  if(!"SummarizedExperiment" %in% class(D)) stop("D is not a SummarizedExperiment object!")
  if(missing(tag_name)) stop("A value for tag_name must be provided!")

  # add status information & plot
  funargs <- mti_funargs()
  D %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = glue::glue("Added reference tag: {tag_name}."),
      output = tag_name
    )

  # return
  D

}
