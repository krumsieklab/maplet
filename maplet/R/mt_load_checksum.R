#' Validate MD5 checksum of file
#'
#' If checksum provided, will check that the file checksum is the same and crash if not. If no checksum provided, stores the checksum
#' of the file in the pipeline log.
#'
#' @param D \code{SummarizedExperiment} input. Missing if first step in pipeline.
#' @param file File path. Does not have to be the same as the dataset loaded in the pipeline.
#' @param checksum Checksum to test for.
#'
#' @return If first step in pipeline, creates \code{SummarizedExperiment} object. Otherwise, does not change the
#'    \code{SummarizedExperiment} object.
#'
#' @examples
#' # first call, to get the checksum (will crash, deliberately)
#' \dontrun{... %>% mt_load_checksum(file="input.xlsx", checksum="") %>% ...}
#'
#' # copy-paste the correct ('actual') checksum from the error message into the call:
#' \dontrun{... %>% mt_load_checksum(file="input.xlsx", checksum="688048bd1eb9c771be0eb49548a6f947") %>% ...}
#'
#' @author JK
#'
#' @export
mt_load_checksum <- function(D, file, checksum) {

  # if first step in pipeline, create SE
  if(missing(D)){
    # create an empty SummarizedExperiment object
    D <- SummarizedExperiment()
  }else{
    # validate argument
    stopifnot("SummarizedExperiment" %in% class(D))
  }

  # throw error if file does not exist
  if(!file.exists(file)) stop(sprintf("File does not exist: %s", file))

  # calculate checksum
  md5 <- tools::md5sum(file)[[1]]
  # crash if wrong
  if (!missing(checksum)) {
    if (checksum != md5) {
      stop(sprintf("Wrong checksum for %s, expected: %s, actual: %s\n", file, checksum, md5))
    }
    logtxt <- sprintf("Correct checksum for %s: %s", file, md5)
  } else {
    # no checksum given by user, just show checksum of file
    logtxt <- sprintf("Checksum for %s: %s", file, md5)
  }

  # add status information
  funargs <- mti_funargs()
  D %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = logtxt
    )

  # return
  D


}
