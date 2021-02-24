# Logging control codes. If not explicitely called, MT will perform standard logging to the console.

# helper function that sends an info message to the "mt" logger
# function also returns the string again, so it can be directly used to store log messages as well
mti_logmsg <- function(msg) { logging::loginfo(mti_escape_percent(msg), logger="mt"); msg }
mti_logstatus <- function(msg) { logging::loginfo(mti_escape_percent(msg), logger="mts"); msg }
mti_logwarning <- function(msg) { logging::loginfo(mti_escape_percent(sprintf("WARNING: %s",msg)), logger="mtw"); msg }

# produce sprintf-safe copy of string (create for new loginfo() behaviour, that interprets strings as formats)
mti_escape_percent <- function(txt) gsub('%','%%',txt)

# main function definition
#' Logging Function
#'
#' Put logging function description here
#'
#' @param D define D here
#' @param console logging active or inactive?
#'
#' @return leaves SummarizedExperiment object unchanged
#'
#' @author whoAuthoredIt?
#'
#' @export
mtm_logging <- function(
  D=NULL,    # SummarizedExperiment, simply being passed through
  console=T  # logging active or inactive?
) {

  # set up logging
  logging::logReset()
  if (console) logging::addHandler(writeToConsole, logger = "mt")
  if (console) logging::addHandler(writeToConsole, logger = "mts")
  if (console) logging::addHandler(writeToConsole, logger = "mtw")

  # return unchanged object
  if(!is.null(D)) D

}



