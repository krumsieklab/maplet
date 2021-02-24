#' Sets and outputs global pipeline settings
#'
#' Call without parameters to show current settings. Call with list of parameters to set settings. Will crash for invalid setting names.
#'
#' @param D \code{SummarizedExperiment} input. Missing if first step in pipeline.
#' @param settings List of settings.
#'
#' @return If first step in pipeline, creates \code{SummarizedExperiment}. Otherwise, does not change the
#'    \code{SummarizedExperiment} object.
#'
#' @examples
#' \dontrun{...  %>%
#'   mt_settings(list(ggplot_fix=F)) %>% ... # deactivate ggplot fixing
#'  }
#'
#' @author JK
#'
#' @export
mt_settings <- function(D, settings) {

  # if first step in pipeline, create SE
  if(missing(D)){
    # create an empty SummarizedExperiment object
    D <- SummarizedExperiment()
  }else{
    # validate argument
    stopifnot("SummarizedExperiment" %in% class(D))
  }

  # make sure that settings exist in metadata of this pipeline
  D %<>% mti_ensure_settings()

  # check if anything needs to be set
  if (!missing(settings)) {
    # get settings list to work with inside this function
    lst <- mti_settings_list()
    # verify that all provided settings are valid
    # if yes, set them
    for (sname in names(settings)) {
      # valid name
      if (!(sname %in% names(lst))) stop(glue::glue("Invalid setting: '{sname}'. Must be one of: {paste0(names(lst),collapse=', ')}"))
      # valid type
      if (!(lst[[sname]]$class %in% class(settings[[sname]]))){stop(glue::glue("Invalid type for parameter '{sname}'. Expected: {lst[[sname]]$class}. Actual: {paste0(class(settings[[sname]]),collapse=', ')}"))}
      # set
      metadata(D)$settings[[sname]] <- settings[[sname]]
    }
  }

  # build log message with all settings
  s <- metadata(D)$settings
  log_changes <- names(settings) %>% purrr::map(~paste0(.,"=",settings[.])) %>% unlist() %>% paste0(collapse = ", ")
  log_final <- names(s) %>% purrr::map(~paste0(.,"=",s[.])) %>% unlist() %>% paste0(collapse = ", ")

  # add status information & plot
  funargs <- maplet:::mti_funargs()
  metadata(D)$results %<>%
    maplet:::mti_generate_result(
      funargs = funargs,
      logtxt = glue::glue("Changes to settings: {log_changes}. Final setting list: {log_final}")
    )

  # return
  D

}

# list of settings, including types and default parameters

# Return default list
# -> This is where MT developers define the parameters
mti_settings_list <- function() {
  list(
    dummy = list(class="numeric", default=5)
  )
}

# Ensures that a pipeline has settings stored in its metadata. Will initialize with default list if not.
# -> Only called internally in this script.
mti_ensure_settings <- function(D) {
  # validate argument
  stopifnot("SummarizedExperiment" %in% class(D))
  # check if settings missing
  if (!("settings" %in% names(metadata(D)))) {
    # build list
    fulllist <- mti_settings_list()
    lst <- names(fulllist) %>% lapply(function(s){fulllist[[s]]$default})
    names(lst) <- names(fulllist)
    # set list
    metadata(D)$settings <- lst
  }
  # return
  D
}

# Retrieve setting.
# Will use default list of no settings have been made so far.
# Important: Does not add a default list to the SummarizedExperiment (i.e., does not affect the SE), but simply access the default argument.
# -> This function should be called internally by all mt_ functions to check a parameter setting.
mti_get_setting <-  function(D, sname) {
  # validate argument
  stopifnot("SummarizedExperiment" %in% class(D))
  # check if settings missing
  if (!("settings" %in% names(metadata(D)))) {
    # no, use default
    # build list
    fulllist <- mti_settings_list()
    lst <- names(fulllist) %>% lapply(function(s){fulllist[[s]]$default})
    names(lst) <- names(fulllist)
  } else {
    # yes, take from metadata
    lst <- metadata(D)$settings
  }
  # error control
  if (!(sname %in% names(lst))){stop(glue::glue("Invalid pipeline setting name: '{sname}'"))}
  # return
  lst[[sname]]
}
