#' Fix ggplot env objects
#'
#' Convert all plots saved as large ggplot env objects to images.
#'
#' @param D \code{SummarizedExperiment} input.
#'
#' @return $results[["plots"]]$output: Convert ggplot env object to image.
#'
#' @examples
#' # example of how to run function
#' \dontrun{... %>%
#'   mt_clean_fix_ggplot() %>%
#' ...}
#'
#' @author KC
#'
#' @import ggplot2
#'
#' @export
mt_clean_fix_ggplot <- function(D){

  # validate arguments
  if(!"SummarizedExperiment" %in% class(D)) stop("D is not a SummarizedExperiment object!")

  plot_list <-  maplet::mtm_res_get_entries(D, c("plots")) %>% purrr::map("output")

  # replace each gg environment with fixed plot
  for(i in 1:length(plot_list)){
    plots <- plot_list[[i]]
    new_plots <- list()
    for(j in 1:length(plots)){
      p <- plots[[j]]
      new_p <- fix_ggplot_env(p)
      new_plots[[j]] <- new_p
    }
    # replace old plots with new plots in metadata
    metadata(D)$results[names(plot_list[i])][[1]]$output <- new_plots

  }

  funargs <- mti_funargs()
  D %<>% 
     mti_generate_result(
       funargs = funargs,
       logtxt = glue::glue("Fixed all ggplot objects.")
     )

  # return
  D

}

#' Fixes the problem of exploding environments when saving ggplots to files
#'
#' Magic code by Mustafa (- JK)
#'
#' @param p ggplot env object
#'
#' @return wrapped plot
#'
#'
#' @author MB (JK)
#'
#' @import ggplot2
#'
#' @noRd
fix_ggplot_env <- function(p) {
  # all the environment for quoted variables leads explosion of the object size
  # problem is also not that simple to just find those variables and clean up the
  # respective environment which I did, which did not solve the problem.
  # best solution so far is to wrap plot, get rid of everything else
  local(
    # transformartion of images into blank panel
    # this is updated
    ggplot2::ggplot(data.frame(x = 0:1, y = 0:1), ggplot2::aes_(x = ~x, y = ~y)) +
      ggplot2::geom_blank() +
      ggplot2::scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
      ggplot2::scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
      ggplot2::annotation_custom(gg_grob, xmin = 0, xmax = 1 , ymin = 0, ymax = 1) +
      ggplot2::theme_void(),
    as.environment(list(gg_grob =ggplotGrob(p))) %>% {parent.env(.)=.GlobalEnv;.})
}
