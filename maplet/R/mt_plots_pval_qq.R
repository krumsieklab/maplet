#' Generate p-value qq plot
#'
#' QQ plot against uniform distribution.
#'
#' @param D \code{SummarizedExperiment} input.
#' @param stat_name Name of the statistical result.
#'
#' @return $result$output: plot, p-value qq plot
#'
#' @examples
#' \dontrun{... %>% mt_plots_pval_qq(stat_name='comp') %>% ...  # for one
#' }
#'
#' @import ggplot2
#'
#' @author JK
#'
#' @export
mt_plots_pval_qq <- function(D, stat_name) {

  # validate argument
  stopifnot("SummarizedExperiment" %in% class(D))

  # trick: access argument so that a missing argument error is thrown from here instead of from inside mtm_get_stat_by_name
  stat_name
  # get statistical result
  res <- mtm_get_stat_by_name(D, stat_name)

  # create plot
  p <- gg_qqplot(res$p.value) + ggtitle(sprintf("P-value QQ plot for '%s'", stat_name))

  # add status information & plot
  funargs <- mti_funargs()
  metadata(D)$results %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = glue::glue("P-value QQ plot for {paste0(stat_name, collapse=', ')}"),
      output = list(p)
    )

  # return
  D

}


# Function from https://slowkow.com/notes/ggplot2-qqplot/
#
#' Create a quantile-quantile plot with ggplot2.
#'
#' Assumptions:
#'   - Expected P values are uniformly distributed.
#'   - Confidence intervals assume independence between tests.
#'     We expect deviations past the confidence intervals if the tests are
#'     not independent.
#'     For example, in a genome-wide association study, the genotype at any
#'     position is correlated to nearby positions. Tests of nearby genotypes
#'     will result in similar test statistics.
#'
#' @param ps Vector of p-values.
#' @param ci Size of the confidence interval, 95% by default.
#'
#' @return A ggplot2 plot.
#'
#' @examples
#' \dontrun{gg_qqplot(runif(1e2)) + theme_grey(base_size = 24)
#' }
#'
#' @noRd
gg_qqplot <- function(ps, ci = 0.95) {
  n  <- length(ps)
  df <- data.frame(
    observed = -log10(sort(ps)),
    expected = -log10(stats::ppoints(n)),
    clower   = -log10(stats::qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
    cupper   = -log10(stats::qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
  )
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  ggplot(df) +
    geom_ribbon(
      mapping = aes(x = expected, ymin = clower, ymax = cupper),
      alpha = 0.1
    ) +
    geom_point(aes(expected, observed), shape = 1, size = 3) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    # geom_line(aes(expected, cupper), linetype = 2, size = 0.5) +
    # geom_line(aes(expected, clower), linetype = 2, size = 0.5) +
    xlab(log10Pe) +
    ylab(log10Po)
}

