#' Output information about statistical results
#'
#' Leaves a log entry starting with:
#'   N (outcome): <group1>: <n>, <group2>: <n>, ...   # or "not applicable" if numeric outcome
#'   N (total): <used> of <all> (dropped: <all - used>)
#' Then appends extras (e.g., per-outcome used/unused counts, group-N summary stats,
#' and optionally how many result rows match a filter).
#'
#' Binary numeric outcomes (<= 2 unique used values) are treated as categorical.
#'
#' @param D \code{SummarizedExperiment} input.
#' @param stat_name Name of the statistical comparison.
#' @param stat_filter Filter formula to display number of results, e.g. p.adj<0.2.
#'
#' @return Does not change the \code{SummarizedExperiment} object (beyond adding a log entry).
#'
#' @examples
#' \dontrun{...  %>%
#'   mt_reporting_stats(stat_name="comp1", stat_filter=p.adj<0.2) %>% ...}
#'
#' @author JK
#'
#' @export
mt_reporting_stats <- function(D, stat_name, stat_filter) {

  # validate argument
  stopifnot("SummarizedExperiment" %in% class(D))

  # trigger missing-arg error here if needed
  stat_name

  # find statistical result (full $output block)
  allstats <- D %>% mtm_res_get_entries(c("stats"))
  statind <- allstats %>%
    purrr::map("output") %>% purrr::map("name") %>%
    purrr::map(~ . == stat_name) %>% unlist() %>% which()
  if (length(statind) == 0) stop(sprintf("comparison '%s' not found", stat_name))
  output <- allstats %>% purrr::map("output") %>% .[[statind[1]]]

  # initialize header
  logtxt <- sprintf("'%s' info<br/><br/>", stat_name)

  # ----- Outcome & usage handling -----
  outcome_vec <- D %>% colData() %>% .[[output$outcome]]
  used <- as.logical(output$samples.used)
  used[is.na(used)] <- FALSE

  n_all  <- length(used)
  n_used <- sum(used, na.rm = TRUE)
  n_drop <- n_all - n_used

  # decide categorical vs numeric (binary numeric -> categorical)
  unique_used_vals <- unique(stats::na.omit(outcome_vec[used]))
  is_categorical <- !is.numeric(outcome_vec) ||
    (is.numeric(outcome_vec) && length(unique_used_vals) <= 2)

  log_lines <- character()

  # Line 1: N (outcome)
  if (is_categorical) {
    Ntab <- table(outcome_vec[used], useNA = "no")
    if (length(Ntab) > 0) {
      outcome_str <- paste(names(Ntab), as.integer(Ntab), sep = ": ", collapse = ", ")
      log_lines <- c(log_lines, sprintf("N (outcome): %s", outcome_str))
    } else {
      log_lines <- c(log_lines, "N (outcome): (none)")
    }
  } else {
    log_lines <- c(log_lines, "N (outcome): [numeric outcome; per-group N not applicable]")
  }

  # Line 2: N (total)
  log_lines <- c(log_lines, sprintf("N (total): %d of %d (dropped: %d)", n_used, n_all, n_drop))

  # ----- Extras -----
  if (is_categorical) {
    # per-outcome used/unused counts
    used_fac <- factor(used, levels = c(FALSE, TRUE))
    tab_full <- table(outcome_vec, used_fac, useNA = "no")
    dn <- dimnames(tab_full)
    if (is.null(dn[[2]]) || !("TRUE" %in% dn[[2]])) {
      outcome_levels <- dn[[1]]
      tab_full <- cbind("FALSE" = tab_full[, "FALSE", drop = TRUE],
                        "TRUE"  = 0L)
      rownames(tab_full) <- outcome_levels
    }
    used_unused <- apply(tab_full, 1L, function(v) sprintf("%d/%d", v["TRUE"], v["FALSE"]))
    if (is.null(names(used_unused))) names(used_unused) <- rownames(tab_full)
    uu_str <- paste(sprintf("%s: %s", names(used_unused), used_unused), collapse = ", ")
    log_lines <- c(log_lines, sprintf("N (used/unused by outcome): %s", uu_str))

    # simple stats of per-group Ns (used only)
    N_vec <- as.integer(tab_full[, "TRUE", drop = TRUE])
    if (length(N_vec) >= 1) {
      mean_n <- mean(N_vec)
      sd_n   <- if (length(N_vec) >= 2) stats::sd(N_vec) else NA_real_
      min_n  <- min(N_vec)
      max_n  <- max(N_vec)
      stats_str <- sprintf("mean=%.2f, sd=%s, min=%d, max=%d",
                           mean_n,
                           ifelse(is.na(sd_n), "NA", sprintf("%.2f", sd_n)),
                           min_n, max_n)
      log_lines <- c(log_lines, sprintf("N (group stats): %s", stats_str))
    }
  } else {
    # numeric outcome summary among USED samples
    y_used <- outcome_vec[used]
    if (n_used > 0) {
      m  <- mean(y_used, na.rm = TRUE)
      s  <- stats::sd(y_used, na.rm = TRUE)
      mn <- suppressWarnings(min(y_used, na.rm = TRUE))
      mx <- suppressWarnings(max(y_used, na.rm = TRUE))
      extras <- sprintf("Outcome summary: mean=%.3g, sd=%s, min=%.3g, max=%.3g",
                        m,
                        ifelse(is.finite(s), sprintf("%.3g", s), "NA"),
                        mn, mx)
      log_lines <- c(log_lines, extras)
    }
  }

  # Optional: result table filter hits
  if (!missing(stat_filter)) {
    s <- dplyr::enquo(stat_filter) %>% as.character()
    s <- gsub("~", "", s)
    hits <- output$table %>% dplyr::filter(!!dplyr::enquo(stat_filter)) %>% nrow()
    log_lines <- c(log_lines, sprintf("features: %s, %d of %d", s, hits, nrow(output$table)))
  }

  # finalize log text
  logtxt <- paste0(logtxt, paste(log_lines, collapse = "<br/>"))

  # record
  funargs <- mti_funargs()
  D %<>% mti_generate_result(
    funargs = funargs,
    logtxt  = logtxt
  )

  # return unchanged object
  D
}
