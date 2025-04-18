% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mt_stats_univ_gam.R
\name{mt_stats_univ_gam}
\alias{mt_stats_univ_gam}
\title{Generalized additive models using mgcv package}
\usage{
mt_stats_univ_gam(
  D,
  outcome,
  outcome_type = "numeric",
  rval = NULL,
  conf_formula,
  int_w_analyte = NULL,
  stat_name,
  samp_filter
)
}
\arguments{
\item{D}{\code{SummarizedExperiment} input.}

\item{outcome}{Outcome to be used for response.}

\item{outcome_type}{Outcome type - numeric/binary/ordinal. Default: numeric}

\item{rval}{For ordinal outcomes. Default: NULL.}

\item{conf_formula}{Confounders to be corrected for.}

\item{int_w_analyte}{Name of covariate that interacts with metabolite. Default: NULL.}

\item{stat_name}{Name under which this comparison will be stored, must be unique to all other statistical results.}

\item{samp_filter}{Term which samples to filter to first (e.g. used if the data contains >2 groups but the user wants
to run a two-group comparison).}
}
\value{
$results$output: Statistics object.
}
\description{
Only supports standard linear models of the form: outcome ~ feature + [covariates].
}
\examples{
\donttest{# run lm with no confounders, "Group" as outcome
# filter to groups "Li_2" and "Li_5"
# name the comparison "Li's"
... \%>\%
 mt_stats_univ_gam(
   conf_formula      = ~ Group,
   samp_filter = (Group \%in\% c("Li_2","Li_5")),
   stat_name         = "Li's"
 ) \%>\% ...
 }

}
\author{
RB
}
