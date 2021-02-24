#' Generate a new data matrix that is one group minus a second group (e.g. difference between two timepoints)
#'
#' Takes the name of a sample factor variable (from colData), and two factor levels. Creates a new data matrix that
#' is made of the second minus the first levels. Will only produce results for cases where both data points were present in the
#' dataset.
#'
#' @param D \code{SummarizedExperiment} input.
#' @param id_col Name of colData column representing sample ID variable  (e.g. patient ID).
#' @param group_col Name of colData column representing sample group variable (e.g. timepoint).
#' @param grp1 Group level 1 (e.g. "day1").
#' @param grp2 Group level 2 (e.g. "day2").
#'
#' @return assay: Replacement containing difference of the two groups.
#' @return colData: Replacement containing only those samples from grp1.
#'
#' @examples
#' \dontrun{# Generate new dataset as differences between day2 and day1
#' ... %>%  mt_modify_diff_groups(id_col = "patientID", group_col = "timepoint", grp1 = "day1", grp2 = "day2")  %>% ... # proceed with statistical analysis
#' }
#'
#' @author JK
#'
#' @export
mt_modify_diff_groups <- function(D, id_col, group_col, grp1, grp2) {

  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))

  # if data is logged, then subtract, otherwise divide
  if (length(mtm_res_get_entries(D, c("pre","trans","log"))) > 0){
    op <- `-`
    opstr <- "-"
  } else {
    op <- `/`
    opstr <- "/"
  }


  # shortcuts
  cd = colData(D)
  X = assay(D)

  # verify that all fields exist
  if (!(id_col %in% colnames(cd))) stop(sprintf("'%s' not found in sample annotations.", id_col))
  if (!(group_col %in% colnames(cd))) stop(sprintf("'%s' not found in sample annotations.", group_col))

  # find matching samples
  s1 = cd[[group_col]] == grp1
  s2 = cd[[group_col]] == grp2

  # extract first group
  ids1 = cd[[id_col]][s1]
  dat1 = X[,s1]
  # extract second group
  ids2 = cd[[id_col]][s2]
  dat2 = X[,s2]

  # map the two groups
  idc = intersect(ids1,ids2)
  m1 = match(idc,ids1)
  m2 = match(idc,ids2)

  # subtract or divide data frames
  diffmat = op(dat2[,m2], dat1[,m1])
  # copy the clinical data of the first group 1
  cd_sub <- cd[s1,]
  cd_new <- cd_sub[match(idc, cd_sub[[id_col]]),]

  # create new SummarizedExperiment object
  D <- SummarizedExperiment(assay    = diffmat,
                            rowData  = rowData(D),
                            colData  = cd_new,
                            metadata = metadata(D))

  # add status information
  funargs <- mti_funargs()
  metadata(D)$results %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = glue('created new dataset, id var: {id_col}, subtracting {group_col}, group {grp2}{opstr}{grp1}')
    )

  # return
  D


}
