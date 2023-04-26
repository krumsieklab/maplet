#' Generate aggregated pathway values
#'
#' @description
#' Takes a pathway annotation column of the features (rowData) and builds one representative sample per pathway for each sample.
#' Also works for overlapping pathway annotations (i.e. where each feature can have >1 pathway).
#'
#' @description
#' Implemented approaches:
#' \enumerate{
#'   \item{Eigenmetabolite/eigengene/eigenvalue PCA based approach. Data matrix cannot have NAs.}
#'   \item{Mean value. Will NOT scale data before. Data matrix can have NAs.}
#' }
#'
#' @param D  \code{SummarizedExperiment} input.
#' @param pw_col Column name from rowData containing pathway annotations.
#' @param method One of: "eigen", "aggmean".
#'
#' @return assay: entirely replaces assay with new data matrix
#' @return rowData: removes original rowData, since now the "features" are pathways
#'
#' @examples
#' \dontrun{%>% mt_modify_agg_pathways(pw_col="SUB_PATHWAY", method="aggmean") %>%  # subpathways from metabolon
#'
#' # add KEGG pathways and use those
#' %>%
#'   mt_anno_pathways_hmdb(in_col = "HMDB", out_col = "kegg_db", pwdb_name = "KEGG", db_dir = system.file("extdata", "precalc/hmdb/", package = "maplet")) %>%
#'   mt_anno_pathways_remove_redundant(feat_id = "HMDB", pw_id = "kegg_db") %>%
#'   mt_modify_agg_pathways(pw_col="kegg_db", method="aggmean") %>%}
#'
#' @author JK
#'
#' @export
mt_modify_agg_pathways <- function(D, pw_col, method) {

  # remove all NAs from a vector
  # alternatively, replaces NAs with a value
  removeNAs <- function(v, replaceWith=NULL) {
    if (!is.null(replaceWith)) {
      v[is.na(v)] <- replaceWith
      v
    } else {
      v[!is.na(v)]
    }
  }


  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  if (!(method %in% c("eigen","aggmean"))) stop("'method' must be either 'eigen' or 'aggmean'")

  # get variable
  if (!(pw_col %in% colnames(rowData(D)))) stop(sprintf("'%s' not found in feature annotations.", pw_col))
  p = rowData(D)[[pw_col]]
  # ensure it's a factor or character vector
  if (!all(sapply(p,class) %in% c("NULL","character"))) stop(sprintf("'%s' has to be a list of character lists", pw_col))

  # collect all pathway names
  up <- unique(unlist(p)) %>% removeNAs()

  # agg calculation
  X = t(assay(D))
  output <- NULL
  if (method=="eigen") {
    # for eigenvalue methods, matrix cannot have any NAs
    if (any(is.na(X))) stop("no NA values allowed for 'eigen' method")
    # calc
    res <- up %>% lapply(function(g) {
      met <- removeNAs(sapply(p, function(v){g %in% v}), replaceWith=F)
      pca = stats::prcomp(as.data.frame(X[,met]))
      list(pc1=pca$x[,1],
           expvar=(pca$sdev)^2 / sum(pca$sdev^2) )
    })
    # assemble into matrix
    M = as.data.frame(sapply(res, function(x){x$pc1}))
    # return explained variances in output
    output$expvar = purrr::map(res, "expvar")

  } else if (method=="aggmean") {
    # aggregated mean
    M <- up %>% sapply(function(g){
      met <- removeNAs(sapply(p, function(v){g %in% v}), replaceWith=F)
      apply(as.data.frame(X[,met]),1,mean, na.rm=T  )
    })

  } else {stop("bug")}

  # prepare matrix
  full.names <- colnames(M)
  colnames(M) <- make.names(full.names)

  # check if data can be copied by grouping the pw_col variable
  res <- try(rowData(D) %>% as.data.frame() %>% tibble::as.tibble() %>% dplyr::group_by_(pw_col), silent = TRUE)
  copyworks <- !(class(res) == "try-error")

  # second check is there to avoid problems caused by nested pathways
  if (all(copyworks) & length(up) == length(unique(rowData(D)[[pw_col]]))) {
    # check which variables can be copied [all of this can probable be done simpler]
    copyover <- sapply(colnames(rowData(D)), function(c) {
      # verify variable, there must be only one value for each instance
      all( (rowData(D) %>% as.data.frame() %>% tibble::as.tibble() %>% dplyr::group_by_(pw_col) %>% dplyr::summarise(dplyr::n_distinct(!!rlang::sym(c))))[[2]] ==1 )
    })
    # generate new rowData
    rd <- rowData(D) %>% as.data.frame() %>% tibble::as.tibble() %>% dplyr::filter(!duplicated(!!rlang::sym(pw_col))) %>%
      dplyr::select(which(copyover)) %>% dplyr::filter(!is.na(!!rlang::sym(pw_col))) %>% dplyr::mutate(name=!!rlang::sym(pw_col)) %>% dplyr::select(name,dplyr::everything())
    # sanity bug check
    stopifnot(all.equal(rd$name, rd[[pw_col]]))

  } else {
    # just create rowdata with ID and name
    rd = data.frame(name= up,
                    pathway_name = unlist(purrr::map(up, ~  (metadata(D)$pathways[[pw_col]] %>% dplyr::filter(ID==.x))$pathway_name)))
    rownames(rd) <- up

  }

  # generate summarized experiment
  assay <- t(M)
  rownames(rd) <- assay %>% row.names()
  newD <- SummarizedExperiment(assay = assay,
                               colData  = colData(D),
                               rowData  = rd,
                               metadata = metadata(D)
  )

  # add status information
  funargs <- mti_funargs()
  newD %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("pathway aggregation by '%s', %d pathway scores generated", pw_col, nrow(newD))
    )

  # return
  newD


}







