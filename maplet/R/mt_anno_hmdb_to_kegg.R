#' Create KEGG identifiers from HMDB identifiers
#'
#' Creates annotations and merges them into current \code{SummarizedExperiment} object. Performs "left-joins", i.e. leaves the
#' original SE unchanged and just adds information where it can be mapped.
#'
#' @param D \code{SummarizedExperiment} input.
#' @param in_col Name of the rowData column containing the HMDB identifiers.
#' @param out_col Name of the new rowData column to be created with the KEGG identifiers. Default="KEGG".
#'
#' @return rowData: New annotation column added.
#'
#' @examples
#' # Load data, sheet with sample annotations, sheet with metabolite annotations and add KEGG identifiers
#' \dontrun{D <-
#'   # load raw data
#'   mt_load_xls(file=file, sheet="data", samples_in_rows=T, id_col="SAMPLE_NAME") %>%
#'   # sample annotations from metabolomics run
#'   mt_anno_xls(file=file, sheet="sampleinfo", anno_type="samples", anno_id_col = "SAMPLE_NAME") %>%
#'   # metabolite annotations
#'   mt_anno_xls(file=file, sheet="metinfo", anno_type="metabolites", anno_id_col="BIOCHEMICAL", data_id_col = "name") %>%
#'   # add KEGG identifiers
#'   mt_anno_hmdb_to_kegg(in_col="HMDb_id", out_col="KEGG") %>%
#'   ...}
#'
#' @author EB
#'
#' @export
mt_anno_hmdb_to_kegg <- function(D, in_col, out_col = "KEGG") {

  ## check input
  stopifnot("SummarizedExperiment" %in% class(D))
  if(missing(in_col))
    stop("in_col must be given")
  # check that in_col is a string
  if(!(class(in_col)=="character"))
    stop("in_col must be a string")
  # check that in_col is a valid column names of the rowData
  if(!(in_col %in% colnames(rowData(D))))
    stop(sprintf("%s is not contained in the rowData of the Summarized Experiment", in_col))
  if(out_col %in% colnames(rowData(D)))
    stop(sprintf("%s already exists, please choose another name", out_col))

  # load look-up table of HMDB to KEGG identifiers
  load(system.file("extdata", "precalc/pathview/MetaboliteMapping.Rds", package = "maplet"))
  # get rowData
  df <- as.data.frame(rowData(D))
  # rename in_col for left_join use
  names(df)[which(names(df)==in_col)] <- "MappingIDs"

  # join with automatic metabolite mapping
  newdf <- df %>%
    dplyr::left_join(MetaboliteMapping[,c("secondary_accessions","kegg_id")], by = c("MappingIDs" = "secondary_accessions")) %>%
    dplyr::left_join(MetaboliteMapping[,c("accession","kegg_id")], by=c("MappingIDs" = "accession")) %>%
    dplyr::distinct()
  # merge KEGG identifiers
  newdf[[out_col]] <- newdf$kegg_id.x
  newdf[[out_col]] [is.na(newdf[[out_col]])] <- newdf$kegg_id.y[is.na(newdf[[out_col]])]
  # drop redundant columns
  newdf <- subset(newdf, select=-c(kegg_id.x,kegg_id.y))

  # go back to original name
  # rename in_col for left_join use
  names(newdf)[which(names(newdf)=="MappingIDs")] <- in_col

  # overwrite entries stored in the manually curated file
  dd <- utils::read.csv2(file=system.file("extdata", "precalc/pathview/MetaboliteMapping_manual.csv", package = "maplet"), sep=",")
  dd$Name <- as.character(dd$Name)
  dd$HMDB <- as.character(dd$HMDB)
  dd$KEGG <- as.character(dd$KEGG)

  # find identifiers to overwrite
  ow <- intersect(newdf[[in_col]],dd$HMDB)
  # subselect
  dd <- dd[dd$HMDB %in% ow,]
  # get order of identifiers in rowData
  dd <- dd[match(newdf[[in_col]][which(newdf[[in_col]] %in% ow)], dd$HMDB),]
  # overwrite identifiers
  newdf[[out_col]][which(newdf[[in_col]] %in% ow)] <- dd$KEGG
  # update rowData
  rowData(D) <- DataFrame(newdf)

  # add status information
  funargs <- mti_funargs()
  metadata(D)$results %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("loaded KEGG annotations for %i out of %i metabolites", length(newdf[[out_col]][!is.na(newdf[[out_col]])]), length(newdf[[out_col]]))
    )

  # return
  D

}
