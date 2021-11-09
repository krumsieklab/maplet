#' Creates KEGG pathways visualization
#'
#' Wrapper function for the \code{pathview::pathview} function. Maps and renders user data on relevant pathway graphs.
#'
#' @param D \code{SummarizedExperiment} input.
#' @param gene_id Name of the rowData column containing the gene identifiers.
#' @param gene_data Either vector (single sample) or a matrix-like data (multiple sample). Vector should be numeric with gene IDs
#'    as names or it may also be character of gene IDs. Character vector is treated as discrete or count data. Matrix-like data
#'    structure has genes as rows and samples as columns. Row names should be gene IDs. Here gene ID is a generic concepts,
#'    including multiple types of gene, transcript and protein uniquely mappable to KEGG gene IDs. KEGG ortholog IDs are also
#'    treated as gene IDs as to handle metagenomic data. Check details for mappable ID types. Default gene_data=NULL.
#' @param gene_idtype The ID type used for the gene_data, case insensitive. Default gene_idtype="entrez", i.e. Entrez Gene,
#'    which are the primary KEGG gene ID for many common model organisms. For other species, gene_idtype should be set to "KEGG"
#'    as KEGG use other types of gene IDs. For the common model organisms (to check the list, do: data(bods); bods), you may also
#'    specify other types of valid IDs. To check the ID list, do: data(gene_idtype.list); gene_idtype.list.
#' @param metab_id Name of the rowData column containing the feature identifiers.
#' @param cpd_data Same as gene_data, excpet named with IDs mappable to KEGG compound IDs. Over 20 types of IDs included in
#'    CHEMBL database can be used here. Check details for mappable ID types. Default cpd_data=NULL. Note that gene_data and cpd_data
#'    can't be NULL simultaneously.
#' @param cpd_idtype The ID type used for the cpd_data. Currently only works with "kegg".
#' @param stat_name Name of the statistics object to apply metab_filter to.
#' @param metab_filter If given, filter will be applied to data and only variables satisfying the condition will be included.
#' @param color_scale If given, this will be used to map colors to a continuous scale.
#' @param color_range A positive numeric value. If given, indicates the color range (-color_range, +color_range). If missing,
#'    color_range will be determined internally.
#' @param show_only_filtered If TRUE generate pathway list only based on filtered variables, otherwise pathways will be generated
#'    based on all variables.
#' @param low,mid,high Each is a list of two colors with "gene" and "cpd" as the names. This argument specifies the color spectra
#'    to code gene_data and cpd_data. Default spectra (low-mid-high) "green"-"gray"-"red" and "yellow"-"gray"-"blue" are used for
#'    gene_data and cpd_data respectively. The values for 'low, mid, high' can be given as color names ('red'), plot color index
#'    (2=red), and HTML-style RGB, ("\#FF0000"=red).
#' @param pathway_id A character vector. The KEGG pathway ID(s), usually 5 digit, may also include the 3 letter KEGG species code.
#'    If missing, the function will find all KEGG pathway annotations for the given KEGG identifiers.
#' @param n_pathways OPTIONAL. Number of pathways to output. Most populated pathway will be plotted first.
#' @param path_database The directory path of KEGG pathway data file (.xml) and image file (.png). If the path does not exist,
#'    the function will create it. Default: "./Pathview_database" (subfolder in the current working directory).
#' @param path_output The directory path of the function output files. If the path does not exist, the function will create it.
#'    Default: "./Pathview_output" (subfolder in the current working directory).
#' @param same_layer Controls if node colors are to be plotted in the same layer as the pathway graph. If FALSE, output
#'    generation will be faster, but output plots will be larger in size.
#' @param out_suffix The suffix to be added after the pathway name as part of the output graph file. Default: "pathview".
#' @param add_pwname_suffix If TRUE will add the pathway name to the output filename. If FALSE, will use what stored
#'    in out_suffix for all files. Default: FALSE.
#' @param db_path path to the pathway database to read annotations from.
#' @param limit A list of two numeric elements with "gene" and "cpd" as the names. This argument specifies the limit values for
#'    gene.data and cpd.data when converting them to pseudo colors.
#' @param \dots  See \code{pathview::pathview} for additional pathview arguments.
#'
#' @return pathview images (external)
#'
#' @examples
#' \dontrun{# plot all pathways with at least one significant feature from the statistical comparison "comp" in them
#' mt_plots_pathview(D = D,
#'                   metab_id="KEGG_mapped",
#'                   stat_name = "comp",
#'                   color_scale = -sign(fc)*log10(p.adj),
#'                   color_range = -log10(0.01),
#'                   metab_filter = p.adj < 0.05,
#'                   show_only_filtered = TRUE,
#'                   path_database = "./Pathview_database",
#'                   path_output = "./results/pathview",
#'                   same_layer = F,
#'                   add_pwname_suffix = T
#'                   ) %>%
#'                   ...}
#'
#' @author EB
#'
#' @import ggplot2
#' @import pathview
#'
#' @export
mt_plots_pathview <- function(D,
                              gene_id = NULL,
                              gene_data = NULL,
                              gene_idtype = "entrez",
                              metab_id = NULL,
                              cpd_data = NULL,
                              cpd_idtype = "kegg",
                              stat_name,
                              metab_filter,
                              color_scale,
                              color_range,
                              show_only_filtered = FALSE,
                              low = list(gene = "green", cpd = "yellow"),
                              mid =list(gene = "gray", cpd = "gray"),
                              high = list(gene = "red", cpd ="blue"),
                              pathway_id,
                              n_pathways,
                              path_database = "./Pathview_database",
                              path_output = "./Pathview_output",
                              same_layer = TRUE,
                              out_suffix = "pathview",
                              add_pwname_suffix = FALSE,
                              db_path = system.file("extdata", "precalc/pathview/KeggPathways.Rds", package = "maplet"),
                              limit = list(gene = 1, cpd = 1),
                              ...) {

  requireNamespace(pathview)
  data("bods", package = "pathview")

  # get argument names from dots
  n <- sapply(as.list(substitute(list(...)))[-1L], deparse)
  dot_args <- names(n)

  # check for defunct argument names
  if ("gene.id" %in% dot_args) stop("You used the old MT naming convention gene.id. Should be: gene_id.")
  if ("gene.data" %in% dot_args) stop("You used the old MT naming convention gene.data. Should be: gene_data.")
  if ("gene.idtype" %in% dot_args) stop("You used the old MT naming convention gene.idtype. Should be: gene_idtype.")
  if ("met.id" %in% dot_args) stop("You used the old MT naming convention met.id. Should be: metab_id.")
  if ("cpd.data" %in% dot_args) stop("You used the old MT naming convention cpd.data. Should be: cpd_data.")
  if ("cpd.idtype" %in% dot_args) stop("You used the old MT naming convention cpd.idtype. Should be: cpd_idtype.")
  if ("metab.filter" %in% dot_args) stop("You used the old MT naming convention metab.filter Should be: metab_filter.")
  if ("color.scale" %in% dot_args) stop("You used the old MT naming convention color.scale Should be: color_scale.")
  if ("color.range" %in% dot_args) stop("You used the old MT naming convention color.range. Should be: color_range.")
  if ("show.only.filtered" %in% dot_args) stop("You used the old MT naming convention show.only.filtered. Should be: show_only_filtered.")
  if ("pathway.id" %in% dot_args) stop("You used the old MT naming convention pathway.id Should be: pathway_id")
  if ("n.pathways" %in% dot_args) stop("You used the old MT naming convention n.pathways. Should be: n_pathways.")
  if ("path.database" %in% dot_args) stop("You used the old MT naming convention path.database. Should be: path_database.")
  if ("path.output" %in% dot_args) stop("You used the old MT naming convention path.output. Should be: path_output.")
  if ("same.layer" %in% dot_args) stop("You used the old MT naming convention same.layer. Should be: same_layer.")
  if ("out.suffix" %in% dot_args) stop("You used the old MT naming convention out.suffix. Should be: out_suffix.")
  if ("add.pwname.suffix" %in% dot_args) stop("You used the old MT naming convention add.pwname.suffix. Should be: add_pwname_suffix.")
  
  ## check input
  stopifnot("SummarizedExperiment" %in% class(D))
  ## if both gene_id and gene_data are given, throw an error
  if(!is.null(gene_id) & !is.null(gene_data))
    stop("The function can only use either gene_data or gene_id. Please provide only one of the two.")
  ## if both metab_id and cpd_data are given, throw an error
  if(!is.null(metab_id) & !is.null(cpd_data))
    stop("The function can only use either cpd_data or metab_id. Please provide only one of the two.")
  ## if gene_id is provided, check that is is a valid column name of the rowData
  if(!is.null(gene_id)) {
    if (!(gene_id %in% colnames(rowData(D))))
      stop(sprintf("%s is not contained in the rowData of the Summarized Experiment", gene_id))}
  ## if metab_id is provided, check that is is a valid column name of the rowData
  if(!is.null(metab_id) ) {
    if(!(metab_id %in% colnames(rowData(D))))
      stop(sprintf("%s is not contained in the rowData of the Summarized Experiment", metab_id))}
  ## if gene_id is provided, check that it is length 1
  if(!is.null(gene_id)) {
    if(length(gene_id)!=1)
      stop(sprintf("%s can only be a single column", gene_id))}
  ## if metab_id is provided, check that it is length 1
  if(!is.null(metab_id)) {
    if(length(metab_id)!=1)
      stop(sprintf("%s can only be a single column", metab_id))}
  ## if metab_filter is given, stat_name must also be given and either metab_id or gene_id must be given as well
  if(!missing(metab_filter)) {
    if(missing(stat_name))
      stop("In order to use metab_filter, stat_name must be given")
    if(is.null(gene_id) & is.null(metab_id))
      stop("In order to use metab_filter, one betweeen gene_id and metab_id must be given")}
  ## if n.pathway is given, it must be numeric
  if(!missing(n_pathways)) {
    if(class(n_pathways)!="numeric")
      stop("n_pathways must be numeric")}
  ## if show_only_filtered is TRUE, metab_filter must be given
  if(show_only_filtered) {
    if(missing(metab_filter))
      stop("show_only_filtered can be TRUE only if metab_filter is given")}
  ## in order for add_pwname_suffix to work when TRUE, pathway_id must be missing
  if(add_pwname_suffix){
    if(!(missing(pathway_id)))
      stop("add_pwname_suffix can only be TRUE if pathway_id is missing")
  }

  ## rowData
  rd <- rowData(D) %>%
    as.data.frame() %>%
    dplyr::mutate(var = rownames(D))

  ## stat
  if(!missing(stat_name)){
    stat <- mtm_get_stat_by_name(D, stat_name) %>%
      dplyr::inner_join(rd, by = "var")
  }else{
    stat <- rd
  }

  ## Add color variable
  if(!missing(color_scale)){
    color_scale_q <- dplyr::enquo(color_scale)
    # add color variable according to input
    stat <- stat %>%
      dplyr::mutate(color=!!color_scale_q)
  } else {
    # if not given, set color to 1
    stat <- stat %>%
      dplyr::mutate(color=1)
  }

  ## FILTER METABOLITES
  if(!missing(metab_filter)){
    metab_filter_q <- dplyr::enquo(metab_filter)
    # filter results
    var <- stat %>%
      dplyr::filter(!!metab_filter_q)
    # collect variable names of filtered results
    var <- var$var
    # set color variable of filtered out variables to 0
    stat$color[!(stat$var %in% var)] <- 0
  }

  # if gene_id is provided, extract identifiers from the rowData
  if(!is.null(gene_id)) {
    # remove duplicated identifiers if they occur
    stat <- stat[!duplicated(stat[[gene_id]]),]
    # create pathview variable
    gene_data <- data.frame(color=stat$color[!is.na(stat[[gene_id]])])
    # remove rows with NAs in the identifiers
    rownames(gene_data) <- stat[[gene_id]][!is.na(stat[[gene_id]])]
  }
  # if metab_id is provided, extract identifiers from the rowData
  if(!is.null(metab_id)) {
    # remove duplicated identifiers if they occur
    stat <- stat[!duplicated(stat[[metab_id]]),]
    # create pathview variable
    cpd_data <- data.frame(color=stat$color[!is.na(stat[[metab_id]])])
    # remove rows with NAs in the identifiers
    rownames(cpd_data) <- stat[[metab_id]][!is.na(stat[[metab_id]])]
  }

  # if path_output is provided, check if the folder exists, otherwise create it
  if (!file.exists(path_output)){
    dir.create(path_output)
  }
  # if path_database is provided, check if the folder exists, otherwise create it
  if (!file.exists(path_database)){
    dir.create(path_database)
  }

  # if no pathway_id list is provided, find annotations for kegg identifiers
  if(missing(pathway_id)) {

    # load KEGG pathway database
    load(db_path)

    # build one big dataframe with all pathway informations
    pwdf <- do.call(rbind, pwdb)

    if (!is.null(gene_data)) {
      if(!is.null(rownames(gene_data))) {
        ids <- rownames(gene_data)
        if(show_only_filtered) {
          ids <- ids[ids %in% rownames(gene_data)[gene_data$color != 0]]
        }
      } else {
        ids <- gene_data
      }

      if(length(ids) != 0) {
        # find gene pathway annotations
        g_anno <- lapply(ids, function(x) {
          pwdf$ID[pwdf$src==x] %>% unique()
        })
        names(g_anno) <- ids
        # build one long list
        g_anno_list <- do.call(c, g_anno)
        # find most common pathway for genes
        pw_gene <- g_anno_list %>% table() %>% as.data.frame()
        colnames(pw_gene) <- c("pathway","Freq")
        # pathway list ordered according to the number of genes with that annotation
        pw_gene <- pw_gene[order(pw_gene$Freq,decreasing = TRUE),]
        # find names of these pathways
        g_pw_names <- lapply(pw_gene$pathway, function(x) {
          pwdf$name[pwdf$ID==x] %>% unique()
        })
        names(g_pw_names) <- pw_gene$pathway
        # build one long list
        pw_names <- do.call(c, g_pw_names)
        # remove ":" from pathway ids for pathview
        pw_gene$pathway <- gsub(":", "", pw_gene$pathway)
        names(pw_names) <- gsub(":", "", names(pw_names))
        pw <- pw_gene
      } else {
        mti_logwarning("Filtering returned an empty matrix")
        pw <- list()
        pw$pathway <- NULL
        gene_data <- NULL
      }
      if(!(is.null(pw$pathway))) {
        # save pathway list only if list of variables to output is not empty
        pathway_id <- pw$pathway
      }
    }

    if (!is.null(cpd_data)) {
      if(!is.null(rownames(cpd_data))) {
        ids <- rownames(cpd_data)
        if(show_only_filtered) {
          ids <- ids[ids %in% rownames(cpd_data)[cpd_data$color != 0]]
        }
      } else {
        ids <- cpd_data
      }

      if(length(ids)!=0) {
        # find feature pathway annotations
        m_anno <- lapply(ids, function(x) {
          pwdf$ID[pwdf$dest==x] %>% unique()
          # cbind(pwdf$ID[pwdf$dest==x] %>% unique(),pwdf$name[pwdf$dest==x] %>% unique())
        })
        names(m_anno) <- ids
        # build one long list
        m_anno_list <- do.call(c, m_anno)
        # find most common pathway for features
        pw_met <- m_anno_list %>% table() %>% as.data.frame()
        colnames(pw_met) <- c("pathway","Freq")
        # pathway list ordered according to the number of features with that annotation
        pw_met <- pw_met[order(pw_met$Freq,decreasing = TRUE),]
        # find names of these pathways
        m_pw_names <- lapply(pw_met$pathway, function(x) {
          pwdf$name[pwdf$ID==x] %>% unique()
        })
        names(m_pw_names) <- pw_met$pathway
        # build one long list
        pw_names <- do.call(c, m_pw_names)
        # remove ":" from pathway ids for pathview
        pw_met$pathway <- gsub(":", "", pw_met$pathway)
        names(pw_names) <- gsub(":", "", names(pw_names))
        pw <- pw_met
      } else {
        mti_logwarning("Filtering returned an empty matrix")
        pw <- list()
        pw$pathway <- NULL
        cpd_data <- NULL
      }
      if(!(is.null(pw$pathway))) {
        # save pathway list only if list of variables to output is not empty
        pathway_id <- pw$pathway
      }
    }

    if (!is.null(gene_data) & !is.null(cpd_data)) {
      # find most common pathway for both genes and features
      pw_list <- c(m_anno_list,g_anno_list)
      pw <- pw_list %>% table() %>% as.data.frame()
      colnames(pw) <- c("pathway","Freq")
      # pathway list ordered according to the number of features/genes with that annotation
      pw <- pw[order(pw$Freq,decreasing = TRUE),]
      # find names of these pathways
      pw_names <- lapply(pw$pathway, function(x) {
        pwdf$name[pwdf$ID==x] %>% unique()
      })
      names(pw_names) <- pw$pathway
      # build one long list
      pw_names <- do.call(c, pw_names)
      # remove ":" from pathway ids for pathview
      pw$pathway <- gsub(":", "", pw$pathway)
      names(pw_names) <- gsub(":", "", names(pw_names))

      if(!(is.null(pw$pathway))) {
        # save pathway list only if list of variables to output is not empty
        pathway_id <- pw$pathway
      }
    }
  }

  ## Set color scale limits
  if(!missing(color_range)) {
    limit = list(gene=color_range, cpd=color_range)
  } else {
    # if(!is.null(metab_id)) {
    #   limit$cpd = max(ceiling(abs(stat$color)), na.rm = TRUE)
    # }
    # if(!is.null(gene_id)) {
    #   limit$gene = max(ceiling(abs(stat$color)), na.rm = TRUE)
    # }
    # limit=list(gene=max(ceiling(abs(stat$color))), cpd=max(ceiling(abs(stat$color))))
    if(!is.null(cpd_data)) {
      if(!is.null(rownames(cpd_data))) {
        limit$cpd = max(ceiling(abs(cpd_data)), na.rm = TRUE)
      }
    }
    if(!is.null(gene_data)) {
      if(!is.null(rownames(gene_data))) {
        limit$gene = max(ceiling(abs(gene_data)), na.rm = TRUE)
      }
    }
  }
  if(limit$gene == 0) {limit$gene = 1}
  if(limit$cpd == 0) {limit$cpd = 1}

  # move working directory to kegg.dir (otherwise some files will be saved in the working directory even if another directory is provided)
  wd <- getwd()
  setwd(path_database)
  save.path <- getwd()
  setwd(wd)
  setwd(path_output)

  if(missing(pathway_id)) {
    file.create(paste0(getwd(),"/NO_RESULTS_AFTER_FILTERING.txt",sep=""))
    pv.out <- NULL
  } else {

    # removing problematic pathways from list if present -> they throw a weird error
    pathway_id <- pathway_id[!(pathway_id %in% c("mmu05206","mmu04666","mmu04723"))]

    if(!missing(n_pathways)) {
      if(n_pathways>length(pathway_id))
        mti_logwarning(sprintf("n.pathway is %i, but there are only %i pathways, so %i pathways will be used", n_pathways, length(pathway_id), length(pathway_id)))
      pathway_id <- pathway_id[1:min(n_pathways,length(pathway_id))]
      pw_names <- pw_names[1:min(n_pathways,length(pathway_id))]
    }

    # print(sprintf("%d pathways detected", length(pathway_id)))

    mti_logstatus(glue::glue("There are {length(pathway_id)} pathways detected. This function may take a few minutes to run."))
    lapply(1:length(pathway_id), function(x) {
      suppressMessages(
        pv.out <- pathview::pathview(gene.data = gene_data, cpd.data = cpd_data, pathway.id = pathway_id[x], kegg.dir = save.path,
                                     cpd.idtype = cpd_idtype, gene.idtype = gene_idtype, limit = limit, low = low,
                                     mid = mid, high = high, same.layer = same_layer, out.suffix = out_suffix, ...)
      )
      # if pathview file not available, pathview returns 0; if pv.out is a double, the file download failed
      if(typeof(pv.out)=="double"){
        mti_logwarning(glue::glue("Download of files for {pathway_id[x]} failed! This pathway may not exist!"))
      }else if(add_pwname_suffix) {
        # add pathway rank and name to filename
        fname <- list.files(".",pattern=pathway_id[x])
        # isolate pathway name from filename
        m <- pw_names[names(pw_names)==substr(fname[length(fname)], 1, 8)]
        m <- gsub('[[:punct:]]+','',m)
        m <- stringr::str_replace_all(m," ","_")
        file.rename(from=fname,to=sub(pattern=sprintf("%s.%s",pathway_id[x], out_suffix),replacement=sprintf("%d_%s.%s",x,m,pathway_id[x]),fname))
      }
    }) %>% invisible()
  }

  setwd(wd)

  n_pw <- ifelse((!missing(pathway_id)),length(pathway_id),0)
  nn <- ifelse((!missing(pathway_id)),pw_names,NA)

  # add status information & plot
  funargs <- mti_funargs()
  D %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("%d Pathway images saved in %s", n_pw, path_output),
      output = NULL,
      output2 = nn
    )

  # return
  D
}
