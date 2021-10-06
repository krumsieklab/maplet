#' Generate feature ratios
#'
#' @description
#' Transforms the dataset into a new dataset where each 'feature' represents a ratio of two features. Warning: For a dataset
#' with originally p features, this will result in p*(p-1) new variables (e.g. 500 features becomes 249500 ratios).
#'
#' @description
#' mt_post_pgain provides a special operation on a ratio data matrix for better interpretation of the resulting p-values.
#'
#' @param D \code{SummarizedExperiment} object.
#' @param stat_name Name of previous network generation call (e.g. \link{mt_stats_cormat_genenet}). Default: None (i.e. no
#'    network-based ratios).
#' @param edge_filter Filter criterion for edge selection, e.g. "p.adj < 0.05", as a term.
#' @param neighborhood Neighborhood degree to use (e.g. first neighbors, second neighbors). Default: 1.

#'
#' @examples
#' \dontrun{# Transform dataset to ratios
#' ... %>%  mt_modify_ratios() %>% ... # proceed with statistical analysis
#' }
#'
#' @return assay: Pairwise ratios from all variables of input.
#' @return rowData: Pairwise ratios from all variables of input.
#'
#' @author JZ, JK
#'
#' @export
mt_modify_ratios <- function(D, stat_name, edge_filter, neighborhood = 1){

    stopifnot("SummarizedExperiment" %in% class(D))

    as <- assay(D)
    p <- nrow(D)
    ## as <- matrix(1:(4*2), nrow = 4, ncol = 2, dimnames = list(letters[23:26], letters[1:2]))

    ## FOLDCHANGE FUNCTION (CONSIDER PREVIOUS LOG)
    op <- "/"
    if (length(maplet::mtm_res_get_entries(D, c("pre","trans","log"))) > 0){
        mti_logstatus("data already logscale, using '-'")
        op <- "-"
    }

    ## CREATE RATIOS
    as_ratio <- purrr::map(1:(nrow(as)-1), ~ sweep(as[(.x+1):nrow(as), , drop = F], 2, as[.x,], op)) %>%
        stats::setNames(rownames(as)[1:(nrow(as)-1)])

    ## CREATE NEW ROWDATA
    rd <- D %>%
        rowData() %>%
        as.data.frame() %>%
        dplyr::mutate(rownames = rownames(D))
    rd_new <- as_ratio %>%
        purrr::imap_dfr(~tibble::tibble(m1 = rownames(.x), m2 = .y)) %>%
        dplyr::mutate(rownames = stringr::str_c(m1, m2, sep = "_")) %>%
        dplyr::left_join(rd %>% dplyr::transmute(m1 = rownames, name1 = name), by = "m1") %>%
        dplyr::left_join(rd %>% dplyr::transmute(m2 = rownames, name2 = name), by = "m2") %>%
        dplyr::mutate(name = stringr::str_c(name1, " / ", name2))
    rd <- rd %>%
        dplyr::mutate(m1 = rownames,
                      m2 = NA,
                      name1 = name,
                      name2 = NA)
    rd <- dplyr::bind_rows(rd, rd_new) %>%
        dplyr::select(rownames, m1, m2, name1, name2, dplyr::everything())
    rownames(rd) <- NULL
    rd %<>% tibble::column_to_rownames("rownames")

    ## COMBINE RATIOS TO SINGLE MATRIX
    as_ratio <- as_ratio %>%
        purrr::imap(~{rownames(.x) <- stringr::str_c(rownames(.x), .y, sep = "_"); .x}) %>%
        purrr::invoke(rbind, .)
    as_ratio <- rbind(as, as_ratio)

    ## CHECK NAMES
    if(!identical(rownames(as_ratio), rownames(rd)))
        stop("something went wrong. check data!")

    ## FILTER BY NETWORK-BASED RATIOS (NBRs)?
    if (!missing(stat_name)) {
        # verify that arguments are given
        if (missing(edge_filter)) stop("If stat_name is given, edge_filter must be supplied.")
        # retrieve statistics object
        res <- D %>% maplet::mtm_get_stat_by_name(stat_name)
        # extract feature pairs according to formula, only keep feature names
        mpairs <- res %>% dplyr::filter(!!dplyr::enquo(edge_filter)) %>% dplyr::select(var1,var2)
        # initialize adjacency matrix
        A <- matrix(0, nrow=nrow(D), ncol=nrow(D))
        colnames(A) <- rownames(A) <- rownames(D)
        # build adjacency matrix using indices of mpairs in matrix
        inds <- data.frame(match(mpairs[,1], colnames(A)),match(mpairs[,2], colnames(A))) %>% as.matrix()
        A[inds] <- 1
        A[inds[,c(2,1)]] <- 1 # symmetric
        diag(A) <- 1
        # get k-neighborhood (A^k matrix multiplication),
        # found this trick on StackExchange, repeated application of %*%, also works with 1
        N <- Reduce("%*%", replicate(neighborhood, A, FALSE)) > 0
        # convert back to pairs
        mneighborpairs <- apply(which(N, arr.ind = T), 2, function(i){colnames(A)[i]}) %>%
            #  build string names as M1_M2
            apply(1, function(row){paste0(row[1],"_",row[2])})
        # match any pair that's in this list AND any single feature
        rn <- as_ratio %>% rownames()
        keep <- (rn %in% mneighborpairs) | !grepl("_", rn) # no "_" in it -> single feature
        # subselect
        as_ratio <- as_ratio[keep,]
        rd <- rd[keep,]
    }

    ## CREATE NEW OBJECT
    D <- SummarizedExperiment(assay    = as_ratio,
                              rowData  = rd,
                              colData  = colData(D),
                              metadata = metadata(D))

    ## add status information & plot
    funargs <- mti_funargs()
    D %<>% 
        mti_generate_result(
            funargs = funargs,
            logtxt = sprintf("created %d feature ratios out of %d features", nrow(D), p),
            output = NULL
        )
    D

}


# #### k-neighborhood test code ----
# # chain of 4, undirected
# A <- matrix(0,nrow=4,ncol=4)
# diag(A) <- 1
# A[1,2] <- A[2,3] <- A[3,4] <- 1
# A[2,1] <- A[3,2] <- A[4,3] <- 1
#
# # k neighborhood (A^k matrix multiplications)
# k=3
# Reduce("%*%", replicate(k, A, FALSE))




