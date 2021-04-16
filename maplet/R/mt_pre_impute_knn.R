#' Impute using KNN method
#'
#' @description
#' Default settings are from the winning approach of our paper.
#' https://link.springer.com/article/10.1007%2Fs11306-018-1420-2
#'
#' @description
#' Multi-core functionality is available by setting the flag use_multicore=T.
#
#' @param D \code{SummarizedExperiment} input.
#' @param method Specific KNN method to use. Default: "knn.obs.euc.sel" (recommendation: do not touch).
#' @param k Number of nearest neighbors to consider. Recommend not changing this parameter. Default: 10 (recommendation: do not touch).
#' @param verbose Output intermediate steps? Default: F.
#' @param use_multicore Use multicore KNN imputation method? Default: F.
#' @param n_cores If use_multicore==T, number of cores to use for imputation. Default: 5.
#'
#' @return assay: Imputed data.
#'
#' @examples
#' \dontrun{# in the context of a SE pipeline
#' ... %>% mt_pre_impute_knn() %>% ...    # standard call
#' }
#'
#' @author JK, PG
#'
#' @export
mt_pre_impute_knn <- function(D,
                              method="knn.obs.euc.sel",
                              k=10,
                              verbose=F,
                              use_multicore=F,
                              n_cores=5) {

  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  stopifnot(k%%1==0) # integer number

  if(use_multicore==T){
    # impute
    assay(D) =
      assay(D) %>%
      t() %>%
      imputeKNN_multicore(K=k, mc_cores = n_cores, verbose=verbose) %>%
      t()
  }else{
    # impute
    assay(D) = t(imputeKNN( t(assay(D)), methods=method, K=k, verbose=verbose) )
  }

  # add status information
  funargs <- mti_funargs()
  D %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = "Imputed via KNN"
    )

  # return
  D

}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# KNN imputation
# K. Do based on Simone Wahls script from Missing values paper

# to use winner of MV paper do: imputeKNN(dat, methods="knn.obs.euc.sel", K=10)

# 18.12.17    - KD,JK
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #



#' KNN imputation
#'
#' Code from Kiki Do's and Simone Wahl's paper
#' https://www.ncbi.nlm.nih.gov/pubmed/30830398
#'
#' To reproduce exact setting of winner method from paper:
#' imputeKNN(dat, methods="knn.obs.euc.sel", K=10)
#'
#' Note: Most parts of this code are currently not used, since we only use the above-mentioned parameter setting in MT.
#'
#' @param dat input data matrix
#' @param methods which variation of KNN imputation to perform
#' @param cor.var.sel correlation threshold for neighbor selection
#' @param K number of neighbors to use
#' @param verbose output status messages?
#'
#' @returns imputed data matrix
#'
#' @noRd
imputeKNN <- function(dat,
                      methods = c("knn.vars","knn.obs","knn.obs.euc","knn.obs.euc.sel"),
                      cor.var.sel = 0.2,
                      K=5,
                      verbose=T) {


  datimp <- dat
  ## knn.variable
  if("knn.vars" %in% methods) {
    incom.vars <- which(apply(dat,2,function(x) any(is.na(x))))
    if(verbose)message(paste0("Number of imcomplete variables: ", length(incom.vars)))

    # DMwR::knnImputation kann nicht damit umgehen, wenn weniger als k verf?gbare (obs) vars
    # also manuell:
    D2 <- as.matrix(stats::dist(t(scale(dat)),upper=T,diag=T)) # eucl
    D2[D2==0] <- NA
    if(verbose)message("Proceeding var... ")
    for (j in incom.vars) {
      if(verbose)cat(paste(j," "))
      comobs <- stats::complete.cases(dat[,j])
      Mean <- mean(dat[,j],na.rm=T)
      SD <- stats::sd(dat[,j],na.rm=T)
      if(any(!is.na(D2[,j]))) {
        KNNvars <- order(D2[,j],na.last=NA)
        KNNvars <- KNNvars[sapply(KNNvars, function(jj) any(!is.na(dat[!comobs,jj])))]
      } else KNNvars <- NULL

      dattmp <- dat
      KNNvars_sel <- KNNvars[1:min(K,length(KNNvars))]
      if(any(!is.na(D2[,j])) & length(KNNvars)>=1) dattmp[!comobs,j]  <-
        sapply(1:length(which(!comobs)),function(co) {
          dattmp_all <- scale(dattmp)[co,KNNvars]
          dattmp_sel <- scale(dattmp)[co,KNNvars_sel]
          if(any(!is.na(dattmp_sel))) ret <- sum((dattmp_sel*SD+Mean)*exp(-D2[KNNvars_sel,j]),na.rm=T)/sum(exp(-D2[KNNvars_sel,j])[!is.na(dattmp_sel)],na.rm=T) else ret <- na.omit(dattmp_all)[1]*SD+Mean
          ret})
      if(any(is.na(dattmp[,j]))) {
        still.incom <-  !stats::complete.cases(dattmp[,j])
        dattmp[still.incom,j] <- mean(dattmp[!still.incom,j],na.rm=T) }
      datimp[,j] <- dattmp[,j]
    }

  } else if("knn.obs" %in% methods) { ## knn.sample: KNN per sample using mahalanobis distance (since dimensions=variables do not have same scale)

    # yaImpute::yai verwendet nur ganz complete obs! --> manuell!
    # dasselbe gilt f?r mahalanobis und mahalanobis.dist[StatMatch], wobei letzteres wenigstens ganze Matrix macht
    #--> manuell
    Sx <- stats::cov(dat, use = "p")
    if(any(is.na(Sx))) {
      submatrix <- dat[-c(which(colSums(is.na(Sx))>0))]
      Sx <- stats::cov(submatrix, use = "p")}
    incom.obs <- which(apply(dat,1,function(x) any(is.na(x))))
    if(verbose)message(paste0("Number of imcomplete observations: ", length(incom.obs)))
    if(verbose)message("Proceeding obs... ")
    for (i in incom.obs){
      if(verbose)cat(paste(i," "))
      comvars <-  stats::complete.cases(as.numeric(dat[i,]))
      datpart <- dat[,comvars]
      obs_i <- as.numeric(datpart[i,])
      Sx_i <- Sx[colnames(datpart),colnames(datpart)]
      D2 <- sapply(1:nrow(dat),function(j) tryCatch(stats::mahalanobis(obs_i[!is.na(datpart[j,])],na.omit(as.numeric(datpart[j,])),Sx_i[!is.na(datpart[j,]),!is.na(datpart[j,])]),error=function(e) NA ))
      D2[D2==0] <- NA
      if(any(!is.na(D2))) {
        KNNids <- order(D2,na.last=NA)
        KNNids <- KNNids[sapply(KNNids,function(j) any(!is.na(dat[j,!comvars])))]
      }  else KNNids <- NULL


      dattmp <-  dat
      if(!is.null(KNNids)) KNNids_sel <- KNNids[1:min(K,length(KNNids))]
      if(any(!is.na(D2)) & length(KNNids)>=1) dattmp[i,!comvars] <-
        sapply(1:length(which(!comvars)), function(co) {
          dattmp_all <- dattmp[KNNids,co]
          dattmp_sel <- dattmp[KNNids_sel,co]
          if(any(!is.na(dattmp_sel))) ret <- sum(dattmp_sel*exp(-D2[KNNids_sel]),na.rm=T)/sum(exp(-D2[KNNids_sel])[!is.na(dattmp_sel)],na.rm=T) else ret <- na.omit(dattmp_all)[1]
          ret})

      if(any(is.na(dattmp[i,]))) {        # weil z.B. D2 NAs enth?lt, ODER weil keine gemeinsam obs vars unter den KNN
        still.incom <-  !stats::complete.cases(as.numeric(dattmp[i,]))
        dattmp[i,still.incom] <- apply(dattmp[,still.incom,drop=F],2,mean,na.rm=T)}   # wenn das nicht geht, mean ?ber alle vars

      datimp[i,] <- dattmp[i,]
    }

  } else if("knn.obs.euc" %in% methods) { ## knn.sample.euc: euclidean distance

    incom.obs <- which(apply(dat,1,function(x) any(is.na(x))))
    if(verbose)message(paste0("Number of imcomplete observations: ", length(incom.obs)))
    if(verbose)message("Proceeding obs... ")
    D2 <- as.matrix(stats::dist(scale(dat),upper=T,diag=T)) # eucl
    D2[D2==0] <- NA
    for (i in incom.obs){
      if(verbose)cat(paste(i," "))
      comvars <-  stats::complete.cases(as.numeric(dat[i,]))
      if(any(!is.na(D2[i,]))) {
        KNNids <- order(D2[i,],na.last=NA)
        KNNids <- KNNids[sapply(KNNids,function(j) any(!is.na(dat[j,!comvars])))]
      }  else KNNids <- KNNids_naomit <- NULL


      dattmp <-  dat
      if(!is.null(KNNids)) KNNids_sel <- KNNids[1:min(K,length(KNNids))]
      if(any(!is.na(D2[i,])) & length(KNNids)>=1) dattmp[i,!comvars] <-
        sapply(1:length(which(!comvars)), function(co) {
          dattmp_all <- dattmp[KNNids,co]
          dattmp_sel <- dattmp[KNNids_sel,co]
          if(any(!is.na(dattmp_sel))) ret <- sum(dattmp_sel*exp(-D2[i,KNNids_sel]),na.rm=T)/sum(exp(-D2[i,KNNids_sel])[!is.na(dattmp_sel)],na.rm=T) else ret <- na.omit(dattmp_all)[1]
          ret})

      if(any(is.na(dattmp[i,]))) {        # weil z.B. D2 NAs enth?lt, ODER weil keine gemeinsam obs vars unter den KNN
        still.incom <-  !stats::complete.cases(as.numeric(dattmp[i,]))
        dattmp[i,still.incom] <- apply(dattmp[,still.incom,drop=F],2,mean,na.rm=T)}   # wenn das nicht geht, mean ?ber alle vars

      datimp[i,] <- dattmp[i,]
    }


  } else if("knn.obs.euc.sel" %in% methods) { ## knn.sample.euc.sel KNN per sample extended to cor cutoff (Tutz)  - use only variables with cor>0.2 for distance computation

    incom.obs <- which(apply(dat,1,function(x) any(is.na(x))))
    incom.vars <- which(apply(dat,2,function(x) any(is.na(x))))
    if(verbose)message(paste0("Number of imcomplete observations: ", length(incom.obs)))
    if(verbose)message("Proceeding obs... ")
    Cor <- stats::cor(dat,use="p")
    D2list <- lapply(incom.vars, function(j) {
      varsel <- which(abs(Cor[j,])>cor.var.sel)    # ist j selbst dabei, ist aber ok
      if(length(varsel)>10) varsel <- order(abs(Cor[j,]),decreasing=T)[1:11]
      if(length(varsel)<5) varsel <- order(abs(Cor[j,]),decreasing=T)[1:6]
      D2 <- as.matrix(stats::dist(scale(dat[,varsel])),upper=T,diag=T)
      if(any(is.na(D2))) {
        D2a <- as.matrix(stats::dist(scale(dat)),upper=T,diag=T)*sqrt(length(varsel)/ncol(dat))
        D2[is.na(D2)] <- D2a[is.na(D2)] }
      diag(D2) <- NA
      D2})
    names(D2list) <- incom.vars
    for (i in incom.obs){
      if(verbose)cat(paste(i," "))
      comvars <-  stats::complete.cases(as.numeric(dat[i,]))
      dattmp <-  dat
      for (j in which(!comvars)) {
        D2 <- D2list[[as.character(j)]]
        if(any(!is.na(D2[i,]))) {
          KNNids <- order(D2[i,],na.last=NA)
          KNNids_naomit <- KNNids[sapply(KNNids,function(ii) any(!is.na(dat[ii,j])))]
        }  else KNNids  <- NULL

        # JK, 10/29/19
        # KNNids cannot actually be NULL, this happens if there's a row with all Infs or NA
        if (is.null(KNNids)) stop("Imputation method could not calculate correlations for some samples. Dataset probably contains samples that are all zero or NA and were then logged.")


        if(!is.null(KNNids)) KNNids_sel <- dplyr::intersect(KNNids[1:min(K,length(KNNids))],KNNids_naomit)
        if(length(KNNids_sel)<1) KNNids_sel <- KNNids_naomit[1:min(floor(K/2),length(KNNids_naomit))] else
          if(length(which(sapply(KNNids_sel,function(ii) !is.na(dat[ii,j])))) < floor(K/2) )  KNNids_sel <- KNNids_naomit[1:min(floor(K/2),length(KNNids_naomit))]
          if(any(!is.na(D2[i,])) & length(KNNids)>=1) {
            dattmp_sel <- dattmp[KNNids_sel,j]
            dattmp[i,j] <- sum(dattmp_sel*exp(-D2[i,KNNids_sel]),na.rm=T)/sum(exp(-D2[i,KNNids_sel])[!is.na(dattmp_sel)],na.rm=T) }
      }
      datimp[i,] <- dattmp[i,]
    }


    datimp <- apply(datimp,2, function(x) {
      if(any(is.na(x))) x[is.na(x)] <- mean(x,na.rm=T)
      x})

  }


  datimp
}




#' kNN multi-core internal method.
#'
#' Default settings are from the winning approach of Kiki's paper paper.
#' https://link.springer.com/article/10.1007%2Fs11306-018-1420-2
#' Specifically "knn.obs.euc.sel"
#'
#' This script creates a tmp folder where distance matrices and imputed
#' values are stored
#
#' @param dat A matrix containing metabolites as columns and samples as rows
#' @param cor.var.sel correlation cutoff for metabolites, default: 0.2
#' @param K Number of nearest neighbors to consider, default is 10 (recommendation: do not touch)
#' @param mc_cores Number of cores to use for imputation, dafualt: 5
#' @param verbose T/F, whether to output intermediate steps, default: F
#'
#' @return assay: imputed data
#'
#' @examples
#' # in the context of a SE pipeline
#' ... %>% mt_pre_impute_knn() %>% ...    # standard call
#'
#' @author PG
#'
#' @importFrom S4Vectors metadata
#'
#' @noRd
imputeKNN_multicore <- function(dat,
                                cor.var.sel = 0.2,
                                K=5,
                                mc_cores = 5,
                                verbose = T) {


  datimp <- dat
  # select row_index of missing observation
  incom_obs <- which(apply(dat,1,function(x) any(is.na(x))))

  # select column_index of missing observation
  incom_vars <- which(apply(dat,2,function(x) any(is.na(x)))) %>% names()


  if(verbose)message(paste0("Number of imcomplete observations: ", length(incom_obs)))

  # calculate correlation matrix using pairwise complete observations pairwise.complete.obs
  Cor <- stats::cor(dat,use="p")


  dist_path <- "tmp/distance_matrices/"
  dir.create(dist_path, recursive = TRUE)

  # Precalculate values that will be reused in the mclapply below. Increases speed...
  scale_dat <- scale(dat)
  dist_scale_dat <- stats::dist(scale_dat)

  # Get a list of sample distance matrices based on each incom_vars and varsel.
  # Basically what this is doing is that, if an incom.var is missing, find
  # metabolites with highest correlation, and based on these metabolites
  # calculate which samples are close.
  if(verbose)message("creating distance matrices")
  junk_collector <- parallel::mclapply(incom_vars, function(j) {

    # select variables (metabolites) with > cor.var.sel
    varsel <- which(abs(Cor[j,])>cor.var.sel)    # ist j selbst dabei, ist aber ok

    # If length varsel > 10, take index of the 11 highest ranking variables
    # (metabolites) correlating with j
    if(length(varsel)>10) varsel <- order(abs(Cor[j,]),decreasing=T)[1:11]

    # If length varsel < 5, take index of the 6 highest ranking variables
    # (metabolites) correlating with j
    if(length(varsel)<5) varsel <- order(abs(Cor[j,]),decreasing=T)[1:6]

    # calculate distance matrix of each sample based on scaled varsel data (column-wise)
    # Convert to a full matrix
    D2 <- as.matrix(stats::dist(scale_dat[,varsel]),upper=T,diag=T)

    # if D2 has missing values, replace them by distance matrix of whole dataframe,
    # scaled by length of varsel and ncol(dat)
    if(any(is.na(D2))) {
      D2a <- as.matrix(dist_scale_dat,upper=T,diag=T)*sqrt(length(varsel)/ncol(dat))
      D2[is.na(D2)] <- D2a[is.na(D2)] }
    diag(D2) <- NA
    saveRDS(D2, paste0(dist_path, j, ".rds"), compress = FALSE)
    NULL
  },
  mc.cleanup = TRUE,
  mc.preschedule = TRUE,
  mc.cores = mc_cores)



  impute_path <- "tmp/imputed_samples/"
  dir.create(impute_path, recursive = TRUE)

  if(verbose)message("calculating imputation values")
  junk_collector <-
    parallel::mclapply(incom_obs,
                       function(i) {

                         # for the observation/sample that has a missing metabolite,
                         # select all metabolites that have values
                         incomvars <-  dat[i, is.na(dat[i,])] %>% names()

                         dattmp <-  dat

                         # iterate over te missing metabolites
                         for (j in incomvars) {

                           # select sample distance matrix calculated for missing metabolite
                           D2_path <- paste0(dist_path, j, ".rds")
                           D2 <- readRDS(D2_path)

                           if(any(!is.na(D2[i,]))) {

                             # order samples by decreasing distance
                             KNNids <- order(D2[i,],na.last=NA)

                             # omit neighbors that also have the same metabolite j missing
                             KNNids_naomit <- KNNids[sapply(KNNids,function(ii) any(!is.na(dat[ii,j])))]

                           }  else KNNids  <- NULL

                           # Select at most 10 nearest neighbors which overlap with KNNids_naomit
                           if(!is.null(KNNids)) KNNids_sel <- dplyr::intersect(KNNids[1:min(K,length(KNNids))],KNNids_naomit)

                           # If there really are no KNNids_sel, get at most floor(K/2) metbaolites from KNNids_naomit
                           # OR if there < floor(K/2) KNNids_sel metabolites, get at most floor(K/2) metbaolites
                           # from KNNids_naomit
                           if(length(KNNids_sel)<1) KNNids_sel <- KNNids_naomit[1:min(floor(K/2),length(KNNids_naomit))] else
                             if(length(which(sapply(KNNids_sel,function(ii) !is.na(dat[ii,j])))) < floor(K/2) )  KNNids_sel <- KNNids_naomit[1:min(floor(K/2),length(KNNids_naomit))]

                           # calculate the imputation value
                           if(any(!is.na(D2[i,])) & length(KNNids)>=1) {
                             dattmp_sel <- dattmp[KNNids_sel,j]
                             dattmp[i,j] <- sum(dattmp_sel*exp(-D2[i,KNNids_sel]),na.rm=T)/sum(exp(-D2[i,KNNids_sel])[!is.na(dattmp_sel)],na.rm=T) }
                         }

                         saveRDS(dattmp[i,], paste0(impute_path, i, ".rds"), compress = FALSE)

                         NULL
                       },
                       mc.cleanup = TRUE,
                       mc.preschedule = TRUE,
                       mc.cores = mc_cores)




  # Fill in imputed values

  for (i in incom_obs){

    datimp[i,] <- readRDS(paste0(impute_path, i, ".rds"))

  }


  # if there are still values missing, just assign mean value
  datimp <- apply(datimp,2, function(x) {
    if(any(is.na(x))) x[is.na(x)] <- mean(x,na.rm=T)
    x})


  datimp
}
