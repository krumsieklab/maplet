#' Run Cross-Validation / Sampling Methods
#'
#' Performs one of the following cross-validation / sampling methods for a provided machine learning
#' algorithm: k-fold cross-validation, (TO-DO: extend to more algorithms later).
#'
#' @param D \code{SummarizedExperiment} input.
#' @param ml_fun The maplet-implemented algorithm to use. See mtml_algorithms for a list of implemented algorithms.
#' @param ml_name Name under which the results will be stored.
#' @param response_col Name of the colData column containing the data labels.
#' @param response_type Indicates the type of learning task. One of "binary" or
#'    "continuous". Default: "binary".
#' @param sampling_method Type of cross-validation / sampling method using. Must be one of: cv (add more later).
#'    Default: cv.
#' @param num_folds For k-fold cross-validation, this represents the number of folds (k).
#' @param rand_seed Value to set random seed.
#' @param mod_args Named list of parameters passed to maplet-implemented machine learning algorithm. See algorithm
#'    implementations for default parameters used. Default: empty list.
#' @param pred_args Named list of parameters passed to predict() for a specific machine learning function. For example,
#'    predict.glmnet accepts the parameter s, which supplies the penalty parameter lambda at which prediction are
#'    required. See details of each function in mtml_algorithms for parameter lists for each algorithm-specific
#'    predict function. See mpalet-implemented ml functions for default parameters used. Default: empty list.
#' @param pos_class Required if response_type = "binary". Indicates which value in the response column
#'    represents the positive class.
#'
#' @return result$output: List containing a list of results for each fold.
#'
#' @author KC
#'
#' @export
mt_ml_repeat <- function(D,
                         ml_fun,
                         ml_name,
                         response_col,
                         response_type="binary",
                         sampling_method = "cv",
                         num_folds,
                         rand_seed,
                         mod_args = list(),
                         pred_args = list(),
                         pos_class){

  # Validate arguments
  ml_fun_name <- deparse(substitute(ml_fun))
  required <- c("D", "ml_fun", "ml_name", "response_col")
  passed <- names(as.list(match.call())[-1])
  if(all(required %in% passed)==F) stop(paste0("Missing one or more required arguments: "),
                                        paste0(required, collapse = ", "))
  if("SummarizedExperiment" %in% class(D) == F) stop("D must be a SummarizedExperiment object.")
  if(response_type == "binary" && missing(pos_class)) stop("A value must be provided for 'pos_class' when response_type is 'binary'.")
  # if(find(ml_fun) != "package:maplet") stop(glue::glue("Could not find {ml_fun} in package maplet."))
  # NOTE TO KELSEY: Add checks later
  # - add check ml_name not used
  # - add check response_col is column in colData
  # - add check sampling_method is valid

  x <- t(assay(D))   # data frame to be split into training / testing
  y <- colData(D)[,response_col]   # corresponding class labels
  if (missing(rand_seed)) rand_seed <- NULL

  # for binary case, ensure y is of type factor and pos_class is the first level
  if(response_type=="binary"){
    y_vals = unique(y)
    pos_idx = which(y_vals==pos_class)

    if(length(pos_idx)==0) stop("The provided value for 'pos_class' was not found in the response column.")
    if(length(y_vals)>2) stop("'response_type' is 'binary' but 'response_col' has more than two values.")

    y <- factor(y, levels=c(y_vals[pos_idx], y_vals[-pos_idx]))
  }else{
    pos_class = NULL
  }

  # --- Choose Sampling Method --- #
  # k-fold cross-validation is only one implemented right now; will add others later

  # perform classic k-fold cross-validation
  if(sampling_method == "cv"){
    if(missing(num_folds)) stop("When performing k-fold cross-validation, a value for argument num_folds must be provided.")
    test_idx <- make_cv_partitions(x, y, num_folds, rand_seed, response_type)
  }

  # call run_train_test, wrapper function for training and testing
  fold_res_list <- lapply(test_idx, run_train_test, ml_fun_name, x, y, response_type, mod_args, pred_args)
  # flip fold list inside-out to separate predictions and indices
  pred_test_idx_list <- purrr::transpose(fold_res_list)

  # add status information
  funargs <- maplet:::mti_funargs()
  D %<>%
    maplet:::mti_generate_result(
      funargs = funargs,
      logtxt = paste(glue::glue("Performed {ml_fun_name}, repeated using {sampling_method} method.")),
      output = list(name = ml_name,
                    response_col = response_col,
                    response_type = response_type,
                    pos_class = pos_class,
                    num_folds = num_folds,
                    pred = pred_test_idx_list$pred,
                    test_idx = pred_test_idx_list$test_idx,
                    object = pred_test_idx_list$object)
    )
  # return SE
  D

}

# Training / Testing function
# This function for subsets the data frames / labels, calls the implemented algorithms, and
#    gets the predicted values.
run_train_test <- function(idx, ml_fun_name, x, y, response_type, mod_args, pred_args){

  # subset training / testing data and training labels
  train_data <- x[-idx,]
  train_label <- y[-idx]
  test_data <- x[idx,]

  # retrieve internal mtml algorithm function
  ml_fun <- utils::getFromNamespace(ml_fun_name, "maplet")

  # train model
  ml_fun_args <- list(train_data = train_data, train_label = train_label, test_data = test_data,
                      response_type=response_type, mod_args = mod_args, pred_args = pred_args)
  pred_args <- do.call(ml_fun, ml_fun_args)

  # get predicted values
  if(is.null(pred_args$type)) pred_args$type <- "response" # get class probabilities by default
  pred <- do.call(predict, pred_args)

  fold_res <- list(test_idx = idx, pred = as.numeric(pred))

  fold_res

}


# Cross-Validation / Sampling Functions
# k-fold cross-validation
make_cv_partitions <- function(x, y, num_folds, rand_seed, response_type){

  # set random seed
  if (!is.null(rand_seed)) set.seed(rand_seed)

  # create folds
  if(response_type=="binary"){
    idx_lists <- caret::createFolds(y, k=num_folds)
  }else{
    idx_lists <- create_fold_splits(y, nfolds=num_folds)
  }

  idx_lists

}

create_fold_splits <- function(y, nfolds){
  df <- cbind.data.frame(index = 1:length(y), value=y)
  df <- df[order(df$value),]
  d <- df$index

  bins <- split(d, ceiling(seq_along(d)/nfolds))
  nbins = length(bins)
  fold_list <- lapply(bins, function(x){
    if(length(x)==nfolds){
      cbind(index=x, fold=sample(seq(1,nfolds)))
    }else{
      cbind(index=x, fold=sample(seq(1,nfolds))[1:length(x)])
    }
  })
  folds <- do.call(rbind.data.frame,fold_list)
  folds <- folds[order(folds$index),] %>% split(., .$fold) %>% lapply(., function(x){x[,"index"]})

  folds
}

# TO-DO: Add additional functions for other cross-validation methods


