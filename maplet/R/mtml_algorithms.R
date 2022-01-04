#' Wrapper function for glmnet::cv.glment
#'
#' Internal maplet wrapper for glment::cv.glmnet. This function called by one of the control functions for
#' machine learning in maplet (e.g. cv) using the wrapper function mt_ml_repeat.
#' 
#' @param train_data Input data frame of dimension nobs x nvars.
#' @param train_label Response variable taken from the colData data frame.
#' @param mod_args List of additional arguments to pass to glmnet::cv.glmnet.
#' @param pred_args List of additional arguments to be passed to predict. If NULL, default arguments will be
#'    returned.
#' 
#' @return List with trained model and predict args.
#' 
#' @author KC
#' 
#' @noRd
mtml_glmnet <- function(train_data, train_label, test_data, mod_args, pred_args){
  
  # get training partitions for x and y
  mod_args$x <- as.matrix(train_data)
  mod_args$y <- train_label
  
  # set alpha if missing and set family to binomial
  if(is.null(mod_args$alpha)) mod_args$alpha <- 0.5
  mod_args$family <- "binomial"   # hard-coded for now; may want to give user options in the future
  
  # train model
  mod <- do.call(glmnet::cv.glmnet, mod_args)
  
  # if predict args missing, use default s
  if(is.null(pred_args)) pred_args = list(s = "lambda.min")
  
  # add model and test data to predict args
  pred_args$object <- mod
  pred_args$newx <- test_data
  
  # return predict args
  pred_args
  
}

#' Wrapper function for gradient boost algorithm
#'
#' Internal maplet wrapper for gradient boost algorithms from the package xgboost. This function is called by 
#' one of the control functions for machine learning in maplet (e.g. mt_ml_cv). Default hyperparameters are used
#' for XGBoost model with the following exceptions: (1) eta=0.1, (2) use cv.xgboost to find best n_rounds.
#' Random grid search hyper-parameter tuning [will be available] in a separate implementation if performance is
#' not satisfactory. Parameter tuning method as recommended by: 
#' https://www.hackerearth.com/practice/machine-learning/machine-learning-algorithms/beginners-tutorial-on-xgboost-parameter-tuning-r/tutorial/.
#'
#' @param train_data 
#' @param train_label
#' @param mod_args 
#' @param pred_args
#' 
#' @return List with trained model and predict args.
#' 
#' @author KC
#' 
#' @noRd
mtml_gradboost <- function(train_data, train_label, test_data, mod_args, pred_args){
  
  # default parameters for xgboost
  DEFAULT_PARAMETERS <- list(booster = "gbtree", objective = "binary:logistic", eta=0.1, gamma=0, 
                             max_depth=6, min_child_weight=1, subsample=1, colsample_bytree=1, eval_metric="error")
  
  # get training partitions for x and y
  train_matrix <- train_data %>% data.matrix()
  train_labels <- as.numeric(train_label)-1
  mod_args$data <- xgboost::xgb.DMatrix(data = train_matrix, label = train_labels)
  
  # use cv.xgboost to find best value for parameter n_rounds
  xgbcv <- xgboost::xgb.cv( params = DEFAULT_PARAMETERS, data = mod_args$data, nrounds = 100, nfold = 5, 
                   showsd = T, early_stopping_rounds = 20, stratified = T, maximize = F, verbose = F)
  
  if(is.null(mod_args$nrounds)) mod_args$nrounds <- xgbcv$best_iteration
  if(is.null(mod_args$params)) mod_args$params <- DEFAULT_PARAMETERS
  
  # train model
  mod <- do.call(xgboost::xgb.train, mod_args)
  
  # add model and formatted test data to predict args
  pred_args$object <- mod
  pred_args$newdata <- xgboost::xgb.DMatrix(data = test_data)
  
  # return return predict args
  pred_args
  
}
