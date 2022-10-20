### -- Example: Average Duplicate Samples & features -- ###
# This script demonstrates how to use the maplet machine learning framework. The framework
# currently consists of the following functions: mt_ml_repeat, mt_plots_evaluate_ml, and the
# internal algorithm functions designated by the prefix 'mtml_'. Machine learning algorithms
# are not called by the user directly. Instead, the selected 'mtml' algorithm is passed as an
# argument to mt_ml_repeat. This flexible framework simplifies running a machine learning analysis
# by containing the execution of the analysis to a single function, while also allowing for the
# continuous development and inclusion of new machine learning algorithms without disrupting the
# existing framework.
# The framework currently supports two types of outcomes: binary (classification) and continuous
# (regression). The framework also allows for the inclusion of covariates of any data type. How to
# run each outcome type and analyses with covariates are demonstrated below.

library(maplet)
library(tidyverse)
library(readxl)

### ---------------- DOWNLOAD AND PREPROCESS DATA ---------------- ###
# download data file
url <- 'https://ndownloader.figshare.com/files/10531342'
localfile <- 'qmdiab_preprocessed.xlsx'
# download and load
download.file(url, localfile)

# read assay
data <- read_excel(localfile, sheet='plasma') %>%
  dplyr::select(-one_of("QMDiab-ID")) %>%
  dplyr::select(-AGE, -GENDER, -BMI, -ETH, -T2D) # only keep metabolites
# read clinical variables, and ensure data.frame is maplet-compatible
clin <- read_excel(localfile, sheet='plasma') %>%
  # only keep clinical variables
  dplyr::select(AGE, GENDER, BMI, ETH, T2D, one_of("QMDiab-ID"))
# rename, not very clean, can't figure out tidyverse-way right now
clin$QMDiab_ID <- clin$`QMDiab-ID`
clin <- clin[,-6]

# read metabolite information, and ensure data.frame is maplet-compatible
metinfo <- read_excel(localfile, sheet='plasma annotations') %>%
  mutate(name = BIOCHEMICAL)

# assemble into SE, and make maplet compatible
D <- SummarizedExperiment(assay=t(data), colData = clin, rowData = metinfo) %>%
  mt_clean_validate_se() %>%
  mt_pre_impute_min()

### ---------------- BINARY OUTCOME ---------------- ###
D1 <- D
D1 %<>% mt_reporting_heading("Binary Outcome") %>%
  mt_reporting_heading("glmnet", lvl=2) %>%
  mt_ml_repeat(ml_fun = mtml_glmnet,
               ml_name = "binary_lasso",
               mod_args = list(alpha = 0.5),
               response_col = "T2D",
               num_folds = 5,
               rand_seed = 225,
               pred_args = list(s = "lambda.min"),
               pos_class=1) %>%
  mt_plots_ml_evaluate(ml_name = "binary_lasso") %>%
  mt_reporting_heading("xgboost", lvl=2) %>%
  mt_ml_repeat(ml_fun = mtml_gradboost,
               ml_name = "binary_xgboost",
               response_col = "T2D",
               num_folds = 5,
               rand_seed = 255,
               pos_class=1) %>%
  mt_plots_ml_evaluate(ml_name = "binary_xgboost")
D1 %>% mt_reporting_html(file = "mt_ml_framework_binary.html",
                         title = "Machine Learning Framework - Binary Example", output_calls = T)


### ---------------- CONTINUOUS OUTCOME ---------------- ###
D2 <- D
D2 %<>% mt_reporting_heading("Continuous Outcome") %>%
  mt_reporting_heading("glmnet", lvl=2) %>%
  mt_ml_repeat(ml_fun = mtml_glmnet,
               ml_name = "continuous_lasso",
               mod_args = list(alpha = 0.5),
               response_col = "BMI",
               response_type = "continuous",
               num_folds = 5,
               rand_seed = 225,
               pred_args = list(s = "lambda.min")) %>%
  mt_plots_ml_evaluate(ml_name = "continuous_lasso") %>%
  mt_reporting_heading("xgboost", lvl=2) %>%
  mt_ml_repeat(ml_fun = mtml_gradboost,
               ml_name = "continuous_xgboost",
               mod_args = list(alpha = 0.5),
               response_col = "BMI",
               response_type = "continuous",
               num_folds = 5,
               rand_seed = 225) %>%
  mt_plots_ml_evaluate(ml_name = "continuous_xgboost")
D2 %>% mt_reporting_html(file = "mt_ml_framework_continuous.html",
                         title = "Machine Learning Framework - Continuous Example", output_calls = T)


### ---------------- WITH COVARIATES ---------------- ###
D3 <- D
D3 %<>% mt_reporting_heading("Covariates Included") %>%
  mt_reporting_heading("Continous outcome - glmnet", lvl=2) %>%
  mt_ml_repeat(ml_fun = mtml_glmnet,
               ml_name = "continuous_lasso",
               mod_args = list(alpha = 0.5),
               response_col = "BMI",
               response_type = "continuous",
               covar_cols = c("T2D"),
               num_folds = 5,
               rand_seed = 225,
               pred_args = list(s = "lambda.min")) %>%
  mt_plots_ml_evaluate(ml_name = "continuous_lasso") %>%
  mt_reporting_heading("binary outcome - xgboost", lvl=2) %>%
  mt_ml_repeat(ml_fun = mtml_gradboost,
               ml_name = "binary_xgboost",
               response_col = "T2D",
               covar_cols = c("BMI", "AGE"),
               num_folds = 5,
               rand_seed = 255,
               pos_class=1) %>%
  mt_plots_ml_evaluate(ml_name = "binary_xgboost")
D3 %>% mt_reporting_html(file = "mt_ml_framework_covariates.html",
                         title = "Machine Learning Framework - Covariates Example", output_calls = T)