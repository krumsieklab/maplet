# Univariate Analysis - Binary Outcome
# This example illustrates one of the simplest and most commonly performed analyses supported by
# the maplet package - a univariate association analysis of a binary outcome. The sections below
# illustrates the typical steps taken from loading the data, through visualization. In this
# example the binary outcome we will be looking at is the sample annotation, Diagnosis.

library(maplet)

file_data <- system.file("extdata", "example_data/simulated_data.xlsx", package = "maplet")
file_output <- "maplet_preprocessing_example.html"

# Load Data ----
D <-
  # load data - this function loads the assay data only
  mt_load_xls(file=file_data, sheet="data", samples_in_row=T, id_col="sample") %>%
  # load metabolite (rowData) annotations
  mt_anno_xls(file=file_data, sheet="metinfo",anno_type="features", anno_id_col="name", data_id_col = "name") %>%
  # load clinical (colData) annotations
  mt_anno_xls(file=file_data, sheet="clin", anno_type="samples", anno_id_col ="sample", data_id_col ="sample") %>%
  # record raw data info
  mt_reporting_data()

# Data Cleaning ----
D %<>%
  mt_reporting_heading(heading = "Pre-processing", lvl=1) %>%
  mt_reporting_heading(heading = "Data Clean-up", lvl = 2) %>%
  # remove samples with missing values for outcome Diagnosis
  mt_modify_filter_samples(filter = !is.na(Diagnosis)) %>%
  # ensure annotaiton column 'Diagnosis' is of type factor
  mt_anno_apply(anno_type = "samples", col_name = "Diagnosis", fun = as.factor) %>%
  # record data info after filtering
  mt_reporting_data()


# Filter Missing Values ----
D %<>%
  mt_reporting_heading(heading = "Filtering", lvl = 2) %>%
  # section text
  mt_reporting_text(text = "Plot percent missingness for each metabolite before filtering, filter out metabolites with >= 50%
                    missingness, plot percent missingness for each metabolite after filtering, add missingness annotation
                    columns to both metabolite and sample annotation data frames.") %>%
  # plot missingness distribution
  mt_plots_missingness(feat_max=0.4) %>%
  # filter metabolites with more than 40% missing values per group
  mt_pre_filter_missingness(feat_max = 0.4) %>%
  # plot missingness distribution after filtering
  mt_plots_missingness(feat_max=0.4) %>%
  # add missingness percentage as annotation to samples (remaining missing)
  mt_anno_missingness(anno_type = "samples", out_col = "missing") %>%
  # add missingness percentage as annotation to metabolites
  mt_anno_missingness(anno_type = "features", out_col = "missing")

# Normalization ----
D %<>%
  mt_reporting_heading(heading = "Normalization", lvl = 2) %>%
  # pre-normalization sample boxplot
  mt_plots_sample_boxplot(title = "Before normalization") %>%
  # normalize abundances using probabilistic quotient
  mt_pre_norm_quot(feat_max = 0.2, ref_samples = Diagnosis==0) %>%
  # post-normalization sample boxplot
  mt_plots_sample_boxplot(title = "After normalization") %>%
  # dilution plot showing dilution factors from quotient normalization
  mt_plots_dilution_factor(boxpl = T, in_col = 'Diagnosis')

# Data Transfomration ----
D %<>%
  mt_reporting_heading(heading = "Data Transformation", lvl=2) %>%
  # Log2 transformation
  mt_pre_trans_log() %>%
  # scale and center data
  mt_pre_trans_scale()

# Imputation ----
D %<>%
    mt_reporting_heading(heading = "Imputation", lvl = 2) %>%
    # pre-imputation sample boxplot
    mt_plots_sample_boxplot(title = "Before imputation") %>%
    mt_pre_impute_knn() %>%
    # post-imputation sample boxplot
    mt_plots_sample_boxplot(title = "After imputation")

# Global Data Overview ----
D %<>% mt_reporting_heading(heading = "Global Statistics", lvl = 1) %>%
    # plot PCA
    mt_plots_pca(scale_data = T, title = "scaled PCA - Diagnosis", color=Diagnosis, size=2.5, ggadd=scale_size_identity()) %>%
    # plot UMAP
    mt_plots_umap(scale_data = T, title = "scaled UMAP - Diagnosis", color=Diagnosis, size=2.5, ggadd=scale_size_identity()) %>%
    # plot heatmap
    mt_plots_heatmap(scale_data = T, annotation_col = c("Diagnosis"), annotation_row = c("SUPER_PATHWAY"),
                     clustering_method = "ward.D2", fontsize = 5, cutree_rows = 3, cutree_cols = 3,
                     color=gplots::bluered(101))

D %>% mt_reporting_html(file=file_output, output_calls = T)
