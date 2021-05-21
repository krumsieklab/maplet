# ------------ maplet Example Pipeline ------------
# This example pipeline demonstrates the use of 79 of the 86 functions available from maplet.
# Users can use `?` to find additional information on each function. Users can also review the User Manual available on
# https://github.com/krumsieklab/maplet/blob/main/guide/maplet_Reference_Guide_markdown.md.
#
# NOTE: This example is divided into sections. Use SHIFT+CTRL+O to see the document outline.
#
# NOTE: One of the main features of maplet is the ability to easily group reusable function steps into meta-functions.
#   However, for the sake of readability, we present the steps of this example pipeline in a linear manner.
#   To see how the maplet can take advantage of these types of meta-functions, see {EXAMPLE_WITH_META_FUNCTIONS}.


library(maplet)
library(tidyverse)

purrr::zap()

file_data <- system.file("extdata", "example_data/simulated_data.xlsx", package = "maplet")

# PART 1 - STARTING A METABOTOOLS PIPELINE ----------------------------------------------------

D <-
  # validate checksum
  mt_load_checksum(file=file_data, checksum = "80afcd72481c6cf3dcf83342e3513699") %>%
  # load data - this function loads the assay data only
  #   alternative loading functions: mt_load_metabolon_v1(), mt_load_metabolon_v2(), mt_load_metabolon_lipidomics(),
  #     mt_load_olink(), mt_load_ucd(), mt_load_wcm(), mt_load_nightingale, mt_load_metabolon_new_format()
  mt_load_xls(file=file_data, sheet="data", samples_in_row=T, id_col="sample") %>%
  # load metabolite (rowData) annotations
  mt_anno_xls(file=file_data, sheet="metinfo",anno_type="features", anno_id_col="name", data_id_col = "name") %>%
  # load clinical (colData) annotations
  mt_anno_xls(file=file_data, sheet="clin", anno_type="samples", anno_id_col ="sample", data_id_col ="sample") %>%
  # # log assay dimensions and number of columns for both metabolite and clincial annotations
  mt_reporting_data() %>%
  # start timing
  mt_reporting_tic() %>%
  {.}
# additional functions used at beginning of pipelines:
#   - mt_settings - set global settings for maplet pipeline
#   - mt_load_flag_logged - for flagging a loaded dataset as already log transformed

# PART 2 - DATA CLEANING ----------------------------------------------------

D <- D %>%
  # heading
  mt_reporting_heading(heading = "Data Clean-up", lvl = 1) %>%
  # section text
  mt_reporting_text(text = "Filter samples that are missing values for Diagnosis,add sample annotation column with log10 of
                    PreBioPSA, convert sample annotaiton column Diagnosis to factors,filter metabolites that are missing values
                    for SUB_PATHWAY, log dataset information for this point of the pipeline.") %>%
  # filter samples
  mt_modify_filter_samples(filter = !is.na(Diagnosis)) %>%
  # create additional variable
  mt_anno_mutate(anno_type = "samples", col_name = "PreBioPSALog", term = log10(PreBioPSA)) %>%
  # modify variable to factor
  mt_anno_apply(anno_type = "samples", col_name = "Diagnosis", fun = as.factor) %>%
  # remove metabolites with no pathway annotation
  mt_modify_filter_features(filter = !is.na(SUB_PATHWAY)) %>%
  # log assay dimensions and number of columns for both metabolite and clinical annotations
  mt_reporting_data() %>%
  {.}
# additional data cleaning function:
#   - mt_pre_zero_to_na - for platforms that represent sub-LOD/missing values as zeros


# PART 3.1 - PREPROCESSING: FILTERING MISSING VALUES ----------------------------------------------------

D <- D %>%
  # heading for html file
  mt_reporting_heading(heading = "Preprocessing", lvl=1) %>%
  # heading for html file
  mt_reporting_heading(heading = "Filtering", lvl = 2) %>%
  # section text
  mt_reporting_text(text = "Plot percent missingness for each metabolite before filtering, filter out metabolites with >= 50%
                    missingness, plot percent missingness for each metabolite after filtering, add missingness annotation
                    columns to both metabolite and sample annotation data frames.") %>%
  # plot missingness distribution
  mt_plots_missingness(feat_max=0.5) %>%
  # filter metabolites with more than 50% missing values per group
  mt_pre_filter_missingness(feat_max = 0.5, group_col = "Diagnosis") %>%
  # plot missingness distribution after filtering
  mt_plots_missingness(feat_max=0.5) %>%
  # add missingness percentage as annotation to samples (remaining missing)
  mt_anno_missingness(anno_type = "samples", out_col = "missing") %>%
  # add missingness percentage as annotation to metabolites
  mt_anno_missingness(anno_type = "features", out_col = "missing") %>%
  {.}


# PART 3.2 - PREPROCESSING: NORMALIZATION ----------------------------------------------------

D <- D %>%
  # heading for html file
  mt_reporting_heading(heading = "Normalization", lvl = 2) %>%
  # section text
  mt_reporting_text(text = "Plot sample boxplots before normalization, apply median batch correction, perform quotient
                    normalization, plot boxplot with dilution factors from quotient normalization, plot sample boxplot after
                    normalization, log transform the data, impute missing data using knn, plot sample boxplot after imputation,
                    detect outliers, log dataset info, write pre-processed data to file.") %>%
  # plot sample boxplots
  mt_plots_sample_boxplot(color=Diagnosis, title = "Original", plot_logged = T) %>%
  # apply batch correction
  #   alternative batch correction function: mt_pre_batch_combat
  mt_pre_batch_median(batch_col = "BOX.NUMBER") %>%
  # plot sample boxplots after batch correction
  mt_plots_sample_boxplot(color=Diagnosis, title = "After batch correction", plot_logged = T) %>%
  # normalize abundances using probabilistic quotient
  #   alternative normalization function: mt_pre_norm_external
  mt_pre_norm_quot(feat_max = 0.2, ref_samples = Diagnosis==0) %>%
  # show dilution plot
  mt_plots_dilution_factor(in_col="Diagnosis") %>%
  # plot sample boxplots after normalization
  mt_plots_sample_boxplot(color=Diagnosis, title = "After normalization", plot_logged = T) %>%
  # log transform
  #   other data transformation functions: mt_pre_trans_exp, mt_pre_trans_relative, mt_pre_trans_scale
  mt_pre_trans_log() %>%
  # impute missing values using knn
  #   alternative imputation functions: mt_pre_impute_min
  mt_pre_impute_knn() %>%
  # plot sample boxplot after imputation
  mt_plots_sample_boxplot(color=Diagnosis, title = "After imputation", plot_logged = T) %>%
  # outlier detection (univariate)
  #   alternative functions: mt_pre_outlier_detection_mahalanobis(), mt_pre_outlier_detection_leverage()
  #   related function: mt_pre_outlier_to_na()
  mt_pre_outlier_detection_univariate() %>%
  # print infos about dataset
  mt_reporting_data() %>%
  # write preprocessed data to Excel file
  #   other writing functions: mt_write_se_rds (save SummarizedExerpiment object)
  mt_write_se_xls(file = "PreprocessedData.xlsx") %>%
  {.}

# Additional pre-processing functions
#   - mt_pre_confounding_correction() - function for correcting confounding variables
#   - mt_pre_confounding_correction_stepwise_aic() - an alterenative function for correcting confounders that uses stepwise aic
# NOTES ON BEST PRACTICES: If incorporated in this pipeline, the above functions would correct for the variable age such that
#     none of the following functions have to take care of those confounders anymore. If this function is included, confounders
#     should not be included in any of the following functions. It is generally agreed that including the confounders in the
#     linear models themselves is preferable to pre-correction.

# PART 4 - GET PATHWAY ANNOTATIONS ----------------------------------------------------

D <- D %>%
  # heading for html file
  mt_reporting_heading(heading = "Get Pathway Annotations", lvl = 1) %>%
  # get KEGG ids from HMDB ids
  mt_anno_hmdb_to_kegg(in_col = "HMDb", out_col = "KEGG_ids") %>%
  # get pathway annotations
  #   alternative functions: mt_anno_pathways_xls, mt_anno_pathways_graphite, mt_anno_pathways_uniprot
  mt_anno_pathways_hmdb(in_col = "HMDb", out_col = "pathway", pwdb_name = "KEGG") %>%
  # remove redundant
  mt_anno_pathways_remove_redundant(feat_col = "KEGG_ids", pw_col = "pathway") %>%
  # write pathway annotations
  mt_write_pathways(file="ExamplePipeline_PathwayAnnotations.xlsx", pw_col = "pathway") %>%
  {.}


# PART 5 - GLOBAL STATISTICS ----------------------------------------------------

D <- D %>%
  # heading for html file
  mt_reporting_heading(heading = "Global Statistics", lvl = 1) %>%
  # plot PCA
  mt_plots_pca(scale_data = T, title = "scaled PCA - Diagnosis", color=Diagnosis, size=2.5, ggadd=scale_size_identity()) %>%
  # plot UMAP
  mt_plots_umap(scale_data = T, title = "scaled UMAP - Diagnosis", color=Diagnosis, size=2.5, ggadd=scale_size_identity()) %>%
  # plot heatmap
  mt_plots_heatmap(scale_data = T, annotation_col = c("Diagnosis"), annotation_row = c("SUPER_PATHWAY"),
                   clustering_method = "ward.D2", fontsize = 5, cutree_rows = 3, cutree_cols = 3, color=gplots::bluered(101)) %>%
                    {.}


# PART 6.1 - STATISTICAL ANALYSIS, OUTCOME: DIAGNOSIS, METHOD: MISSINGNESS ANALYSIS ---------------------------------------

#create another SE object for first analysis branch (missingness & metabolites)
D1 <- D

D1 <- D1 %>%
  # heading for html file
  mt_reporting_heading(heading = "Missingness analysis", lvl = 1) %>%
  # section text
  mt_reporting_text(text = "Perform missingness analysis to determine if NAs significantly accumulate in one of the Diagnosis
                    groups. Adjust output of test using multiple testing correction.") %>%
  # compute Fisher's exact test
  mt_stats_univ_missingness(in_col="Diagnosis", stat_name="missingness") %>%
  # create p-value qq plot
  mt_plots_pval_qq(stat_name = "missingness") %>%
  # apply multiple testing correction
  #   alternative function: mt_post_multtest_effdim
  mt_post_multtest(stat_name="missingness", method="BH") %>%
  {.}


# PART 6.2 - STATISTICAL ANALYSIS, OUTCOME: AGE, METHOD: LINEAR REGRESSION -----------------------------------------------

D1 <- D1 %>%
  # heading for html file
  mt_reporting_heading(heading = "Statistical Analysis", lvl = 1) %>%
  # heading for html file
  mt_reporting_heading(heading = "Age analysis", lvl = 2) %>%
  # Pearson correlation
  mt_stats_univ_cor(method = "pearson",
                    in_col = "Age",
                    stat_name = "Age met")%>%
  # create p-value qq plot
  mt_plots_pval_qq(stat_name = "Age met") %>%
  # add multiple testing correction
  mt_post_multtest(stat_name = "Age met", method = "BH") %>%
  # add stats logging
  mt_reporting_stats(stat_name = "Age met", stat_filter = p.adj < 0.05) %>%
  # write statistical results to file
  mt_write_stats(file = "AgeAnalysis.xlsx", feat_col = "BOCHEMICAL") %>%
  # pvalue histogram
  mt_plots_pval_hist(stat_list = "Age met") %>%
  # volcano plot as overview of results
  mt_plots_volcano(stat_name = "Age met",
                   x = statistic,
                   feat_filter = p.adj < 0.05,
                   color = p.adj < 0.05) %>%
  # scatter plot
  mt_plots_box_scatter(stat_name = "Age met",
                       x = Age,
                       plot_type = "scatter",
                       feat_filter = p.adj < 1E-10, # made very small because otherwise would take an eternity to generate all plots
                       feat_sort = p.value,
                       annotation = "{sprintf('P-value: %.2e', p.value)}\nP.adj: {sprintf('%.2e', p.adj)}") %>%
  {.}


# PART 6.3 - STATISTICAL ANALYSIS, OUTCOME: DIAGNOSIS, METHOD: LINEAR REGRESSION (t-test) ----------------------------

D1 <- D1 %>%
  # heading for html file
  mt_reporting_heading(heading = "Diagnosis analysis", lvl = 2) %>%
  # linear model for binary function (equivalent to t-test)
  #   alternative functions: mt_stats_univ_wilcox, mt_stats_univ_lm_matrixeqtl
  mt_stats_univ_lm(formula = ~ Diagnosis,
                   stat_name = "Diagnosis met") %>%
  # add fold change
  mt_post_fold_change(stat_name = "Diagnosis met") %>%
  # add multiple testing correction
  mt_post_multtest(stat_name = "Diagnosis met", method = "BH") %>%
  # add stats logging
  mt_reporting_stats(stat_name = "Diagnosis met", stat_filter = p.adj < 0.05) %>%
  # write statistical results to file
  mt_write_stats(file = "DiagnosisAnalysis.xlsx") %>%
  # pvalue histogram
  mt_plots_pval_hist(stat_list = "Diagnosis met") %>%
  # volcano plot as overview of results
  mt_plots_volcano(stat_name = "Diagnosis met",
                   x = fc,
                   feat_filter = p.adj < 0.05,
                   colour       = p.adj < 0.05) %>%
  # boxplot
  mt_plots_box_scatter(stat_name          ="Diagnosis met",
                       x                  = Diagnosis,
                       fill               = Diagnosis,
                       plot_type          = "box",
                       feat_filter       = p.adj < 0.05,
                       feat_sort         = p.value,
                       annotation         = "{sprintf('P-value: %.2e', p.value)}\nPadj: {sprintf('%.2e', p.adj)}") %>%
                       {.}


# PART 7.1 - STATISTICAL RESULTS PRESENTATION: STATS PATHWAY BAR PLOT & PATHVIEW ----------------------------------------------------

D1 <- D1 %>%
  # heading for html file
  mt_reporting_heading(heading = "Statisitcal Results Presentation", lvl = 1) %>%
  # heading for html file
  mt_reporting_heading(heading = "Barplot", lvl = 2) %>%
  # create statsbarplots
  mt_plots_stats_pathway_bar(stat_list = c("Age met", "Diagnosis met"),
                     feat_filter = p.adj < 0.05,
                     group_col = "SUB_PATHWAY",
                     color_col = "SUPER_PATHWAY",
                     y_scale = "count",
                     sort_by_y = T,
                     assoc_sign_col = "statistic") %>%
  mt_plots_stats_pathway_bar(stat_list = c("Age met", "Diagnosis met"),
                     feat_filter = p.adj < 0.05,
                     group_col = "SUB_PATHWAY",
                     color_col = "SUPER_PATHWAY",
                     y_scale = "fraction",
                     sort_by_y = T,
                     assoc_sign_col = "statistic") %>%
  # NOTE: THIS FUNCTION CAN TAKE SEVERAL MINUTES TO RUN
  mt_plots_pathview(metab_id="KEGG_ids",
                    # n.pathways = 5,
                    # take results from statistical analysis called "Age met"
                    stat_name = "Age met",
                    # color scale function
                    color_scale = sign(statistic),
                    # metabolite filtering condition
                    metab_filter = p.adj < 0.05,
                    # get pathway list only from filtered metabolites
                    show_only_filtered = TRUE,
                    # kegg pathway files will be created in a folder called "Pathview_database" inside the current working directory
                    path_database = "./Pathview_database",
                    # output will be created in a folder called "Pathview_output" inside the current working directory
                    path_output = "./Pathview_output",
                    # set to false to speed-up (output files will be bigger in size)
                    same_layer = FALSE,
                    add_pwname_suffix = TRUE) %>%
                    {.}


# PART 7.2 - STATISTICAL RESULT PRESENTATION: MULTIPLE STATISTICS HEATMAP --------------------------------------------------

D1 <- D1 %>%
  # heading for html file
  mt_reporting_heading(heading = "Multiple Statistics Heatmap", lvl = 2) %>%
  # heatmap of all statistical results
  mt_plots_multstats_heatmap(cutoff = 0.05) %>%
  {.}


# PART 7.3 - STATISTICAL RESULT PRESENTATION: RESULT COMPARISON ----------------------------------------------------

D1 <- D1 %>%
  # heading for html file
  mt_reporting_heading(heading = "Result Comparison", lvl = 2) %>%
  # comparison plot
  mt_plots_stats_compare(stat1 = "Age met", filter1 = p.adj < 0.05,
                         D2 = D1, stat2 = "Diagnosis met", filter2 = p.adj < 0.05,
                         filter_op = "OR") %>%
                         {.}




# PART 8 - PARTIAL CORRELATION NETWORK ----------------------------------------------------

D1 <- D1 %>%
  # heading for html file
  mt_reporting_heading(heading = "Partial Correlation Network", lvl = 2) %>%
  # compute partial correlation matrix
  mt_stats_cormat_genenet(stat_name = "GGM") %>%
  # add multiple testing correction
  mt_post_multtest(stat_name = "GGM", method = "BH") %>%
  # plot network and color according to age analysis
  # NOTE: THIS FUNCTION CAN TAKE SEVERAL MINUTES TO RUN
  mt_plots_net(stat_name = "GGM", cor_filter = p.adj < 0.05, node_coloring = "Age met") %>%
  {.}


# PART 9 - PATHWAY AGGREGATION ANALYSIS ----------------------------------------------------
# NOTE ON BEST PRACTICES: This is now first aggregating the metabolite matrix into pathways creating a new matrix of pathway
#   concentration values, and then repeating the parts of the same pipeline as above. In a real scenario, you would not
#   perform both kinds of analysis within the same pipeline.

# create another SE object for second analysis branch (pathways)
D2 <- D

D2 <- D2 %>%
  # aggregate metabolites in the same pathways
  mt_modify_agg_pathways(pw_col = "pathway", method = "aggmean") %>%

  # STATISTICAL ANALYSIS, OUTCOME: AGE
  # heading for html file
  mt_reporting_heading(heading = "Pathway Aggregation Analysis", lvl = 1) %>%
  # heading for html file
  mt_reporting_heading(heading = "Age analysis", lvl = 2) %>%
  # Pearson correlation
  mt_stats_univ_cor(method = "pearson",
                    in_col = "Age",
                    stat_name = "Age pw") %>%
  # add multiple testing correction
  mt_post_multtest(stat_name = "Age pw", method = "BH") %>%
  # add stats logging
  mt_reporting_stats(stat_name = "Age pw", stat_filter = p.adj < 0.05) %>%
  # pvalue histogram
  mt_plots_pval_hist(stat_list = "Age pw") %>%
  # volcano plot as overview of results
  mt_plots_volcano(stat_name = "Age pw",
                   x = statistic,
                   feat_filter = p.adj < 0.05,
                   colour = p.adj < 0.05) %>%
  # scatter plot
  mt_plots_box_scatter(stat_name = "Age pw",
                       x = Age,
                       plot_type = "scatter",
                       feat_filter = p.adj < 1E-10,
                       feat_sort = p.value,
                       annotation = "{sprintf('P-value: %.2e', p.value)}\nP.adj: {sprintf('%.2e', p.adj)}") %>%

  # STATISTICAL ANALYSIS, OUTCOME: DIAGNOSIS
  # heading for html file
  mt_reporting_heading(heading = "Diagnosis analysis", lvl = 2) %>%
  # linear model for binary function (equivalent to t-test)
  mt_stats_univ_lm(formula = ~ Diagnosis,
                   stat_name = "Diagnosis pw") %>%
  # add fold change
  mt_post_fold_change(stat_name = "Diagnosis pw") %>%
  # add multiple testing correction
  mt_post_multtest(stat_name = "Diagnosis pw", method = "BH") %>%
  # add stats logging
  mt_reporting_stats(stat_name = "Diagnosis pw", stat_filter = p.adj < 0.05) %>%
  # pvalue histogram
  mt_plots_pval_hist(stat_list = "Diagnosis pw") %>%
  # volcano plot as overview of results
  mt_plots_volcano(stat_name = "Diagnosis pw",
                   x = fc,
                   feat_filter = p.adj < 0.05,
                   colour       = p.adj < 0.05) %>%
  # boxplot
  mt_plots_box_scatter(stat_name          ="Diagnosis pw",
                        x                  = Diagnosis,
                        fill               = Diagnosis,
                        plot_type          = "box",
                        feat_filter       = p.adj < 0.05,
                        feat_sort         = p.value,
                        annotation         = "{sprintf('P-value: %.2e', p.value)}\nPadj: {sprintf('%.2e', p.adj)}")

D2 <- D2 %>%
  # MULTIPLE STATISTICS HEATMAP
  # heading for html file
  mt_reporting_heading(heading = "Statistical Results Presentation", lvl = 2) %>%
  # heading for html file
  mt_reporting_heading(heading = "Multiple Statistics Heatmap", lvl = 2) %>%
  # heatmap of all statistical results
  mt_plots_multstats_heatmap(cutoff = 0.05) %>%

  # COMPARE STATISTICAL RESULTS
  # heading for html file
  mt_reporting_heading(heading = "Result Comparison", lvl = 2) %>%
  # comparison plot from pw analysis
  mt_plots_stats_compare(stat1 = "Age pw", filter1 = p.adj < 0.05,
                         D2 = D2, stat2 = "Diagnosis pw", filter2 = p.adj < 0.05,
                         filter_op = "OR")


# PART 10 - SUB PATHWAY ANALYSIS ----------------------------------------------------
# create another SE object for third analysis branch (sub-pathways)
D3 <- D

D3 <- D3 %>%
  # aggregate metabolites in the same pathways
  mt_modify_agg_pathways(pw_col = "SUB_PATHWAY", method = "aggmean") %>%

  # STATISTICAL ANALYSIS, OUTCOME: AGE
  # heading for html file
  mt_reporting_heading(heading = "Sub Pathway Aggregation Analysis", lvl = 1) %>%
  # heading for html file
  mt_reporting_heading(heading = "Age analysis", lvl = 2) %>%
  # Pearson correlation
  mt_stats_univ_cor(method = "pearson",
                    in_col = "Age",
                    stat_name = "Age sub pw") %>%
  # add multiple testing correction
  mt_post_multtest(stat_name = "Age sub pw", method = "BH") %>%
  # add stats logging
  mt_reporting_stats(stat_name = "Age sub pw", stat_filter = p.adj < 0.05) %>%
  # pvalue histogram
  mt_plots_pval_hist(stat_list = "Age sub pw") %>%
  # volcano plot as overview of results
  mt_plots_volcano(stat_name = "Age sub pw",
                   x = statistic,
                   feat_filter = p.adj < 0.05,
                   colour = p.adj < 0.05) %>%
  # scatter plot
  mt_plots_box_scatter(stat_name = "Age sub pw",
                       x = Age,
                       plot_type = "scatter",
                       feat_filter = p.adj < 1E-10,
                       feat_sort = p.value,
                       annotation = "{sprintf('P-value: %.2e', p.value)}\nP.adj: {sprintf('%.2e', p.adj)}") %>%

  # STATISTICAL ANALYSIS, OUTCOME: DIAGNOSIS
  # heading for html file
  mt_reporting_heading(heading = "Diagnosis analysis", lvl = 2) %>%
  # linear model for binary function (equivalent to t-test)
  mt_stats_univ_lm(formula = ~ Diagnosis,
                   stat_name = "Diagnosis sub pw") %>%
  # add fold change
  mt_post_fold_change(stat_name = "Diagnosis sub pw") %>%
  # add multiple testing correction
  mt_post_multtest(stat_name = "Diagnosis sub pw", method = "BH") %>%
  # add stats logging
  mt_reporting_stats(stat_name = "Diagnosis sub pw", stat_filter = p.adj < 0.05) %>%
  # pvalue histogram
  mt_plots_pval_hist(stat_list = "Diagnosis sub pw") %>%
  # volcano plot as overview of results
  mt_plots_volcano(stat_name = "Diagnosis sub pw",
                   x         = fc,
                   feat_filter = p.adj < 0.05,
                   colour      = p.adj < 0.05) %>%
  # boxplot
  mt_plots_box_scatter(stat_name          ="Diagnosis sub pw",
                       x                  = Diagnosis,
                       fill               = Diagnosis,
                       plot_type          = "box",
                       feat_filter        = p.adj < 0.05,
                       feat_sort          = p.value,
                       annotation         = "{sprintf('P-value: %.2e', p.value)}\nPadj: {sprintf('%.2e', p.adj)}")

D3 <- D3 %>%
  # MULTIPLE STATISTICS HEATMAP
  # heading for html file
  mt_reporting_heading(heading = "Statistical Results Presentation", lvl = 2) %>%
  # heading for html file
  mt_reporting_heading(heading = "Multiple Statistics Heatmap", lvl = 2) %>%
  # heatmap of all statistical results
  mt_plots_multstats_heatmap(cutoff = 0.05) %>%

  # COMPARE STATISTICAL RESULTS
  # heading for html file
  mt_reporting_heading(heading = "Result Comparison", lvl = 2) %>%
  # comparison plot from pw analysis
  mt_plots_stats_compare(stat1 = "Age sub pw", filter1 = p.adj < 0.05,
                         D2 = D3, stat2 = "Diagnosis sub pw", filter2 = p.adj < 0.05,
                         filter_op = "OR") %>%
  {.}




# PART 11 - SUB PATHWAY ANALYSIS  ----------------------------------------------------
# create another SE object for third analysis branch (super-pathways)
D4 <- D

D4 <- D4 %>%
  # aggregate metabolites in the same pathways
  mt_modify_agg_pathways(pw_col = "SUPER_PATHWAY", method = "aggmean") %>%

  # STATISTICAL ANALYSIS, OUTCOME: AGE
  # heading for html file
  mt_reporting_heading(heading = "Super Pathway Aggregation Analysis", lvl = 1) %>%
  # add tag
  mt_reporting_tag(tag_name = "Beginning of Super PW Age Analysis") %>%
  # heading for html file
  mt_reporting_heading(heading = "Age analysis", lvl = 2) %>%
  # Pearson correlation
  mt_stats_univ_cor(method = "pearson",
                    in_col = "Age",
                    stat_name = "Age super pw") %>%
  # add multiple testing correction
  mt_post_multtest(stat_name = "Age super pw", method = "BH") %>%
  # add stats logging
  mt_reporting_stats(stat_name = "Age super pw", stat_filter = p.adj < 0.05) %>%
  # pvalue histogram
  mt_plots_pval_hist(stat_list = "Age super pw") %>%
  # volcano plot as overview of results
  mt_plots_volcano(stat_name = "Age super pw",
                   x = statistic,
                   feat_filter = p.adj < 0.05,
                   colour = p.adj < 0.05) %>%
  # scatter plot
  mt_plots_box_scatter(stat_name = "Age super pw",
                       x = Age,
                       plot_type = "scatter",
                       feat_filter = p.adj < 1E-10,
                       feat_sort = p.value,
                       annotation = "{sprintf('P-value: %.2e', p.value)}\nP.adj: {sprintf('%.2e', p.adj)}") %>%

  # STATISTICAL ANALYSIS, OUTCOME: DIAGNOSIS
  # add tag
  mt_reporting_tag(tag_name = "Beginning of Super PW Diagnosis Analysis") %>%
  # heading for html file
  mt_reporting_heading(heading = "Diagnosis analysis", lvl = 2) %>%
  # linear model for binary function (equivalent to t-test)
  mt_stats_univ_lm(formula = ~ Diagnosis,
                   stat_name = "Diagnosis super pw") %>%
  # add fold change
  mt_post_fold_change(stat_name = "Diagnosis super pw") %>%
  # add multiple testing correction
  mt_post_multtest(stat_name = "Diagnosis super pw", method = "BH") %>%
  # add stats logging
  mt_reporting_stats(stat_name = "Diagnosis super pw", stat_filter = p.adj < 0.05) %>%
  # pvalue histogram
  mt_plots_pval_hist(stat_list = "Diagnosis super pw") %>%
  # volcano plot as overview of results
  mt_plots_volcano(stat_name = "Diagnosis super pw",
                   x =  fc,
                   feat_filter = p.adj < 0.05,
                   colour       = p.adj < 0.05) %>%
  # boxplot
  mt_plots_box_scatter(stat_name          ="Diagnosis super pw",
                       x                  = Diagnosis,
                       fill               = Diagnosis,
                       plot_type          = "box",
                       feat_filter       = p.adj < 0.05,
                       feat_sort         = p.value,
                       annotation         = "{sprintf('P-value: %.2e', p.value)}\nPadj: {sprintf('%.2e', p.adj)}")

D4 <- D4 %>%
  # MULTIPLE STATISTICS HEATMAP
  # heading for html file
  mt_reporting_heading(heading = "Statistical Results Presentation", lvl = 2) %>%
  # heading for html file
  mt_reporting_heading(heading = "Multiple Statistics Heatmap", lvl = 2) %>%
  # heatmap of all statistical results
  mt_plots_multstats_heatmap(cutoff = 0.05) %>%

  # COMPARE STATISTICAL RESULTS
  # heading for html file
  mt_reporting_heading(heading = "Result Comparison", lvl = 2) %>%
  # comparison plot from pw analysis
  mt_plots_stats_compare(stat1 = "Age super pw", filter1 = p.adj < 0.05,
                         D2 = D4, stat2 = "Diagnosis super pw", filter2 = p.adj < 0.05,
                         filter_op = "OR") %>%
  {.}

# PART 12 - COMPARE SUPER- AND SUB- PATHWAY ANALYSES  ----------------------------------------------------
D4 <- D4 %>% mt_plots_equalizer(stat1 = "Age super pw",
                                D2 = D3,
                                stat2 = "Age sub pw",
                                legend_fine="sub pathway",
                                legend_coarse='super pathway',
                                vline_fine = p.adj < 0.1,
                                vline_coarse = p.adj < 0.1) %>%
  # end timing
  mt_reporting_toc() %>%
  {.}

# PART 13 - CREATE ANALYSIS REPORTS ----------------------------------------------------

# OPTIONAL - this function will remove all plot entries in the pipeline
#   - mt_clean_remove_results(remove = "plots")

# metabolite analysis html report
D1 %>% mt_reporting_html(file = "Example_Pipeline_Metabolite_Analysis.html",
                         title = "Example Pipeline - Statistical Analysis")
# pathway analysis html report
D2 %>% mt_reporting_html(file = "Example_Pipeline_Pathway_Analysis.html",
                         title = "Example Pipeline - Pathway Aggregation Analysis")

# sub-pathway analysis html report
D3 %>% mt_reporting_html(file = "Example_Pipeline_Sub_Pathway_Analysis.html",
                         title = "Example Pipeline - Sub Pathway Aggregation Analysis")

# super-pathway analysis html report
D4 %>% mt_reporting_html(file = "Example_Pipeline_Super_Pathway_Analysis.html",
                         title = "Example Pipeline - Super Pathway Aggregation Analysis",
                         start_after = "Beginning of Super PW Age Analysis")

# create a combined report with both analsyses
mt_reporting_html_nonlinear(D_list = list(D1, D2, D3, D4), file = "ExamplePipeline_CombinedReport.html",
                            title = "Combined Report")
