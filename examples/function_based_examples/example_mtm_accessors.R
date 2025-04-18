### -- Example: Access plots and stats tables -- ###
# This script demonstrates how to use the mtm accessor functions to access plots and statistical tables.

library(maplet)

D <-
  # load data
  mt_load_metabolon_v1(file = system.file("extdata", "example_data/sampledata.xlsx", package = "maplet"), sheet = "OrigScale") %>%
  # set global settings to save all versions of assays
  mt_settings(settings = c("save_all_assays"=T)) %>%
  # timing start
  mt_reporting_tic() %>%
  
  ###
  # heading
  mt_reporting_heading(heading = "Preprocessing") %>%
  mt_reporting_heading(heading = "Part 1", lvl=2) %>%
  # sample boxplot
  mt_plots_sample_boxplot() %>%
  # missingness plot
  mt_plots_missingness() %>%
  # filter metabolites with >20% missing values, then samples with >10% missing values
  mt_pre_filter_missingness(feat_max=0.2) %>%
  mt_pre_filter_missingness(samp_max=0.1) %>%
  # batch correction by variable BATCH_MOCK
  mt_pre_batch_median(batch_col = "BATCH_MOCK") %>%
  # heading
  mt_reporting_heading(heading = "Part 2", lvl=2) %>%
  # quotient normalization
  mt_pre_norm_quot() %>%
  # check if there is any correlation between normalization factors and outcomes (bad sign if so)
  mt_plots_dilution_factor(in_col="num1") %>%
  mt_plots_dilution_factor(in_col="Group") %>%
  # logging
  mt_pre_trans_log() %>%
  # KNN imputation
  mt_pre_impute_knn() %>%
  # outlier detection (multivariate) & visualization
  # mt_pre_outlier(method="mahalanobis", pval=0.01, reduce.dim = T) %>%
  mt_plots_pca(color='outlier_mahalanobis') %>%
  # final sample boxplot
  mt_plots_sample_boxplot(color=Group, title = 'final') %>%
  # PCA, colored by some rowData() fields... this function shows 2 PCs
  mt_plots_pca(color=Group, shape=BATCH_MOCK, size=NUM_MOCK) %>%
  # heatmap
  mt_plots_heatmap(scale_data = T) %>%
  
  ###
  mt_reporting_heading(heading = "Statistics") %>%
  # linear model, differential test on Group
  mt_stats_univ_lm(
    formula      = ~ Group,
    samp_filter = (Group %in% c("treatment1","treatment2")),
    stat_name         = "comp",
    n_cores     = 1
  ) %>%
  # add fold changes to result tables
  mt_post_fold_change(stat_name = "comp") %>%
  # add multiple testing correction
  mt_post_multtest(stat_name = "comp", method = "BH") %>%
  # p-value histogram
  mt_plots_pval_hist() %>%
  # Volcano plot as overview of results
  mt_plots_volcano(stat_name     = "comp",
                   feat_filter = p.adj < 0.1,
                   colour       = p.value < 0.05) %>%
  
  ###
  # heading
  mt_reporting_heading("All boxplots") %>%
  # boxplots
  mt_plots_box_scatter(stat_name           = "comp",
                       plot_type = "box",
                       x                  = Group,
                       fill               = Group,
                       correct_confounder = ~BATCH_MOCK,
                       feat_filter       = p.value<0.01,
                       feat_sort         = p.value,
                       annotation         = "{sprintf('P-value: %.1e', p.value)}\nStatistic: {sprintf('%.2f', statistic)}") %>%
  # final timing
  mt_reporting_toc() %>%
  
  {.}

# Access Results ----------------------------

# the function mtm_res_get_entries can be used to return a list of results from any given namespace For example, to retrieve the 
# results lists from the batch correction steps pass the list of terms c("pre", "batch"). The terms must be provided in the order 
# they appear in the function name.
pre_list <- D %>% mtm_res_get_entries(c("pre"))   # get list of results from all functions beginning with "mt_pre"
pre_filter_list <- D %>% mtm_res_get_entries(c("pre", "filter"))  # get results from all functions beginning with "mt_pre_filter"

# the function mtm_get_stat_by_name returns a stats table given the "stat_name" value of a mt_stats function call
lm_res_table <- D %>% mtm_get_stat_by_name("comp")  # return stats table for statistical analysis "comp"

# the function mtm_res_get_plots
all_plots_list <- D %>% mtm_res_get_plots()  # return list of all plots in order they were generated

# the function mtm_plot_all_tofile takes a list of plots and 
mtm_plot_all_tofile(all_plots_list[1:4], file="example_plots.pdf")  # output a list of plots to a file

# Access Assay Versions ----------------------------
# When the global setting save_all_assays is TRUE, maplet saves a version of the assay every time it is changed. Each assay version
# has a "pointer" and each function records which assay version it is associated with. Users have two available methods for
# accessing assay versions - directly by assay pointer or indirectly by function uuid.

uuid_df <- D %>% mtm_res_get_uuids() # get uuids for all function calls
fun_uuid <- uuid_df[uuid_df$call_order==12,"UUIDs"] # get uuid for function call mt_pre_norm_quot (function call #12)
uuid_sel_assay <- D %>% mtm_get_assay_by_id(uuid = fun_uuid) # 
  
assay_ptr_df <- D %>% mtm_res_get_ptrs() # get assay pointers for all function calls
assay_ptr <- assay_ptr_df[assay_ptr_df$call_order==16,"assay_ptrs"] # get assay pointer for function call mt_pre_impute_knn (function call #16)
ptr_sel_assay <- D %>% mtm_get_assay_by_ptr(ptr = assay_ptr)

# Access UUIDs for Reporting ----------------------------
# The function mtm_res_get_uuids() is also useful for selecting a point to start an HTML report. The reporting function
# mt_reporting_html accepts an argument "start_after" which expects either a tag name or a UUID.

# start HTML report after the pipeline call mt_reporting_heading(heading="Statistics") (pipeline step #21)
start_uuid <- uuid_df[uuid_df$call_order==21,"UUIDs"]
D %>% mt_reporting_html(file="Start_After_UUID.html", start_after = start_uuid)
