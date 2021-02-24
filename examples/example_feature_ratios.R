### -- Example: Replace features with Ratios -- ###
# This script demonstrates the functionality of mt_modify_ratios() and mt_post_pgain().

library(maplet)

D <-
  # load data
  mt_load_metabolon_v1(file = system.file("extdata", "example_data/sampledata.xlsx", package = "maplet"), sheet = "OrigScale") %>%
  # timing start
  mt_reporting_tic() %>%

  ### Preprocessing ---------
  # heading
  mt_reporting_heading(heading = "Preprocessing") %>%
  mt_reporting_heading(heading = "Part 1", lvl=2) %>%
  # sample boxplot
  mt_plots_sample_boxplot() %>%
  # missingness plot
  mt_plots_missingness() %>%
  # filter features with >20% missing values, then samples with >10% missing values
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



  # Modify data frame to represent features as ratios ---------
  # GGM
  mt_stats_cormat_genenet(stat_name ="GGM") %>%
  mt_post_multtest(stat_name = "GGM", method="fdr") %>%
  mt_modify_ratios(stat_name = "GGM", edge_filter = p.adj<0.5, neighborhood = 2) %>% # liberal cutoff because this is a mock dataset

  # Perform analysis on feature ratios ---------
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
  mt_plots_pval_hist(stat_list = "comp") %>%
  # pgains
  mt_post_pgain(stat_name = "comp") %>%
  # Volcano plot as overview of results
  mt_plots_volcano(stat_name     = "comp",
                   feat_filter = p.adj < 0.1 & pgain > 2,
                   colour       = p.adj < 0.1 & pgain > 2) %>%

                   {.}


# Generate HTML report ---------

D %>% mt_reporting_html(file="example_ratios.html", output_calls = T)

