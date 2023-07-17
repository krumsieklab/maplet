### -- Example: Reorder Samples Annotations for Plotting -- ###
# This script demonstrates the functionality of mt_anno_reorder_factor().

library(maplet)

file_data <- system.file("extdata", "example_data/sampledata.xlsx", package = "maplet")

D <-
  # load data
  mt_load_metabolon_v1(file = file_data, sheet = "OrigScale") %>%
  # timing start
  mt_reporting_tic() %>%

  # Preprocessing ---------
  # heading
  mt_reporting_heading(heading = "Sample Annotation Re-ordering") %>%
  mt_reporting_text(text = "Users can reorder specific sample annotations using the function mt_anno_reorder_factor. Annotations remain
                    in this new order for the remainder of the pipeline.") %>%
  # sample boxplot
  mt_plots_sample_boxplot(color=Group) %>%
  # filter metabolites with >20% missing values, then samples with >10% missing values
  mt_pre_filter_missingness(feat_max=0.2) %>%
  mt_pre_filter_missingness(samp_max=0.1) %>%
  # batch correction by variable BATCH_MOCK
  mt_pre_batch_median(batch_col = "BATCH_MOCK") %>%
  # heading
  mt_reporting_heading(heading = "Part 2", lvl=2) %>%
  # quotient normalization
  mt_pre_norm_quot() %>%
  mt_plots_dilution_factor(in_col="Group") %>%

  # Reorder groups ---------
  mt_anno_reorder_factor(col_name = "Group", new_order = c("Vehicle", "treatment1", "treatment2")) %>%
  # dilution plot after reordering
  mt_plots_dilution_factor(in_col = "Group") %>%
  # final sample boxplot (with reodering)
  mt_plots_sample_boxplot(color=Group, title = 'final')


# Generate HTML Report ---------
D %>% mt_reporting_html(file="example_factor_reorder.html", output_calls = T)
