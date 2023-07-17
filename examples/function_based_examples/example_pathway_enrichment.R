### -- Example: Pathway Enrichment Analysis -- ###
# This script demonstrates the functionality of mt_stats_pathway_enrichment().


library(maplet)

file_data <- system.file("extdata", "example_data/sampledata.xlsx", package = "maplet")

# Preprocess dataset ------------------------------------------------------
D_pre <-
  mt_load_metabolon_v1(file=file_data, sheet="OrigScale") %>%
  mt_anno_pathways_hmdb(in_col = "HMDb_ID",
                        out_col = "kegg_db",
                        pwdb_name = "KEGG",
                        db_dir = system.file("extdata", "precalc/hmdb", package = "maplet")) %>%
  mt_anno_pathways_remove_redundant(feat_col = "HMDb_ID", pw_col = "kegg_db") %>%
  mt_pre_filter_missingness(feat_max=0.2) %>%
  mt_pre_filter_missingness(samp_max=0.1) %>%
  # batch correction by variable BATCH_MOCK
  mt_pre_batch_median(batch_col = "BATCH_MOCK") %>%
  # quotient normalization
  mt_pre_norm_quot() %>%
  # logging
  mt_pre_trans_log() %>%
  # KNN imputation
  mt_pre_impute_knn() %>%
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
  mt_post_multtest(stat_name = "comp", method = "BH")




# Apply mt_stats_pathway_enrichment ---------------------------------------

D_pw <-
  D_pre %>%
  mt_stats_pathway_enrichment(pw_col = "kegg_db",
                              stat_name = "comp",
                              cutoff = 0.4)


# Show results ---------------------------------------
metadata(D_pw)$pathways$enrichment_results
