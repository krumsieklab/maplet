# ------------ maplet Meta-Functions Example Pipeline ------------
#  One of the main features of maplet is the ability to easily group reusable function steps into meta-functions.
#  In this example, we demonstrate how to create these meta-functions and how these can be used in
#  a workflow. This example uses meta-function to produce a workflow identical to that demonstrated for
#  pathway analyses in the example_pipeline.R file, lines 347-614.

# NOTE: This example is divided into sections. Use SHIFT+CTRL+O to see the document outline.


library(maplet)
library(tidyverse)

# Users may notice that the last 3 pathway analyses in the example pipeline (lines 347-614) are identical. 
# When this is the case, it is convenient to group maplet functions into groups of meta functions, 
# like the examples age_analysis, diagnosis_analysis, and plot_stats below. When these meta-functions are
# passed D as the first argument and return D, they can be seamlessly chained to maplet functions.

# META-FUNCTION DEFINITIONS ----------------------------------------------------

age_analysis <- function(D, stat_name){
 
  D %<>%
  # STATISTICAL ANALYSIS, OUTCOME: AGE
  # heading for html file
  mt_reporting_heading(heading = "Pathway Aggregation Analysis", lvl = 1) %>%
    # heading for html file
    mt_reporting_heading(heading = "Age analysis", lvl = 2) %>%
    # Pearson correlation
    mt_stats_univ_cor(method = "pearson",
                      in_col = "Age",
                      stat_name = stat_name) %>%
    # add multiple testing correction
    mt_post_multtest(stat_name = stat_name, method = "BH") %>%
    # add stats logging
    mt_reporting_stats(stat_name = stat_name, stat_filter = p.adj < 0.05) %>%
    # pvalue histogram
    mt_plots_pval_hist(stat_list = stat_name) %>%
    # volcano plot as overview of results
    mt_plots_volcano(stat_name = stat_name,
                     x = statistic,
                     feat_filter = p.adj < 0.05,
                     colour = p.adj < 0.05) %>%
    # scatter plot
    mt_plots_box_scatter(stat_name = stat_name,
                         x = Age,
                         plot_type = "scatter",
                         feat_filter = p.adj < 1E-10,
                         feat_sort = p.value,
                         annotation = "{sprintf('P-value: %.2e', p.value)}\nP.adj: {sprintf('%.2e', p.adj)}") 
  
  D
}

diagnosis_analysis <- function(D, stat_name){
  
  D %<>% 
  # STATISTICAL ANALYSIS, OUTCOME: DIAGNOSIS
  # heading for html file
  mt_reporting_heading(heading = "Diagnosis analysis", lvl = 2) %>%
    # linear model for binary function (equivalent to t-test)
    mt_stats_univ_lm(formula = ~ Diagnosis,
                     stat_name = stat_name) %>%
    # add fold change
    mt_post_fold_change(stat_name = stat_name) %>%
    # add multiple testing correction
    mt_post_multtest(stat_name = stat_name, method = "BH") %>%
    # add stats logging
    mt_reporting_stats(stat_name = stat_name, stat_filter = p.adj < 0.05) %>%
    # pvalue histogram
    mt_plots_pval_hist(stat_list = stat_name) %>%
    # volcano plot as overview of results
    mt_plots_volcano(stat_name = stat_name,
                     x = fc,
                     feat_filter = p.adj < 0.05,
                     colour       = p.adj < 0.05) %>%
    # boxplot
    mt_plots_box_scatter(stat_name          =stat_name,
                         x                  = Diagnosis,
                         fill               = Diagnosis,
                         plot_type          = "box",
                         feat_filter       = p.adj < 0.05,
                         feat_sort         = p.value,
                         annotation         = "{sprintf('P-value: %.2e', p.value)}\nPadj: {sprintf('%.2e', p.adj)}")
  
  D
  
}

plot_stats <- function(D, stat1, stat2){
  
  D %<>%
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
    mt_plots_stats_compare(stat1 = stat1, filter1 = p.adj < 0.05,
                           D2 = D, stat2 = stat2, filter2 = p.adj < 0.05,
                           filter_op = "OR")
  
  D
  
}

# MAIN WORKFLOW ----------------------------------------------------
# GET D BY RUNNING LINES 22-171 FROM example_pipeline.R

# create copies of SE objects for pathway analyses
D2 <- D
D3 <- D
D4 <- D

# Pathway analysis
D2 <- D2 %>%
  # aggregate metabolites in the same pathways
  mt_modify_agg_pathways(pw_col = "pathway", method = "aggmean") %>%
  # run age analysis
  age_analysis(stat_name = "Age pw") %>%
  # run diagnosis analysis
  diagnosis_analysis(stat_name = "Diagnosis pw") %>%
  # plot stats results
  plot_stats(stat1 = "Age pw", stat2 = "Diagnosis pw") %>%
  # pathway analysis html report
  mt_reporting_html(file = "Example_Pipeline_Pathway_Analysis.html",
                    title = "Example Pipeline - Pathway Aggregation Analysis")


# Sub-pathway analysis
D3 <- D3 %>%
  # aggregate metabolites in the same pathways
  mt_modify_agg_pathways(pw_col = "SUB_PATHWAY", method = "aggmean") %>%
  # run age analysis
  age_analysis(stat_name = "Age sub pw") %>%
  # run diagnosis analysis
  diagnosis_analysis(stat_name = "Diagnosis sub pw") %>%
  # plot stats results
  plot_stats(stat1 = "Age sub pw", stat2 = "Diagnosis sub pw") %>%
  # sub=pathway analysis html report
  mt_reporting_html(file = "Example_Pipeline_Sub_Pathway_Analysis.html",
                    title = "Example Pipeline - Sub Pathway Aggregation Analysis")


# Super-pathway analysis
D4 <- D4 %>%
  # aggregate metabolites in the same pathways
  mt_modify_agg_pathways(pw_col = "SUPER_PATHWAY", method = "aggmean") %>%
  # run age analysis
  age_analysis(stat_name = "Age super pw") %>%
  # run diagnosis analysis
  diagnosis_analysis(stat_name = "Diagnosis super pw") %>%
  # plot stats results
  plot_stats(stat1 = "Age super pw", stat2 = "Diagnosis super pw") %>%
  # super-pathway analysis html report
  mt_reporting_html(file = "Example_Pipeline_Super_Pathway_Analysis.html",
                    title = "Example Pipeline - Super Pathway Aggregation Analysis",
                    start_after = "Beginning of Super PW Age Analysis")
  
