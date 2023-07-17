### -- Example: Average Duplicate Samples & features -- ###
# This script demonstrates the functionality of mt_modify_averagesample() and mt_modify_avg_features().

library(maplet)

# Average Dupcliate Samples -----------------
# Load toy SE object with duplicated samples
data("D_ex1")

# 100 samples and 100 features
D_ex1 %>% dim()
# 5 samples are duplicates according to Sample_ID
colData(D_ex1)$Sample_ID %>% duplicated() %>% sum()

# use maplet to average + combine samples
D_ex1 <- D_ex1 %>%
  mt_modify_avg_samples(group_col = "Sample_ID")

# 95 unique samples and 100 features
D_ex1 %>% assay() %>% as.data.frame() %>% dim()
colData(D_ex1)$Sample_ID %>% duplicated() %>% sum()


# Average Duplciate features -----------------
# Load toy SE object with duplicated features
data("D_ex2")

# 50 samples and 50 features
D_ex2 %>% dim()
# 3 features are duplicated according to COMP_ID
rowData(D_ex2)$COMP_ID %>% duplicated() %>% sum()

# use maplet to average + combine features
D_ex2 <- D_ex2 %>%
  mt_modify_avg_features(group_col = "COMP_ID")

# 47 unique features and 50 samples
D_ex2 %>% assay() %>% as.data.frame() %>% dim()
rowData(D_ex2)$COMP_ID %>% duplicated() %>% sum()

