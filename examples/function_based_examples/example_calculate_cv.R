### -- Example: Calculate Coefficient of Variation (CV) -- ###
# This script demonstrates the functionality of mt_modify_cv().

library(maplet)
data("D_ex1")

# D_ex1 has two types of QC pools
colData(D_ex1)$QCID %>% table()

# Add cv from both QC pools to rowData
D_ex1 <- D_ex1 %>%
  mt_pre_cv(qc_samples = QCID == "QC_1", out_col = "QC1_cv") %>%
  mt_pre_cv(qc_samples = QCID == "QC_2", out_col = "QC2_cv")

rowData(D_ex1) %>% as.data.frame() %>% View()
