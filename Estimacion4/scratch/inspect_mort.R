
library(dplyr)
q_f <- readRDS('results/20260507_ESTIMATE_V14/raw_data/res_spec_linear_s4_freeze.rds')
cat('Names in mort_anchor_pred_detail (pred_base):\n')
print(names(q_f$resM$mort_anchor_pred_detail))

# We can't see mort_data_cond directly from RDS, but we can see mort_anchor_data_cond
cat('\nNames in mort_anchor_data_cond:\n')
print(names(q_f$resM$mort_anchor_data_cond))

# Check if age can be reconstructed
cat('\nIs age present in pred_detail? ', 'age' %in% names(q_f$resM$mort_anchor_pred_detail), '\n')
cat('Is cohort present in pred_detail? ', 'cohort' %in% names(q_f$resM$mort_anchor_pred_detail), '\n')
