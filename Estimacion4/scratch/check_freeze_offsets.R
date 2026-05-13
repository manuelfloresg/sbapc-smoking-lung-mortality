library(tidyverse)
r_f <- readRDS("d:/Dropbox/Investigacion/Bloomberg_2025/Estimacion4/results/20260505_STABLE_V1/raw_data/res_spec_linear_s4_freeze.rds")
df_f <- r_f$resM$diag$z_prev_future %>% filter(period == 2023, age == 60) %>% select(period, q_eff, offset_prev_rr, coef_fc_offset_I, coef_fc_offset_I_epi, coef_fc_offset_I_apc)
print(df_f)
