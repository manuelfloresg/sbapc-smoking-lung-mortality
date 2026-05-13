library(tidyverse)
r_f <- readRDS("d:/Dropbox/Investigacion/Bloomberg_2025/Estimacion4/results/20260505_STABLE_V1/raw_data/res_spec_linear_s4_freeze.rds")
r_q <- readRDS("d:/Dropbox/Investigacion/Bloomberg_2025/Estimacion4/results/20260505_STABLE_V1/raw_data/res_spec_linear_s4_quit.rds")

df_f <- r_f$resM$diag$z_prev_future %>% filter(period >= 2020, period <= 2030, age == 60) %>% select(period, q_eff, offset_prev_rr, coef_fc_offset_I)
df_q <- r_q$resM$diag$z_prev_future %>% filter(period >= 2020, period <= 2030, age == 60) %>% select(period, q_eff, offset_prev_rr, coef_fc_offset_I)

print("FREEZE (Age 60):")
print(df_f)
print("QUIT (Age 60):")
print(df_q)

print("Difference in Total Offset (Quit - Freeze):")
print(df_q$coef_fc_offset_I - df_f$coef_fc_offset_I)
