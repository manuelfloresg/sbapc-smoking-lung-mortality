
library(dplyr)
# New results directory
DIR <- 'results/20260507_ESTIMATE_V14/raw_data/'
q_f <- readRDS(paste0(DIR, 'res_spec_linear_s4_freeze.rds'))
q_q <- readRDS(paste0(DIR, 'res_spec_linear_s4_quit.rds'))

cat('--- Sensitivity Check (ESTIMATE MODE) ---\n')
cat('Estimated Theta (M): ', q_f$resM$inc_fit$beta_P_eff, '\n')
cat('Estimated Theta (F): ', q_f$resF$inc_fit$beta_P_eff, '\n')

# Truth rates for M
tr_f <- q_f$inc_truth_grid %>% dplyr::filter(sex == "M")
tr_q <- q_q$inc_truth_grid %>% dplyr::filter(sex == "M")

# Informed rates for M
inf_f <- q_f$resM$inc_fit$rates_all_full
inf_q <- q_q$resM$inc_fit$rates_all_full

df_tr <- tr_f %>% dplyr::select(period, age, truth_f = rateI_base_true) %>%
         dplyr::inner_join(tr_q %>% dplyr::select(period, age, truth_q = rateI_scen_true), by=c("period","age"))

df_inf <- inf_f %>% dplyr::select(period, age, inf_f = rate_hat) %>%
          dplyr::inner_join(inf_q %>% dplyr::select(period, age, inf_q = rate_hat), by=c("period","age"))

df_comp <- df_tr %>% dplyr::inner_join(df_inf, by=c("period","age"))

cat('\nYear 2070 Sample (Age 50):\n')
row <- df_comp[df_comp$period == 2070 & df_comp$age == 50, ]
print(row)

cat('\nTruth Ratio:    ', row$truth_q / row$truth_f, '\n')
cat('Informed Ratio: ', row$inf_q / row$inf_f, '\n')

cat('\n--- Aggregated 2070 check ---\n')
tr_ratio <- sum(df_comp$truth_q[df_comp$period==2070]) / sum(df_comp$truth_f[df_comp$period==2070])
inf_ratio <- sum(df_comp$inf_q[df_comp$period==2070]) / sum(df_comp$inf_f[df_comp$period==2070])
cat('Truth 2070 (M) Ratio: ', tr_ratio, '\n')
cat('Informed 2070 (M) Ratio: ', inf_ratio, '\n')

# Final comparison of z_prev
cat('\n--- z_prev Check (2070) ---\n')
cat('Truth z_prev (Freeze): ', df_tr$truth_f[df_tr$period==2070 & df_tr$age==50], ' (wait, this is rate)\n')
cat('Informed z_prev (Freeze): ', q_f$resM$inc_fit$rates_all_full$z_prev[q_f$resM$inc_fit$rates_all_full$period==2070 & q_f$resM$inc_fit$rates_all_full$age==50], '\n')
cat('Informed z_prev (Quit):   ', q_q$resM$inc_fit$rates_all_full$z_prev[q_q$resM$inc_fit$rates_all_full$period==2070 & q_q$resM$inc_fit$rates_all_full$age==50], '\n')
