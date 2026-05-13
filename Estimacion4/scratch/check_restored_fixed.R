
library(dplyr)
DIR <- 'results/20260507_ESTIMATE_V14/raw_data/'
q_f <- readRDS(paste0(DIR, 'res_spec_linear_s4_freeze.rds'))
q_q <- readRDS(paste0(DIR, 'res_spec_linear_s4_quit.rds'))

cat('--- Sensitivity Check (FIXED MODE - RESTORED) ---\n')
cat('Beta Mode: ', q_f$meta$beta_mode, '\n')

# Incidence check
inf_f_inc <- q_f$resM$inc_fit$rates_all %>% filter(period == 2070, age == 50)
inf_q_inc <- q_q$resM$inc_fit$rates_all %>% filter(period == 2070, age == 50)
cat('Incidence Rate (Freeze, 2070, Age 50): ', inf_f_inc$rate_hat, '\n')
cat('Incidence Rate (Quit, 2070, Age 50):   ', inf_q_inc$rate_hat, '\n')
cat('Incidence Ratio (Quit/Freeze):        ', inf_q_inc$rate_hat / inf_f_inc$rate_hat, '\n')

# Mortality check
inf_f_mort <- q_f$resM$annual_anchor %>% filter(period == 2070)
inf_q_mort <- q_q$resM$annual_anchor %>% filter(period == 2070)
cat('\nMortality Deaths (Freeze, 2070): ', inf_f_mort$deaths_hat, '\n')
cat('Mortality Deaths (Quit, 2070):   ', inf_q_mort$deaths_hat, '\n')
cat('Mortality Ratio (Quit/Freeze):   ', inf_q_mort$deaths_hat / inf_f_mort$deaths_hat, '\n')
