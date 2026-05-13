
q_f <- readRDS('results/20260507_ESTIMATE_V14/raw_data/res_spec_linear_s4_freeze.rds')
cat('beta_P_eff (M): ', q_f$resM$inc_fit$beta_P_eff, '\n')
cat('beta_P_rule (M): ', q_f$resM$inc_fit$beta_P_rule, '\n')
cat('beta_P (M):      ', q_f$resM$inc_fit$beta_P, '\n')
