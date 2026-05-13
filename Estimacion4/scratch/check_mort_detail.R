
library(dplyr)
q_f <- readRDS('results/20260507_ESTIMATE_V14/raw_data/res_spec_linear_s4_freeze.rds')
q_q <- readRDS('results/20260507_ESTIMATE_V14/raw_data/res_spec_linear_s4_quit.rds')

cat('--- Mortality Projection Detail Check ---\n')
m_f <- q_f$resM$mort_anchor_pred_detail %>% filter(period == 2070, cohort_true == 2020)
m_q <- q_q$resM$mort_anchor_pred_detail %>% filter(period == 2070, cohort_true == 2020)

cat('Freeze mu_hat (2070, Cohort 2020): ', m_f$mu_hat, '\n')
cat('Quit mu_hat (2070, Cohort 2020):   ', m_q$mu_hat, '\n')

cat('\nNames in pred_detail:\n')
print(names(m_q))

cat('\nSum of mu_hat (2070):\n')
cat('Freeze: ', sum(q_f$resM$mort_anchor_pred_detail$mu_hat[q_f$resM$mort_anchor_pred_detail$period == 2070], na.rm=T), '\n')
cat('Quit:   ', sum(q_q$resM$mort_anchor_pred_detail$mu_hat[q_q$resM$mort_anchor_pred_detail$period == 2070], na.rm=T), '\n')
