
library(dplyr)
q_f <- readRDS('results/20260507_0153_STABLE_V14/raw_data/res_spec_linear_s4_freeze.rds')
q_q <- readRDS('results/20260507_0153_STABLE_V14/raw_data/res_spec_linear_s4_quit.rds')

get_total_cases <- function(obj) {
  df <- obj$resM$inc_fit$rates_all_full
  df %>% dplyr::group_by(period) %>% dplyr::summarise(total = sum(rate_hat * E, na.rm = TRUE), .groups = "drop")
}

t_f <- get_total_cases(q_f)
t_q <- get_total_cases(q_q)

cat('--- Total Predicted Incidence (M) ---\n')
cat('\nYear 2070:\n')
cat(sprintf('Freeze: %.2f\n', t_f$total[t_f$period == 2070]))
cat(sprintf('Quit:   %.2f\n', t_q$total[t_q$period == 2070]))
cat(sprintf('Ratio:  %.4f\n', t_q$total[t_q$period == 2070] / t_f$total[t_f$period == 2070]))
