
q_f <- readRDS('results/20260507_0153_STABLE_V14/raw_data/res_spec_linear_s4_freeze.rds')
q_q <- readRDS('results/20260507_0153_STABLE_V14/raw_data/res_spec_linear_s4_quit.rds')

cat('--- rate_hat Check ---\n')
df_r <- data.frame(
  period = q_f$resM$inc_fit$rates_all_full$period,
  age = q_f$resM$inc_fit$rates_all_full$age,
  r_freeze = q_f$resM$inc_fit$rates_all_full$rate_hat,
  r_quit = q_q$resM$inc_fit$rates_all_full$rate_hat,
  o_freeze = q_f$resM$inc_fit$rates_all_full$offset_prev_rr,
  o_quit = q_q$resM$inc_fit$rates_all_full$offset_prev_rr
)

cat('\nYear 2050 Sample (Age 50):\n')
print(df_r[df_r$period == 2050 & df_r$age == 50, ])

cat('\nImplied ratio (r_quit / r_freeze):\n')
row <- df_r[df_r$period == 2050 & df_r$age == 50, ]
cat(sprintf('Actual Ratio: %.10f\n', row$r_quit / row$r_freeze))
cat(sprintf('Expected Ratio (exp(o_quit - o_freeze)): %.10f\n', exp(row$o_quit - row$o_freeze)))
