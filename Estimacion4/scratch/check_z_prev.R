
q_f <- readRDS('results/20260507_0153_STABLE_V14/raw_data/res_spec_linear_s4_freeze.rds')
q_q <- readRDS('results/20260507_0153_STABLE_V14/raw_data/res_spec_linear_s4_quit.rds')

cat('--- z_prev (Informed log-offset) Check ---\n')
z_f <- q_f$resM$inc_fit$rates_all_full$z_prev
z_q <- q_q$resM$inc_fit$rates_all_full$z_prev

df_z <- data.frame(
  period = q_f$resM$inc_fit$rates_all_full$period,
  age = q_f$resM$inc_fit$rates_all_full$age,
  z_freeze = z_f,
  z_quit = z_q
)

# Look at 2050
cat('\nYear 2050 Sample (Age 50):\n')
print(df_z[df_z$period == 2050 & df_z$age == 50, ])

cat('\nMean difference in z_prev (Future only):\n')
fut_idx <- df_z$period > 2022
cat(sprintf('Mean(z_freeze - z_quit) in future: %.10f\n', mean(df_z$z_freeze[fut_idx] - df_z$z_quit[fut_idx])))

cat('\n--- Predicted Rates Check ---\n')
r_f <- q_f$resM$inc_fit$rates_all_full$rate_hat
r_q <- q_q$resM$inc_fit$rates_all_full$rate_hat
cat(sprintf('Mean Absolute Difference in rate_hat (Future): %.10f\n', mean(abs(r_f[fut_idx] - r_q[fut_idx]))))
