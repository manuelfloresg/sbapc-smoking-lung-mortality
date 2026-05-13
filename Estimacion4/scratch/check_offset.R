
q_f <- readRDS('results/20260507_0153_STABLE_V14/raw_data/res_spec_linear_s4_freeze.rds')
q_q <- readRDS('results/20260507_0153_STABLE_V14/raw_data/res_spec_linear_s4_quit.rds')

cat('--- offset_prev_rr Check ---\n')
o_f <- q_f$resM$inc_fit$rates_all_full$offset_prev_rr
o_q <- q_q$resM$inc_fit$rates_all_full$offset_prev_rr

df_o <- data.frame(
  period = q_f$resM$inc_fit$rates_all_full$period,
  age = q_f$resM$inc_fit$rates_all_full$age,
  o_freeze = o_f,
  o_quit = o_q
)

cat('\nYear 2050 Sample (Age 50):\n')
print(df_o[df_o$period == 2050 & df_o$age == 50, ])

cat('\nMean difference in offset_prev_rr (Future only):\n')
fut_idx <- df_o$period > 2022
cat(sprintf('Mean(o_freeze - o_quit) in future: %.10f\n', mean(df_o$o_freeze[fut_idx] - df_o$o_quit[fut_idx])))

cat('\n--- q_eff Check ---\n')
qeff_f <- q_f$resM$inc_fit$rates_all_full$q_eff
qeff_q <- q_q$resM$inc_fit$rates_all_full$q_eff
cat(sprintf('Mean(qeff_freeze - qeff_quit) in future: %.10f\n', mean(qeff_f[fut_idx] - qeff_q[fut_idx])))

cat('\n--- p_cur Check ---\n')
p_f <- q_f$resM$inc_fit$rates_all_full$p_cur
p_q <- q_q$resM$inc_fit$rates_all_full$p_cur
cat(sprintf('Mean(p_freeze - p_quit) in future: %.10f\n', mean(p_f[fut_idx] - p_q[fut_idx])))
