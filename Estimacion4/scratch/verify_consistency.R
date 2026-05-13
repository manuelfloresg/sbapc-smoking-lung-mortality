
q_f <- readRDS('results/20260507_0153_STABLE_V14/raw_data/res_spec_linear_s4_freeze.rds')
q_q <- readRDS('results/20260507_0153_STABLE_V14/raw_data/res_spec_linear_s4_quit.rds')

cat('--- Combined Benchmark Check (Freeze vs Quit) ---\n')

# Check annual_bapc
b_f <- q_f$combined$annual_bapc
b_q <- q_q$combined$annual_bapc
diff_bapc <- mean(abs(b_f$deaths_hat - b_q$deaths_hat))
cat(sprintf('Mean Absolute Difference in annual_bapc: %.10f\n', diff_bapc))

# Check annual_anchor_noP
n_f <- q_f$combined$annual_anchor_noP
n_q <- q_q$combined$annual_anchor_noP
diff_noP <- mean(abs(n_f$deaths_hat - n_q$deaths_hat))
cat(sprintf('Mean Absolute Difference in annual_anchor_noP: %.10f\n', diff_noP))

# Check sex-specific objects for completeness
cat('\n--- Sex-Specific Benchmark Check (M) ---\n')
m_f <- q_f$resM$annual_bapc
m_q <- q_q$resM$annual_bapc
cat(sprintf('Mean Absolute Difference in resM$annual_bapc: %.10f\n', mean(abs(m_f$deaths_hat - m_q$deaths_hat))))

cat('\n--- Sex-Specific Benchmark Check (F) ---\n')
f_f <- q_f$resF$annual_bapc
f_q <- q_q$resF$annual_bapc
cat(sprintf('Mean Absolute Difference in resF$annual_bapc: %.10f\n', mean(abs(f_f$deaths_hat - f_q$deaths_hat))))
