library(dplyr)
rb <- readRDS('results/20260506_1250_FINAL_DIAG/raw_data/res_spec_linear_s4_quit.rds')

# Los datos están en inc_fit$z_hist e inc_fit$z_future
z_h <- rb$resM$inc_fit$z_hist
z_f <- rb$resM$inc_fit$z_future

print("--- FRONTIER OFFSET DIAGNOSTIC ---")
h22 <- z_h %>% filter(period == 2022) %>% summarise(z = mean(z_prev, na.rm=TRUE), q = mean(q_eff, na.rm=TRUE))
f23 <- z_f %>% filter(period == 2023) %>% summarise(z = mean(z_prev, na.rm=TRUE), q = mean(q_eff, na.rm=TRUE))

print(paste("2022 (History): q_eff =", round(h22$q, 4), "| z_prev =", round(h22$z, 4)))
print(paste("2023 (Future) : q_eff =", round(f23$q, 4), "| z_prev =", round(f23$z, 4)))

print(paste("DELTA z_prev:", round(f23$z - h22$z, 4)))
