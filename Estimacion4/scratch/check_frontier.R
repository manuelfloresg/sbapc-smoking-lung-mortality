library(dplyr)
rb <- readRDS('results/20260506_1250_FINAL_DIAG/raw_data/res_spec_linear_s4_quit.rds')
diag_m <- rb$resM$diag

# z_hist y z_future son dataframes (deberían serlo por mi último cambio)
z_h <- diag_m$z_hist
z_f <- diag_m$z_future

print("--- OFFSET FRONTIER CHECK (Mean across ages) ---")
h22 <- z_h %>% filter(period == 2022) %>% summarise(z = mean(z_prev, na.rm=TRUE), q = mean(q_eff, na.rm=TRUE))
f23 <- z_f %>% filter(period == 2023) %>% summarise(z = mean(z_prev, na.rm=TRUE), q = mean(q_eff, na.rm=TRUE))

print(paste("2022 (Hist): z =", h22$z, "| q =", h22$q))
print(paste("2023 (Fut) : z =", f23$z, "| q =", f23$q))

# Ver si hay un salto en q_eff
print(paste("Jump in q_eff:", f23$q - h22$q))
print(paste("Jump in z_prev:", f23$z - h22$z))
