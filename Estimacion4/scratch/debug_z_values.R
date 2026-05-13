library(dplyr)
rb <- readRDS('results/20260506_STABLE_V2_PURE/raw_data/res_spec_linear_s4_quit.rds')
diag_m <- rb$resM$diag

print("Names in diag_m:")
print(names(diag_m))

# Extract z_prev from future
z_f <- diag_m$z_prev_future
if (is.data.frame(z_f)) {
  summ_z <- z_f %>% group_by(period) %>% summarise(z = mean(z_prev, na.rm=TRUE), q = mean(q_eff, na.rm=TRUE))
  print("Summary of z_prev (Estimator):")
  print(summ_z %>% filter(period %in% c(2023, 2030, 2040, 2060)))
} else {
  print("z_prev_future is not a data frame, it is a:")
  print(class(z_f))
}
