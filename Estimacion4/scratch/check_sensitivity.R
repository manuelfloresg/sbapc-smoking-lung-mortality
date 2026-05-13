
library(dplyr)
# Cargar un resultado de la carpeta final
res_path <- "d:/Dropbox/Investigacion/Bloomberg_2025/Estimacion4/results/batches/20260504_0044_FINAL_PROD/res_misspec_tanh_s11_quit.rds"
res_both <- readRDS(res_path)
inc <- res_both$combined$annual_anchor

cat("\n--- COMPARISON OF Z VALUES (MALE, AGE 65) ---\n")
inc_sub <- inc %>% filter(sex == "M", age == 65) %>% arrange(period)

print(inc_sub %>% select(period, z_hat, zI_true_used) %>% head(10))
print(inc_sub %>% select(period, z_hat, zI_true_used) %>% tail(10))

cat("\n--- RANGE OF Z ---\n")
print(summary(inc_sub$z_hat))
print(summary(inc_sub$zI_true_used))

cat("\n--- SENSITIVITY CHECK (QUIT SCENARIO) ---\n")
# El drop en Z del estimador vs DGP
z_drop_hat <- inc_sub$z_hat[inc_sub$period == 2050] - inc_sub$z_hat[inc_sub$period == 2022]
z_drop_true <- inc_sub$zI_true_used[inc_sub$period == 2050] - inc_sub$zI_true_used[inc_sub$period == 2022]

cat("Z drop (Estimator): ", z_drop_hat, "\n")
cat("Z drop (DGP):       ", z_drop_true, "\n")
