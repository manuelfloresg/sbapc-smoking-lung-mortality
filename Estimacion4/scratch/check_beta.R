
res <- readRDS("d:/Dropbox/Investigacion/Bloomberg_2025/Estimacion4/results/batches/20260504_0044_FINAL_PROD/res_spec_linear_s11_quit.rds")
cat("\n--- ESTIMATOR COEFFICIENTS (resM) ---\n")
print(res$resM$inc_fit$beta_z)
cat("\n--- ESTIMATOR COEFFICIENTS (resF) ---\n")
print(res$resF$inc_fit$beta_z)
