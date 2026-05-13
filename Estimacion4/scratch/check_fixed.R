
res <- readRDS("d:/Dropbox/Investigacion/Bloomberg_2025/Estimacion4/results/batches/20260504_0044_FINAL_PROD/res_spec_linear_s11_quit.rds")
cat("\n--- FIXED EFFECTS (MALE) ---\n")
if (!is.null(res$resM$inc_fit$summary_fixed)) {
  print(res$resM$inc_fit$summary_fixed)
} else {
  cat("summary_fixed is NULL\n")
}
