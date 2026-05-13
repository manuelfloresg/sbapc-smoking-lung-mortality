
res <- readRDS("d:/Dropbox/Investigacion/Bloomberg_2025/Estimacion4/results/batches/20260504_0044_FINAL_PROD/res_spec_linear_s11_quit.rds")
fit <- res$resM$inc_fit %||% res$resM$inc_fit_bapc
if (!is.null(fit$summary_fixed)) {
  print(fit$summary_fixed)
} else {
  cat("summary_fixed NOT FOUND in resM$inc_fit or inc_fit_bapc\n")
}
