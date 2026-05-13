
res <- readRDS("d:/Dropbox/Investigacion/Bloomberg_2025/Estimacion4/results/batches/20260504_0044_FINAL_PROD/res_spec_linear_s11_quit.rds")
if (!is.null(res$params)) {
  cat("\nparams names:\n")
  print(names(res$params))
}
