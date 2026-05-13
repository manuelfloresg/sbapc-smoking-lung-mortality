
batch_id <- "20260504_0044_FINAL_PROD"
seed <- 1
dgp <- "spec_linear"
scen <- "quit"

path <- file.path("d:/Dropbox/Investigacion/Bloomberg_2025/Estimacion4/results/batches", batch_id, 
                  paste0("res_", dgp, "_s", seed, "_", scen, ".rds"))

if (file.exists(path)) {
  obj <- readRDS(path)
  if (!is.null(obj$resM$diag)) {
    cat("Names of obj$resM$diag:", paste(names(obj$resM$diag), collapse=", "), "\n")
  } else {
    cat("obj$resM$diag is NULL\n")
  }
}
