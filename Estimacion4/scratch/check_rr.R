library(tidyverse)
res <- readRDS("d:/Dropbox/Investigacion/Bloomberg_2025/Estimacion4/results/20260505_STABLE_V1/raw_data/res_spec_linear_s4_freeze.rds")
print("resM meta:")
print(res$resM$params)
print("resF meta:")
print(res$resF$params)

# Check rr_inc used in the result
print("RR_I from resM:")
print(res$resM$inc_fit$rr_inc) # Wait, where is it stored?
