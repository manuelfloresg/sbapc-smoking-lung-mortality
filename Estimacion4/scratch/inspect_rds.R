rb <- readRDS('results/20260506_1250_FINAL_DIAG/raw_data/res_spec_linear_s4_quit.rds')
print("Names of rb:")
print(names(rb))

print("Names of rb$resM:")
print(names(rb$resM))

if (!is.null(rb$resM$inc_fit)) {
  print("Names of rb$resM$inc_fit:")
  print(names(rb$resM$inc_fit))
} else {
  print("rb$resM$inc_fit is NULL")
}
