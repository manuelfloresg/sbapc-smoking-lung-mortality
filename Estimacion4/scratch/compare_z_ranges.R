
res <- readRDS("d:/Dropbox/Investigacion/Bloomberg_2025/Estimacion4/results/batches/20260504_0044_FINAL_PROD/res_spec_linear_s11_quit.rds")
# El z del estimador está en rates_all_full$z
inc_hat <- res$resM$inc_fit$rates_all_full
# El z del DGP (una semilla)
library(dplyr)
# Cargar la verdad (necesitamos simularla o buscarla en los audit_plots)
# Pero podemos verla en el objeto res si lo guardamos (no lo guardamos completo)
# Usaremos el script simulate_PIM_data directamente para ver qué genera con esa semilla
source("runs/_source_all.R")
sim_true <- simulate_PIM_data(seed = 11, dgp = "spec_linear", scenario_name = "quit")

cat("\n--- ESTIMATOR Z RANGE ---\n")
print(summary(inc_hat$z))

cat("\n--- DGP Z RANGE ---\n")
print(summary(sim_true$inc_truth_grid$zI_true_used))
