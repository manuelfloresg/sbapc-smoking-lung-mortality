library(dplyr)
library(ggplot2)

# 1. Cargar el resultado de la estimación (Informed SBAPC)
rds_scen <- "results/20260506_STABLE_V2_PURE/raw_data/res_spec_linear_s4_quit.rds"
rb <- readRDS(rds_scen)

# Offset del Estimador (z_prev futuro)
z_est <- rb$resM$diag$z_prev_future %>%
  group_by(period) %>%
  summarise(z_est = mean(z_prev, na.rm = TRUE), .groups = "drop") %>%
  mutate(source = "Estimator")

# 2. Cargar la 'Verdad' desde el objeto de simulación (si está disponible)
# O recalcularlo usando la lógica del DGP con los mismos parámetros
# En los RDS de simulación, solemos guardar el meta-data.

# Vamos a intentar sacar el z_true de los diagnósticos
# En compare_pipeline_to_truth comparamos rates, pero no guardamos el z del DGP directamente.
# Voy a extraerlo del objeto 'sim' si existe en el entorno o recrearlo.

# Para simplificar, compararé el 'q_eff' (proporción de abandono efectiva) 
# que es la base del offset.
q_est <- rb$resM$diag$z_prev_future %>%
  group_by(period) %>%
  summarise(q_eff_est = mean(q_eff, na.rm = TRUE), .groups = "drop")

# Ver los valores
print("--- Comparación de q_eff (Estimador) ---")
print(head(q_est))

# También voy a ver el valor del RR usado en la estimación
print(paste("RR_I usado en estimación (M):", rb$resM$inc_fit$rr_inc))
