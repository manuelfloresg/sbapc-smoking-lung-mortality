# audit_pilot_parity.R
source("runs/_source_all.R")
source("R/31_diagnostics_against_truth.R")

batch_id <- "20260504_1426_PILOT_PARITY"
batch_dir <- file.path(BAPC_PATHS$results, "batches", batch_id)

message(">>> Auditando batch: ", batch_id)

# 1. Agrupar resultados
res_files <- list.files(batch_dir, pattern = "\\.rds$", full.names = TRUE)
audit_list <- list()

for (f in res_files) {
  nm <- basename(f)
  parts <- strsplit(nm, "_")[[1]]
  seed <- as.integer(gsub("s", "", parts[4]))
  scn  <- gsub("\\.rds", "", parts[5])
  
  res <- readRDS(f)
  sim <- res$meta$sim_data
  
  # Verdad del DGP
  z_true_df <- sim$z_scen_true %>% dplyr::filter(period > 2022, sex == "M")
  z_true_mean <- mean(z_true_df$z_prev, na.rm = TRUE)
  
  # Estimación del Modelo
  # (Se encuentra en res$resM$inc_fit$rates_all_full o similar)
  res_sex <- res$resM
  rates_hat <- res_sex$inc_fit$rates_all_full
  z_hat_mean <- if (!is.null(rates_hat)) {
    mean(rates_hat$z_prev[rates_hat$period > 2022], na.rm = TRUE)
  } else NA_real_

  # Diagnóstico de incidencia (rel_error)
  diag <- compare_pipeline_to_truth(res, sim = sim, out_dir = NULL)
  bias_inc <- mean(diag$inc$rel_error[diag$inc$period > 2022], na.rm = TRUE)

  # Audit específico para edades avanzadas (para verificar que el 'quit' aplica a carry-states)
  rates_hat_old <- rates_hat %>% dplyr::filter(period > 2022, age >= 75, age <= 80)
  z_true_old <- res$meta$sim_data$z_scen_true %>% dplyr::filter(period > 2022, age >= 75, age <= 80)
  
  p_hat_old <- mean(rates_hat_old$p_cur, na.rm = TRUE)
  p_true_old <- mean(z_true_old$q_eff, na.rm = TRUE) # en DGP q_eff es el driver de incidencia
  
  message(sprintf("Scenario: %-10s | Bias Inc: %6.2f%% | z_true: %6.4f | z_hat: %6.4f | p_hat_75-80: %6.4f", 
                  scn, bias_inc * 100, z_true_mean, z_hat_mean, p_hat_old))

  audit_list[[nm]] <- tibble::tibble(
    seed = seed, scenario = scn,
    z_true = z_true_mean, z_hat = z_hat_mean,
    bias_inc = bias_inc,
    p_hat_old = p_hat_old, p_true_old = p_true_old
  )
}

audit_df <- dplyr::bind_rows(audit_list)

message(">>> RESUMEN DE SESGO POR ESCENARIO:")
print(audit_df %>% 
  dplyr::group_by(scenario) %>% 
  dplyr::summarise(
    z_true = mean(z_true, na.rm = TRUE),
    z_hat = mean(z_hat, na.rm = TRUE),
    bias_inc = mean(bias_inc, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ))

# Guardar auditoría completa
readr::write_csv(audit_df, file.path(batch_dir, "audit_parity_summary.csv"))
message(">>> Auditoría guardada en: ", file.path(batch_dir, "audit_parity_summary.csv"))
