# runs/run_audit_simulation.R
# Canonical simulation audit for SBAPC production pipeline.
# Generates diagnostics for multiple seeds, DGPs, and scenarios.

source("runs/_runtime_setup.R")
source("runs/_source_all.R")
source("adapters/build_inputs_sim.R")

library(dplyr)
library(ggplot2)
library(readr)

# --- Configuration ---
seeds <- c(4, 11, 16)
dgps <- c("spec_linear", "misspec_tanh")
scenarios <- c("freeze", "up1pc", "down1pc", "quit")
cause_id <- "lung"
# ---------------------

results_list <- list()
deltas_list <- list()
support_list <- list()

for (seed in seeds) {
  for (dgp in dgps) {
    message("\n", paste0(rep("=", 60), collapse=""))
    message(">>> SEED: ", seed, " | DGP: ", dgp)
    message(paste0(rep("=", 60), collapse=""))
    
    # 1) Build common inputs for this Seed/DGP (using freeze as base for history)
    sim_base <- simulate_PIM_data(cause_id = cause_id, seed = seed, dgp = dgp, scenario_name = "freeze")
    inputs <- build_inputs_sim(sim_base, cause_id = cause_id)
    
    cfg_row <- tibble::tibble(
      cause_id = cause_id,
      AGE_M_MIN = sim_base$meta$age_min, AGE_M_MAX = sim_base$meta$age_max,
      AGE_P_MIN = sim_base$meta$age_min_p %||% sim_base$meta$age_min,
      AGE_P_MAX = sim_base$meta$age_max_p %||% sim_base$meta$age_max,
      AGE_I_MIN = sim_base$meta$age_min, AGE_I_MAX = sim_base$meta$age_max,
      L_I_MAX_YEARS = 3L,
      MORT_SHOCK_YEARS = list(integer(0)),
      DOWNWEIGHT_F = list(integer(0))
    )

    # Store freeze results to compute Deltas later
    freeze_mort <- NULL

    for (scen in scenarios) {
      message("\nSCENARIO: ", scen)
      
      # 2) Simulate scenario-specific truth
      sim_scen <- simulate_PIM_data(cause_id = cause_id, seed = seed, dgp = dgp, scenario_name = scen)
      
      # 3) Run Pipeline
      prev_cfg <- get_prev_config(scenario = scen)
      res_both <- run_pipeline_both_from_inputs(inputs = inputs, cfg_row = cfg_row, prev_cfg = prev_cfg)
      
      # 4) Generate Diagnostics
      out_dir <- file.path("results", "diagnostics", paste0("audit_", dgp), paste0("seed", seed, "_", scen))
      if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
      
      diag_res <- compare_pipeline_to_truth(res_both, sim_scen, out_dir = out_dir)
      
      # 5) Collect summary metrics
      bias_inc <- diag_res$inc %>%
        dplyr::filter(period <= 2022) %>%
        dplyr::group_by(sex) %>%
        dplyr::summarise(hist_bias = mean(rel_error) * 100, .groups = "drop")
      
      proj_bias <- diag_res$inc %>%
        dplyr::filter(period > 2022) %>%
        dplyr::group_by(sex) %>%
        dplyr::summarise(proj_bias = mean(rel_error) * 100, .groups = "drop")
      
      metrics <- bias_inc %>%
        dplyr::left_join(proj_bias, by = "sex") %>%
        dplyr::mutate(seed = seed, dgp = dgp, scenario = scen)
      
      results_list[[paste(seed, dgp, scen, sep="_")]] <- metrics

      # 6) Collect support metrics for consolidation
      if (!is.null(diag_res$support)) {
        support_list[[paste(seed, dgp, scen, sep="_")]] <- diag_res$support %>%
          dplyr::mutate(seed = seed, dgp = dgp, scenario = scen)
      }

      # 7) Delta calculation (Scenario vs Freeze)
      current_mort <- diag_res$mort %>% dplyr::select(period, sex, deaths_hat)
      if (scen == "freeze") {
        freeze_mort <- current_mort %>% dplyr::rename(deaths_freeze = deaths_hat)
      } else if (!is.null(freeze_mort)) {
        delta_df <- current_mort %>%
          dplyr::left_join(freeze_mort, by = c("period", "sex")) %>%
          dplyr::mutate(
            delta_deaths = deaths_hat - deaths_freeze,
            rel_delta = delta_deaths / pmax(deaths_freeze, 1e-12),
            seed = seed, dgp = dgp, scenario = scen
          )
        deltas_list[[paste(seed, dgp, scen, sep="_")]] <- delta_df
      }
      
      message("Completed: ", seed, "_", dgp, "_", scen)
    }
  }
}

# Consolidate and Save Results
final_metrics <- dplyr::bind_rows(results_list)
write.csv(final_metrics, "results/diagnostics/audit_full_summary.csv", row.names = FALSE)

final_deltas <- dplyr::bind_rows(deltas_list)
write.csv(final_deltas, "results/diagnostics/audit_deltas_summary.csv", row.names = FALSE)

final_support <- dplyr::bind_rows(support_list)
write.csv(final_support, "results/diagnostics/audit_support_summary.csv", row.names = FALSE)

message("\nAudit complete.")
message(" - Metrics saved to: results/diagnostics/audit_full_summary.csv")
message(" - Deltas saved to: results/diagnostics/audit_deltas_summary.csv")
message(" - Support saved to: results/diagnostics/audit_support_summary.csv")

print(as.data.frame(final_metrics))

# Optional Beep
if (requireNamespace("beepr", quietly = TRUE)) beepr::beep(2) else cat("\a")
