# Autopsia de la Semilla 1 con TRACEBACK
options(error = function() {
  cat("\n--- STACK TRACE ---\n")
  traceback(20)
  q("no", status = 1)
})

source("runs/_source_all.R")

seed <- 1
dgp  <- "spec_linear"
scn  <- "freeze"

message(">>> Re-simulando datos...")
sim <- simulate_PIM_data(cause_id = dgp, seed = seed, dgp = dgp, scenario_name = scn)
inputs <- build_inputs_sim(sim, cause_id = dgp)
prev_cfg <- make_prev_config(scenario = scn)

message(">>> Lanzando pipeline (TRACEBACK MODE)...")
res <- run_pipeline_both(
  mort_hist_tbl = inputs$mort_hist, pop_all_tbl = inputs$pop_all,
  inc_hist_tbl = inputs$inc_hist, prev_micro_df = inputs$prev_data,
  cause_id_override = dgp, beta_mode = "fixed_rr_offset", prev_cfg = prev_cfg,
  emit_prev_diag_write = FALSE
)

message(">>> FIN DE AUTOPSIA")
