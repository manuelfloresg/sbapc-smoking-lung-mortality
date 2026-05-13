`%||%` <- function(x, y) if (is.null(x)) y else x

# Main project loader for the current replication-ready workflow.
if (!exists("project_root") || !dir.exists(project_root)) {
  project_root <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

# 1) Core environment setup
# =============================================================
# This script must be sourced first to define BAPC_PATHS and load libraries
if (!exists("BAPC_PATHS")) {
  source(file.path(project_root, "runs", "_runtime_setup.R"))
}

# Core defaults and helpers
source(file.path(project_root, "R", "00_defaults.R"))
source(file.path(project_root, "R", "01_core_helpers.R"))
source(file.path(project_root, "R", "01b_engine_consistency_helpers.R"))
source(file.path(project_root, "R", "01c_prediction_rebuild_helpers.R"))

# Data-input builders
source(file.path(project_root, "adapters", "build_inputs_real.R"))
source(file.path(project_root, "adapters", "build_inputs_sim.R"))

# Main modeling pipeline
source(file.path(project_root, "R", "02_stage_models.R"))
source(file.path(project_root, "R", "03_pipeline_sex.R"))
source(file.path(project_root, "R", "04_pipeline_both.R"))
source(file.path(project_root, "R", "05_postprocess.R"))
source(file.path(project_root, "R", "06_qc.R"))

# 3) UI/Plotting and Diagnostics helpers
# =============================================================
source(file.path(project_root, "R", "31_diagnostics_against_truth.R"))


message(">>> SBAPC Environment Loaded Successfully.")
