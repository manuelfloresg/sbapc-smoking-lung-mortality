# runs/run_seed4_fix.R
source("runs/_runtime_setup.R")
source("runs/_source_all.R")
source("runs/replication_diagnostics.R")

# Override output directory to avoid overwriting EVERYTHING
OUT_BASE <- "results/20260515_FIX_SEED4"
OUT_RAW  <- file.path(OUT_BASE, "raw_data")
OUT_SEC4 <- file.path(OUT_BASE, "section4")
dir.create(OUT_RAW, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_SEC4, recursive = TRUE, showWarnings = FALSE)

# Run Seed 4
message("Running Seed 4 with FREEZE method...")
run_single_seed_replication(seed = 4, dgp = "spec_linear", scens = c("freeze", "up1pc", "quit"), force_rerun = TRUE)

# Generate Section 4 plots for Seed 4
message("Generating plots...")
# Helper to point plot functions to the new OUT_RAW
# (Internal to replication_diagnostics.R, but we need to ensure it uses our OUT_BASE)
# Since replication_diagnostics.R uses global variables, we just update them.
# Note: In a real script we would pass paths, but here we leverage the orchestrator's structure.

# 1. Deconstruction
g1 <- plot_deconstruction_figure(seed = 4, dgp = "spec_linear", scen = "up1pc")
ggsave(file.path(OUT_SEC4, "fig_deconstruction_seed4_up1pc_FIX.png"), g1, width = 10, height = 6, bg = "white")

# 2. Sensitivity
g_sens <- plot_scenario_sensitivity_informed(seed = 4, dgp = "spec_linear")
ggsave(file.path(OUT_SEC4, "fig_sensitivity_seed4_FIX.png"), g_sens, width = 10, height = 6, bg = "white")

message("Done. Results in results/20260515_FIX_SEED4/section4")
