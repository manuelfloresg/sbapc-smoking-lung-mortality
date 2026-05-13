# runs/run_v7_seed4.R
source("runs/_runtime_setup.R")
source("runs/_source_all.R")
source("adapters/build_inputs_sim.R")
source("R/09_figures_maintext.R")
source("runs/replication_diagnostics.R")

# Override output path for V14
OUT_BASE <<- "results/20260507_0153_STABLE_V14"
OUT_SEC4 <<- file.path(OUT_BASE, "section4")
OUT_APPENDIX <<- file.path(OUT_BASE, "appendixC")
OUT_RAW <<- file.path(OUT_BASE, "raw_data")

dir.create(OUT_SEC4, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_APPENDIX, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_RAW, recursive = TRUE, showWarnings = FALSE)

# Run Seed 4 only
run_single_seed_replication(seed = 4, dgp = "spec_linear", force_rerun = TRUE)

# Generate Parity plots (Fig 2 equivalent)
message("Generating Parity plots...")
# Parity for each scenario
for (scen in c("freeze", "up1pc", "down1pc", "quit")) {
  g <- plot_deconstruction_figure(seed = 4, dgp = "spec_linear", scen = scen)
  ggsave(file.path(OUT_SEC4, sprintf("fig_parity_s4_%s.png", scen)), g, width = 10, height = 6, bg = "white")
}

# Generate Waterfall
g2 <- plot_transmission_waterfall(seed = 4, dgp = "spec_linear", scen = "quit")
ggsave(file.path(OUT_SEC4, "fig_waterfall_seed4.png"), g2, width = 8, height = 10, bg = "white")

message("V7 Seed 4 diagnostic completed in: ", OUT_BASE)
