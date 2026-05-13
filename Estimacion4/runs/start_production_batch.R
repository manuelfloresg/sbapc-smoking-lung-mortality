
# runs/start_production_batch.R
source("runs/replication_diagnostics.R")
# Run 50 seeds in parallel
run_simulation_replication(seeds = 1:50, n_cores = 8, force_rerun = TRUE)
