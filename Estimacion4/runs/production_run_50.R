# production_run_50.R
source("runs/_runtime_setup.R")
source("runs/replication_diagnostics.R")

# Ensure results directory is clean for a fresh run
OUT_BASE    <- "results/20260515_FINAL_PROD"
OUT_RAW     <- file.path(OUT_BASE, "raw_data")
# dir.create(OUT_RAW, recursive = TRUE, showWarnings = FALSE) # replication_diagnostics creates it

message(">>> STARTING FULL PRODUCTION RUN (50 SEEDS) <<<")
message("Time: ", Sys.time())

# We use 6 cores as recommended in the summary
run_simulation_replication(n_cores = 6, seeds = 1:50)

message(">>> PRODUCTION RUN COMPLETED <<<")
message("Time: ", Sys.time())
