# Final 200-seed simulation run for Section 4 and Appendix C.
#
# Intended RStudio use:
#   setwd("D:/Git/Bloomberg_2025/Estimacion4")
#   source("runs/run_final_simulation_200.R")
#
# The script writes only under Estimacion4/results and keeps INLA temporaries
# in C:/tmp_inla, never in Dropbox.

Sys.setenv(
  BAPC_OUT_BASE = Sys.getenv("BAPC_OUT_BASE", "results/20260521_FINAL_200SEEDS"),
  BAPC_N_SEEDS = Sys.getenv("BAPC_N_SEEDS", "200"),
  BAPC_FIG_FORMAT = Sys.getenv("BAPC_FIG_FORMAT", "both"),
  INLA_TMPDIR = "C:/tmp_inla",
  TMPDIR = "C:/tmp_inla",
  TMP = "C:/tmp_inla",
  TEMP = "C:/tmp_inla",
  OMP_NUM_THREADS = "1"
)

dir.create("C:/tmp_inla", recursive = TRUE, showWarnings = FALSE)

source("runs/replication_diagnostics.R")

final_seeds <- seq_len(as.integer(Sys.getenv("BAPC_N_SEEDS", "200")))
final_n_cores <- as.integer(Sys.getenv("BAPC_FINAL_N_CORES", "4"))
if (!is.finite(final_n_cores) || final_n_cores < 1L) final_n_cores <- 4L

replicate_final_simulations(
  seeds = final_seeds,
  n_cores = final_n_cores,
  force_rerun = FALSE,
  run_oracle = TRUE,
  run_misspec = TRUE
)

