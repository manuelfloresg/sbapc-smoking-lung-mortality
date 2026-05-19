# Run missing simulation diagnostics for Section 4 / Appendix C products.
#
# This script intentionally keeps INLA parallelism moderate. It is meant for
# topping up an existing production-candidate run, not for force-overwriting it.

Sys.setenv(BAPC_OUT_BASE = Sys.getenv("BAPC_OUT_BASE", "results/20260518_PROD_CANDIDATE"))
source("runs/replication_diagnostics.R")

message("Output base: ", OUT_BASE)
message("INLA temp dir: ", Sys.getenv("INLA_TMPDIR", unset = Sys.getenv("TMPDIR", unset = "")))

missing_n_cores <- suppressWarnings(as.integer(Sys.getenv("BAPC_MISSING_N_CORES", "1")))
if (!is.finite(missing_n_cores) || missing_n_cores < 1L) missing_n_cores <- 1L
message("Top-up worker count: ", missing_n_cores)

oracle_missing <- setdiff(CANONICAL_SEEDS, available_result_seeds(OUT_RAW_ORACLE, dgp = "spec_linear"))
misspec_missing <- setdiff(CANONICAL_SEEDS, available_result_seeds(OUT_RAW, dgp = "misspec_tanh"))

message("Missing full-support spec_linear seeds: ", paste(oracle_missing, collapse = ", "))
message("Missing misspec_tanh seeds: ", paste(misspec_missing, collapse = ", "))

if (length(oracle_missing)) {
  run_simulation_replication(
    seeds = oracle_missing,
    dgps = "spec_linear",
    scens = CANONICAL_SCENS,
    force_rerun = FALSE,
    n_cores = missing_n_cores,
    information_set = "oracle",
    raw_dir = OUT_RAW_ORACLE
  )
}

if (length(misspec_missing)) {
  run_simulation_replication(
    seeds = misspec_missing,
    dgps = "misspec_tanh",
    scens = CANONICAL_SCENS,
    force_rerun = FALSE,
    n_cores = missing_n_cores,
    information_set = "realistic",
    raw_dir = OUT_RAW
  )
}

message("Regenerating Section 4 and Appendix C products after missing runs.")
data_main <- extract_all_metrics(dgps = "spec_linear")
generate_scenario_effect_products(data_main)
generate_appendix_c()
write_figure_titles_notes("section4")
write_float_inventories()

message("Missing simulation products completed.")
