# =============================================================
# Run empirical pipeline for lung cancer only
# Same estimation architecture; distinct scenario set lives in config only
# =============================================================

`%||%` <- function(x, y) if (is.null(x)) y else x
source(file.path(dirname(sys.frame(1)$ofile %||% getwd()), "_source_all.R"))

.extract_scalar <- function(x, default = NULL) {
  if (is.null(x)) return(default)
  if (is.list(x) && length(x) == 1L) return(x[[1]])
  x
}

.extract_intvec <- function(x) {
  if (is.null(x)) return(integer(0))
  if (is.list(x) && length(x) == 1L) x <- x[[1]]
  if (length(x) == 0) return(integer(0))
  as.integer(x)
}

scenario_paths_lung <- function(run_cfg, scenario_name) {
  root <- file.path(run_cfg$out_root, scenario_name)
  shared_root <- file.path(run_cfg$out_root, "_shared")
  list(
    root = root,
    result_root = root,
    master = file.path(root, "master"),
    diagnostics = file.path(root, "diagnostics"),
    shared_root = shared_root,
    shared_common = file.path(shared_root, "common"),
    shared_result_root = shared_root
  )
}

cfg <- dplyr::filter(causes, .data$cause_id == "lung")
if (nrow(cfg) != 1) stop("No pude identificar una única fila para 'lung' en causes.")

run_cfg <- list(
  causes_tbl         = cfg,
  scenario_set       = SCENARIOS_METHOD_LUNG,
  default_scenario   = "freeze",
  out_root           = file.path(BAPC_PATHS$results, "lung_method"),
  aggregate_outputs  = FALSE
)

run_pipeline_both_cause_lung <- function(cfg_row, prev_cfg, run_cfg = NULL) {
  if (is.null(run_cfg)) run_cfg <- get("run_cfg", envir = parent.frame())
  cfg_row <- tibble::as_tibble(cfg_row)
  stopifnot(nrow(cfg_row) == 1)

  inputs <- build_inputs_real_cause(cfg_row)

  run_pipeline_both_from_inputs(
    inputs = inputs,
    cfg_row = cfg_row,
    prev_cfg = prev_cfg
  )
}

run_empirical_lung <- function(run_cfg = NULL, scenario_name = NULL, prev_cfg = NULL) {
  if (is.null(run_cfg)) run_cfg <- get("run_cfg", envir = .GlobalEnv)
  if (is.null(scenario_name)) scenario_name <- run_cfg$default_scenario
  scenario_name <- normalize_prev_scenario_name(as.character(scenario_name)[1])
  if (is.null(prev_cfg)) prev_cfg <- get_prev_config(scenario = scenario_name)
  paths <- scenario_paths_lung(run_cfg, scenario_name)
  invisible(lapply(unname(paths[c("root","result_root","master","diagnostics","shared_root","shared_common","shared_result_root")]),
                   dir.create, recursive = TRUE, showWarnings = FALSE))

  cfg <- run_cfg$causes_tbl[1, ]
  cat(">>> [", scenario_name, "] Running cause:", cfg$cause_id[[1]], "-", cfg$label[[1]], "\n", sep = "")
  res <- run_pipeline_both_cause_lung(cfg, prev_cfg = prev_cfg, run_cfg = run_cfg)

  attr(res, "cause_id") <- cfg$cause_id[[1]]
  attr(res, "label")    <- cfg$label[[1]]
  attr(res, "scenario") <- scenario_name

  maybe_save_prev_apc_global(res, paths$shared_common)

  save_all_outputs(
    res,
    cfg$cause_id[[1]],
    cfg$label[[1]],
    out_base = paths$result_root,
    static_out_base = paths$shared_result_root,
    flatten_single_cause = TRUE
  )
  params_tbl <- pack_params(res, cfg$cause_id[[1]], cfg$label[[1]]) |> dplyr::mutate(scenario = scenario_name, .before = 1)
  proj_tbl   <- pack_proj(res,   cfg$cause_id[[1]], cfg$label[[1]]) |> dplyr::mutate(scenario = scenario_name, .before = 1)
  horizon_tbl <- pack_horizon(res, cfg$cause_id[[1]], cfg$label[[1]]) |> dplyr::mutate(scenario = scenario_name, .before = 1)

  readr::write_csv(params_tbl, file.path(paths$master, "params_lung.csv"))
  readr::write_csv(proj_tbl,   file.path(paths$master, "projections_lung.csv"))
  readr::write_csv(horizon_tbl, file.path(paths$master, "projection_horizon_lung.csv"))

  qc_tbl <- write_qc_outputs(params_tbl, proj_tbl,
                             out_file = file.path(paths$diagnostics, "qc_flags_lung.csv"),
                             print_top = TRUE,
                             top_n = 6)

  invisible(list(scenario = scenario_name, prev_cfg = prev_cfg, res = res, params_tbl = params_tbl, proj_tbl = proj_tbl, horizon_tbl = horizon_tbl, qc_tbl = qc_tbl, paths = paths))
}

run_empirical_lung_all_scenarios <- function(run_cfg = NULL, scenarios = NULL) {
  if (is.null(run_cfg)) run_cfg <- get("run_cfg", envir = .GlobalEnv)
  if (is.null(scenarios)) scenarios <- run_cfg$scenario_set
  out <- setNames(vector("list", length(scenarios)), scenarios)
  for (scn in scenarios) out[[scn]] <- run_empirical_lung(run_cfg = run_cfg, scenario_name = scn)
  invisible(out)
}

if (sys.nframe() == 0L) {
  run_empirical_lung(run_cfg = run_cfg, scenario_name = run_cfg$default_scenario)
}
