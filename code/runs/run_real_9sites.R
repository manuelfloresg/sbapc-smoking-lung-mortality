`%||%` <- function(x, y) if (is.null(x)) y else x

# =============================================================
# Run empirical pipeline for the full set of causes
# Third cut: explicit scenario config + scenario-aware outputs
# =============================================================

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

scenario_paths_9sites <- function(run_cfg, scenario_name) {
  root <- file.path(run_cfg$out_root, scenario_name)
  shared_root <- file.path(run_cfg$out_root, "_shared")
  list(
    root = root,
    by_cause = file.path(root, "by_cause"),
    master = file.path(root, "master"),
    aggregate = file.path(root, "aggregate"),
    diagnostics = file.path(root, "diagnostics"),
    shared_root = shared_root,
    shared_common = file.path(shared_root, "common"),
    shared_by_cause = file.path(shared_root, "by_cause")
  )
}

run_cfg <- list(
  causes_tbl         = causes,
  scenario_set       = SCENARIOS_REAL_9SITES,
  default_scenario   = "freeze",
  out_root           = file.path(BAPC_PATHS$results, "real_9sites"),
  aggregate_outputs  = TRUE
)

run_pipeline_both_cause_9sites <- function(cfg_row, prev_cfg, run_cfg = NULL) {
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


validate_external_kernel_result <- function(res, cause_id = NA_character_) {
  .check_one <- function(res_sex, sex_tag) {
    if (is.null(res_sex)) return(invisible(NULL))
    mort_link_mode <- tryCatch(res_sex$params$mort_link_mode, error = function(e) NA_character_)
    if (!identical(mort_link_mode, "external_kernel")) {
      stop(sprintf("[%s|%s] mort_link_mode inesperado: %s", cause_id, sex_tag, as.character(mort_link_mode)), call. = FALSE)
    }

    max_lag <- tryCatch(res_sex$params$mort_kernel_max_lag, error = function(e) NA_integer_)
    if (!isTRUE(is.finite(max_lag)) || as.integer(max_lag) != 5L) {
      stop(sprintf("[%s|%s] mort_kernel_max_lag inesperado: %s", cause_id, sex_tag, as.character(max_lag)), call. = FALSE)
    }

    kernel_tbl <- tryCatch(res_sex$diag$mort_kernel, error = function(e) NULL)
    if (!is.data.frame(kernel_tbl) || !nrow(kernel_tbl)) {
      stop(sprintf("[%s|%s] diag$mort_kernel ausente o vacío.", cause_id, sex_tag), call. = FALSE)
    }

    invisible(NULL)
  }

  .check_one(tryCatch(res$resM, error = function(e) NULL), "M")
  .check_one(tryCatch(res$resF, error = function(e) NULL), "F")
  invisible(TRUE)
}

run_empirical_9sites <- function(run_cfg = NULL, scenario_name = NULL, prev_cfg = NULL) {
  if (is.null(run_cfg)) run_cfg <- get("run_cfg", envir = .GlobalEnv)
  if (is.null(scenario_name)) scenario_name <- run_cfg$default_scenario
  scenario_name <- normalize_prev_scenario_name(as.character(scenario_name)[1])
  if (is.null(prev_cfg)) prev_cfg <- get_prev_config(scenario = scenario_name)
  paths <- scenario_paths_9sites(run_cfg, scenario_name)
  invisible(lapply(unname(paths[c("root","by_cause","master","aggregate","diagnostics","shared_root","shared_common","shared_by_cause")]),
                   dir.create, recursive = TRUE, showWarnings = FALSE))

  params_tbl <- tibble::tibble()
  proj_tbl   <- tibble::tibble()
  horizon_tbl <- tibble::tibble()

  for (i in seq_len(nrow(run_cfg$causes_tbl))) {
    cfg <- run_cfg$causes_tbl[i, ]
    cat(">>> [", scenario_name, "] Running cause:", cfg$cause_id[[1]], "-", cfg$label[[1]], "\n", sep = "")
    res <- run_pipeline_both_cause_9sites(cfg, prev_cfg = prev_cfg, run_cfg = run_cfg)

    attr(res, "cause_id") <- cfg$cause_id[[1]]
    attr(res, "label")    <- cfg$label[[1]]
    attr(res, "scenario") <- scenario_name

    validate_external_kernel_result(res, cause_id = cfg$cause_id[[1]])

    if (i == 1L) maybe_save_prev_apc_global(res, paths$shared_common)

    save_all_outputs(res, cfg$cause_id[[1]], cfg$label[[1]], out_base = paths$by_cause, static_out_base = paths$shared_by_cause)
    params_tbl <- dplyr::bind_rows(params_tbl, pack_params(res, cfg$cause_id[[1]], cfg$label[[1]]) |> dplyr::mutate(scenario = scenario_name, .before = 1))
    proj_tbl   <- dplyr::bind_rows(proj_tbl,   pack_proj(res,   cfg$cause_id[[1]], cfg$label[[1]]) |> dplyr::mutate(scenario = scenario_name, .before = 1))
    horizon_tbl <- dplyr::bind_rows(horizon_tbl, pack_horizon(res, cfg$cause_id[[1]], cfg$label[[1]]) |> dplyr::mutate(scenario = scenario_name, .before = 1))
  }

  readr::write_csv(params_tbl, file.path(paths$master, "params_by_cause.csv"))
  readr::write_csv(proj_tbl,   file.path(paths$master, "projections_by_cause.csv"))
  readr::write_csv(horizon_tbl, file.path(paths$master, "projection_horizon_by_cause.csv"))

  if (isTRUE(run_cfg$aggregate_outputs) && nrow(run_cfg$causes_tbl) > 1) {
    proj_total_cancers <- agregar_todas_causas(proj_tbl, method = "normal")
    readr::write_csv(proj_total_cancers, file.path(paths$aggregate, "projections_total_cancers.csv"))
  }

  qc_tbl <- write_qc_outputs(params_tbl, proj_tbl,
                             out_file = file.path(paths$diagnostics, "qc_flags_by_cause.csv"),
                             print_top = TRUE)

  invisible(list(scenario = scenario_name, prev_cfg = prev_cfg, params_tbl = params_tbl, proj_tbl = proj_tbl, horizon_tbl = horizon_tbl, qc_tbl = qc_tbl, paths = paths))
}

run_empirical_9sites_all_scenarios <- function(run_cfg = NULL, scenarios = NULL) {
  if (is.null(run_cfg)) run_cfg <- get("run_cfg", envir = .GlobalEnv)
  if (is.null(scenarios)) scenarios <- run_cfg$scenario_set
  out <- setNames(vector("list", length(scenarios)), scenarios)
  for (scn in scenarios) out[[scn]] <- run_empirical_9sites(run_cfg = run_cfg, scenario_name = scn)
  invisible(out)
}

if (sys.nframe() == 0L) {
  run_empirical_9sites(run_cfg = run_cfg, scenario_name = run_cfg$default_scenario)
}
