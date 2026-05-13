
# =============================================================
# Runner for method-paper main outputs (first pass)
# =============================================================
options(error = function() {
  traceback(2)
  q("no", status = 1)
})

`%||%` <- function(x, y) if (is.null(x)) y else x
source(file.path(if(sys.nframe() > 0) dirname(sys.frame(1)$ofile) else getwd(), "_source_all.R"))
source(file.path(if(sys.nframe() > 0) dirname(sys.frame(1)$ofile) else getwd(), "run_real_lung.R"))

method_paper_paths <- function(run_cfg = NULL) {
  if (is.null(run_cfg)) run_cfg <- get("run_cfg", envir = .GlobalEnv)
  root <- file.path(BAPC_PATHS$results, "method_paper")
  list(
    root = root,
    figures = file.path(root, "figures", "main_text"),
    tables = file.path(root, "tables", "main_text"),
    registry = file.path(root, "figures", "figure_registry.csv")
  )
}

run_method_paper_empirical <- function(run_cfg = NULL,
                                       scenarios = NULL,
                                       show_ci = FALSE,
                                       export_panels = TRUE,
                                       export_figure = TRUE,
                                       export_diagnostics = TRUE,
                                       ref_scenario = "freeze") {
  if (is.null(run_cfg)) run_cfg <- get("run_cfg", envir = .GlobalEnv)
  if (is.null(scenarios)) scenarios <- run_cfg$scenario_set
  paths <- method_paper_paths(run_cfg)
  invisible(lapply(unname(paths[c("root", "figures", "tables")]), dir.create, recursive = TRUE, showWarnings = FALSE))
  runs <- run_empirical_lung_all_scenarios(run_cfg = run_cfg, scenarios = scenarios)
  panel_out <- NULL
  if (isTRUE(export_panels)) {
    panel_out <- export_empirical_lung_panel_set(
      run_list = runs,
      out_dir = paths$figures,
      scenarios = scenarios,
      show_ci = show_ci,
      registry_file = paths$registry,
      stem = "panel_empirical_lung"
    )
  }
  fig_out <- NULL
  if (isTRUE(export_figure)) {
    fig_out <- export_empirical_lung_main_figure(
      run_list = runs,
      out_dir = paths$figures,
      scenarios = scenarios,
      show_ci = show_ci,
      registry_file = paths$registry,
      stem = "fig_main_empirical_lung"
    )
  }
  diag_out <- NULL
  if (isTRUE(export_diagnostics)) {
    diag_out <- export_scenario_transmission_diagnostics(
      run_list = runs,
      out_dir = paths$tables,
      ref_scenario = ref_scenario,
      scenarios = scenarios,
      stem = "empirical_lung_scenario_transmission"
    )
  }
  invisible(list(runs = runs, paths = paths, panels = panel_out, figure = fig_out, diagnostics = diag_out))
}

if (sys.nframe() == 0L) {
  run_method_paper_empirical(show_ci = FALSE, export_figure = TRUE)
}
