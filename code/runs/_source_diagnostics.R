`%||%` <- function(x, y) if (is.null(x)) y else x
source(file.path(dirname(sys.frame(1)$ofile %||% getwd()), "_source_all.R"))

project_root <- normalizePath(
  file.path(dirname(sys.frame(1)$ofile %||% getwd()), ".."),
  winslash = "/",
  mustWork = FALSE
)
if (!dir.exists(project_root)) project_root <- getwd()

diag_files <- c(
  "11_diagnostics_mortality_channel.R",
  "12_diagnostics_mortality_design.R",
  "13_diagnostics_mortality_wiring.R",
  "14_diagnostics_mortality_objectpath.R",
  "15_diagnostics_mortality_age.R",
  "16_diagnostics_mortality_celltrace.R",
  "17_diagnostics_LI_path.R",
  "18_diagnostics_mortality_predictor_decomp.R",
  "19_diagnostics_full_chain_trace.R",
  "20_diagnostics_incidence_eta_trace.R",
  "21_diagnostics_mortality_truth_trace.R",
  "22_diagnostics_mortality_kernel_scale.R",
  "23_inspect_mortality_input_construction.R",
  "24_inspect_attach_external_mortality_offset.R",
  "25_inspect_get_postdx_kernel.R",
  "26_inspect_incidence_to_cases_scale.R",
  "27_audit_units_chain.R",
  "28_patch_mort_truth_from_inc_kernel.R"
)

for (ff in diag_files) {
  path_ff <- file.path(project_root, "R", ff)
  if (file.exists(path_ff)) {
    try(source(path_ff), silent = TRUE)
  }
}
