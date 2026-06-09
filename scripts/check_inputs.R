repo_root <- normalizePath(file.path(dirname(sys.frame(1)$ofile), ".."), winslash = "/", mustWork = FALSE)

defaults <- c(
  BAPC_PATH_MORT_CSV = file.path(repo_root, "data", "analysis_ready", "uruguay_mortality_smooth_cancer.csv"),
  BAPC_PATH_POP_DTA  = file.path(repo_root, "data", "analysis_ready", "uruguay_population_1950_2070.dta"),
  BAPC_PATH_PREV_DTA = file.path(repo_root, "data", "analysis_ready", "uruguay_smoking_prevalence_aggregated.csv"),
  BAPC_PATH_INC_CSV  = file.path(repo_root, "data", "analysis_ready", "uruguay_incidence_smooth_1998_2022.csv")
)

status <- data.frame(
  variable = names(defaults),
  path = vapply(names(defaults), function(nm) {
    val <- Sys.getenv(nm, unset = "")
    if (nzchar(val)) val else defaults[[nm]]
  }, character(1)),
  stringsAsFactors = FALSE
)
status$exists <- file.exists(status$path)

print(status, row.names = FALSE)
if (!all(status$exists)) {
  stop("Some empirical inputs are missing. Place files in data/analysis_ready/ or set BAPC_PATH_* variables.")
}
