# =============================================================
# Runtime setup for project execution
# =============================================================

# Make user-level R libraries visible across R 4.x patch/minor versions.
# R.version$minor is "5.3" for R 4.5.3, so do not test it as "4.5".
.candidate_user_libs <- unique(c(
  Sys.getenv("BAPC_R_LIB", unset = ""),
  Sys.getenv("R_LIBS_USER", unset = ""),
  file.path(Sys.getenv("LOCALAPPDATA", unset = ""), "R", "win-library", c("4.6", "4.5", "4.4"))
))
.candidate_user_libs <- .candidate_user_libs[nzchar(.candidate_user_libs)]
.candidate_user_libs <- .candidate_user_libs[dir.exists(.candidate_user_libs)]
if (length(.candidate_user_libs)) {
  .libPaths(unique(c(.candidate_user_libs, .libPaths())))
}

if (!exists("project_root")) {
  project_root <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

# Prefer a project-local package library when present. This keeps replication
# dependencies out of Dropbox and out of system-level R folders.
.r_minor <- paste(R.version$major, strsplit(R.version$minor, ".", fixed = TRUE)[[1]][1], sep = ".")
.project_lib_roots <- unique(c(
  file.path(project_root, ".Rlib"),
  file.path(dirname(project_root), ".Rlib")
))
.project_libs <- as.vector(outer(.project_lib_roots, c(.r_minor, R.version$major), file.path))
.project_libs <- .project_libs[dir.exists(.project_libs)]
if (length(.project_libs)) {
  .libPaths(unique(c(.project_libs, .libPaths())))
}

# 0) Configuration and Paths
# =============================================================

default_results <- file.path(project_root, "results")
results_override <- Sys.getenv("BAPC_RESULTS_DIR", unset = "")
if (!nzchar(results_override)) {
  results_override <- getOption("BAPC_RESULTS_DIR", default_results)
}

inla_tmp_override <- Sys.getenv("BAPC_INLA_TMPDIR", unset = "")
if (!nzchar(inla_tmp_override)) {
  inla_tmp_override <- getOption("BAPC_INLA_TMPDIR", "C:/tmp_inla")
}

BAPC_PATHS <- list(
  project_root = project_root,
  runtime      = file.path(project_root, "runtime"),
  results      = results_override,
  inla_tmp     = inla_tmp_override
)

# Only the truly project-local folders live inside the project.
# INLA tmp is external on purpose.
invisible(lapply(
  unname(BAPC_PATHS[c("runtime", "results")]),
  dir.create, recursive = TRUE, showWarnings = FALSE
))
dir.create(BAPC_PATHS$inla_tmp, recursive = TRUE, showWarnings = FALSE)

# Force all temp-related paths away from Dropbox before loading INLA.
Sys.setenv(
  INLA_TMPDIR = BAPC_PATHS$inla_tmp,
  TMPDIR      = BAPC_PATHS$inla_tmp,
  TMP         = BAPC_PATHS$inla_tmp,
  TEMP        = BAPC_PATHS$inla_tmp
)

# 1) Packages and options
# =============================================================

required_pkgs <- c(
  "sn","dplyr","tidyr","readr","tibble","purrr","forecast","urca","tseries",
  "stringr","haven","mgcv","ggplot2","scales","patchwork","memoise","INLA",
  "stringi","future","future.apply","svglite"
)
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) {
  stop(
    sprintf(
      "Missing required packages: %s. Install them before running the project. Current .libPaths(): %s",
      paste(missing_pkgs, collapse = ", "),
      paste(.libPaths(), collapse = " | ")
    )
  )
}

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(tibble); library(purrr)
  library(forecast); library(urca); library(tseries); library(stringr); library(haven)
  library(mgcv); library(ggplot2); library(scales); library(patchwork); library(memoise)
  library(INLA); library(stringi)
})

# ---------------------------
# GLOBAL INLA CONFIGURATION
# ---------------------------
# Setting this here ensures every script uses the thread-safe shared tmp.
options(INLA.tmpdir = BAPC_PATHS$inla_tmp)
INLA::inla.setOption(working.directory = BAPC_PATHS$inla_tmp)

# Avoid multi-threading issues in parallel workers.
INLA::inla.setOption(num.threads = 1L)
Sys.setenv(OMP_NUM_THREADS = "1")

# 2) Real-data path defaults
# =============================================================
# Simulation replication is self-contained and does not require empirical
# Uruguay files. Empirical replication uses analysis-ready inputs in the public
# repository by default. Alternatively, set the BAPC_PATH_* environment
# variables or R options before sourcing this file.
.repo_root <- normalizePath(file.path(project_root, ".."), winslash = "/", mustWork = FALSE)
.public_data_dir <- file.path(.repo_root, "data", "analysis_ready")
if (is.null(getOption("BAPC_PATH_MORT_CSV"))) {
  options(BAPC_PATH_MORT_CSV = file.path(.public_data_dir, "uruguay_mortality_smooth_cancer.csv"))
}
if (is.null(getOption("BAPC_PATH_POP_DTA"))) {
  options(BAPC_PATH_POP_DTA = file.path(.public_data_dir, "uruguay_population_1950_2070.dta"))
}
if (is.null(getOption("BAPC_PATH_PREV_DTA"))) {
  options(BAPC_PATH_PREV_DTA = file.path(.public_data_dir, "uruguay_smoking_prevalence_aggregated.csv"))
}
if (is.null(getOption("BAPC_PATH_INC_CSV"))) {
  options(BAPC_PATH_INC_CSV = file.path(.public_data_dir, "uruguay_incidence_smooth_1998_2022.csv"))
}

message("BAPC project_root : ", project_root)
message("BAPC runtime      : ", BAPC_PATHS$runtime)
message("BAPC results      : ", BAPC_PATHS$results)
message("BAPC INLA tmpdir  : ", BAPC_PATHS$inla_tmp)
message("INLA work dir     : ", INLA::inla.getOption("working.directory"))
