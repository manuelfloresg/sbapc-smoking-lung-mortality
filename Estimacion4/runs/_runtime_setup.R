# =============================================================
# Runtime setup for project execution
# =============================================================
# .libPaths(c("C:/Users/Manuel/AppData/Local/R/win-library/4.5", .libPaths()))
# Ensure we include the 4.5 library if it exists and we are on 4.5.x
if (R.version$major == "4" && grepl("^4\\.5", R.version$minor)) {
  u_lib <- "C:/Users/Manuel/AppData/Local/R/win-library/4.5"
  if (dir.exists(u_lib)) .libPaths(unique(c(u_lib, .libPaths())))
}


if (!exists("project_root")) {
  project_root <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

# 0) Configuration and Paths
# =============================================================

# Prioritize Dropbox for results if available
dropbox_results <- "d:/Dropbox/Investigacion/Bloomberg_2025/Estimacion4/results"
default_results <- file.path(project_root, "results")

BAPC_PATHS <- list(
  project_root = project_root,
  runtime      = file.path(project_root, "runtime"),
  results      = if (dir.exists(dropbox_results)) dropbox_results else default_results,
  inla_tmp     = "C:/tmp_inla"
)

# Only the truly project-local folders live inside the project.
# INLA tmp is external on purpose.
invisible(lapply(
  unname(BAPC_PATHS[c("runtime", "results")]),
  dir.create, recursive = TRUE, showWarnings = FALSE
))
dir.create(BAPC_PATHS$inla_tmp, recursive = TRUE, showWarnings = FALSE)

# 1) Packages and options
# =============================================================

required_pkgs <- c(
  "sn","dplyr","tidyr","readr","tibble","purrr","forecast","urca","tseries",
  "stringr","haven","mgcv","ggplot2","scales","patchwork","memoise","INLA","stringi"
)
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) {
  stop(
    sprintf(
      "Missing required packages: %s. Install them before running the project.",
      paste(missing_pkgs, collapse = ", ")
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
# Setting this here ensures every script uses the thread-safe shared tmp
options(INLA.tmpdir = BAPC_PATHS$inla_tmp)
INLA::inla.setOption(working.directory = BAPC_PATHS$inla_tmp)

# Avoid multi-threading issues in parallel workers
INLA::inla.setOption(num.threads = 1L)
# For the DGP tanh logic
Sys.setenv(OMP_NUM_THREADS = "1")

# 2) Environment variables (Persistent Defaults)
# =============================================================
# Hard-coded defaults for this workspace to ensure zero-config execution
options(BAPC_PATH_MORT_CSV = "d:/Dropbox/Investigacion/Bloomberg_2025/Mortalidad/muertes_suavizadas_cancer.csv")
options(BAPC_PATH_POP_DTA  = "d:/Dropbox/Investigacion/Bloomberg_2025/Base de datos/Proyecciones población/poblacion_1950_2070_empalmada.dta")
options(BAPC_PATH_PREV_DTA = "d:/Dropbox/Investigacion/Bloomberg_2025/Base de datos/base_completa.dta")
options(BAPC_PATH_INC_CSV  = "d:/Dropbox/Investigacion/Bloomberg_2025/Resultados/incidencia_suavizada_1998_2022.csv")

message("BAPC project_root : ", project_root)
message("BAPC runtime      : ", BAPC_PATHS$runtime)
message("BAPC results      : ", BAPC_PATHS$results)
message("BAPC INLA tmpdir  : ", BAPC_PATHS$inla_tmp)
message("INLA work dir     : ", INLA::inla.getOption("working.directory"))
