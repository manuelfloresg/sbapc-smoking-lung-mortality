repo_root <- normalizePath(file.path(dirname(sys.frame(1)$ofile), ".."), winslash = "/", mustWork = FALSE)
project_dir <- file.path(repo_root, "code")

if (!dir.exists(project_dir)) stop("Cannot find code directory: ", project_dir)

if (!nzchar(Sys.getenv("BAPC_OUT_BASE", unset = ""))) {
  Sys.setenv(BAPC_OUT_BASE = "results/20260521_FINAL_200SEEDS")
}
if (!nzchar(Sys.getenv("BAPC_N_SEEDS", unset = ""))) {
  Sys.setenv(BAPC_N_SEEDS = "200")
}
if (!nzchar(Sys.getenv("BAPC_FINAL_N_CORES", unset = ""))) {
  Sys.setenv(BAPC_FINAL_N_CORES = "4")
}
if (!nzchar(Sys.getenv("BAPC_FIG_FORMAT", unset = ""))) {
  Sys.setenv(BAPC_FIG_FORMAT = "both")
}
if (!nzchar(Sys.getenv("BAPC_INLA_TMPDIR", unset = ""))) {
  Sys.setenv(BAPC_INLA_TMPDIR = "C:/tmp_inla")
}

old_wd <- getwd()
on.exit(setwd(old_wd), add = TRUE)
setwd(project_dir)
source("runs/run_final_simulation_200.R")
