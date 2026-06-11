repo_root <- normalizePath(file.path(dirname(sys.frame(1)$ofile), ".."), winslash = "/", mustWork = FALSE)
project_dir <- file.path(repo_root, "code")

if (!dir.exists(project_dir)) stop("Cannot find code directory: ", project_dir)

if (!nzchar(Sys.getenv("BAPC_URUGUAY_OUT_BASE", unset = ""))) {
  Sys.setenv(BAPC_URUGUAY_OUT_BASE = "results/20260520_URUGUAY_CANDIDATE")
}
if (!nzchar(Sys.getenv("BAPC_INLA_TMPDIR", unset = ""))) {
  Sys.setenv(BAPC_INLA_TMPDIR = "C:/tmp_inla")
}

source(file.path(repo_root, "scripts", "check_inputs.R"))

old_wd <- getwd()
on.exit(setwd(old_wd), add = TRUE)
setwd(project_dir)
source("runs/uruguay_products.R")
replicate_uruguay_empirical(save_raw_rds = FALSE, run_multisite = TRUE)
