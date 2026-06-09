# Required packages for the SBAPC replication repository.

cran_packages <- c(
  "sn",
  "dplyr",
  "tidyr",
  "readr",
  "tibble",
  "purrr",
  "forecast",
  "urca",
  "tseries",
  "stringr",
  "haven",
  "mgcv",
  "ggplot2",
  "scales",
  "patchwork",
  "memoise",
  "stringi",
  "future",
  "future.apply",
  "svglite"
)

inla_package <- "INLA"
inla_repository <- "https://inla.r-inla-download.org/R/stable"

missing_cran <- cran_packages[!vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_cran)) {
  install.packages(missing_cran, dependencies = TRUE)
}

if (!requireNamespace(inla_package, quietly = TRUE)) {
  install.packages(inla_package, repos = c(getOption("repos"), INLA = inla_repository))
}
