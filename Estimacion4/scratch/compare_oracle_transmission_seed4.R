setwd(normalizePath(file.path(getwd(), "Estimacion4"), winslash = "/", mustWork = TRUE))

source("runs/replication_diagnostics.R")

summarize_map <- function(raw_dir, label) {
  g <- plot_transmission_map(
    seed = 4,
    dgp = "spec_linear",
    sex_lab = "M",
    raw_dir = raw_dir,
    title_suffix = label
  )

  tibble::as_tibble(g$data) |>
    dplyr::mutate(period_block = dplyr::if_else(period <= 2022, "historical", "future")) |>
    tidyr::pivot_wider(names_from = series, values_from = value) |>
    dplyr::mutate(
      diff = .data[["Informed SBAPC"]] - .data[["Truth"]],
      rel = diff / pmax(abs(.data[["Truth"]]), 1e-9),
      label = label
    ) |>
    dplyr::group_by(label, metric, period_block) |>
    dplyr::summarise(
      diff_mean = mean(diff, na.rm = TRUE),
      rel_mean = mean(rel, na.rm = TRUE),
      truth_mean = mean(.data[["Truth"]], na.rm = TRUE),
      informed_mean = mean(.data[["Informed SBAPC"]], na.rm = TRUE),
      .groups = "drop"
    )
}

raw_realistic <- file.path(OUT_BASE, "raw_data_realistic_rngfixed")
raw_oracle <- file.path(OUT_BASE, "raw_data_oracle_rngfixed")

out <- dplyr::bind_rows(
  summarize_map(raw_realistic, "window_limited"),
  summarize_map(raw_oracle, "full_support")
)

readr::write_csv(
  out,
  file.path(OUT_SEC4, "transmission_map_oracle_vs_realistic_seed4_M.csv")
)

print(out)
