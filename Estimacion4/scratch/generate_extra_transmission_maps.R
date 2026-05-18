setwd(normalizePath(file.path(getwd(), "Estimacion4"), winslash = "/", mustWork = TRUE))

source("runs/replication_diagnostics.R")

extra_seeds <- c(4, 1, 13, 21, 32, 36, 50)
sex_lab <- "M"

summaries <- lapply(extra_seeds, function(seed) {
  message(sprintf("Generating transmission map for seed %d (%s)", seed, sex_lab))
  g <- plot_transmission_map(seed = seed, dgp = "spec_linear", sex_lab = sex_lab)
  save_paper_plot(
    g,
    file.path(OUT_SEC4, sprintf("fig_transmission_map_seed%d_%s", seed, sex_lab)),
    width = 13,
    height = 9,
    bg = "white"
  )

  g$data |>
    dplyr::mutate(
      seed = seed,
      period_block = dplyr::if_else(period <= 2022, "historical", "future")
    ) |>
    tidyr::pivot_wider(names_from = series, values_from = value) |>
    dplyr::mutate(
      diff_informed_minus_truth = `Informed SBAPC` - Truth,
      rel_diff_informed_minus_truth = diff_informed_minus_truth / pmax(abs(Truth), 1e-9)
    ) |>
    dplyr::group_by(seed, scenario, metric, period_block) |>
    dplyr::summarise(
      truth_mean = mean(Truth, na.rm = TRUE),
      informed_mean = mean(`Informed SBAPC`, na.rm = TRUE),
      diff_mean = mean(diff_informed_minus_truth, na.rm = TRUE),
      rel_diff_mean = mean(rel_diff_informed_minus_truth, na.rm = TRUE),
      .groups = "drop"
    )
})

summary_df <- dplyr::bind_rows(summaries)
readr::write_csv(
  summary_df,
  file.path(OUT_SEC4, sprintf("transmission_map_seed_summary_%s.csv", sex_lab))
)

print(summary_df)
