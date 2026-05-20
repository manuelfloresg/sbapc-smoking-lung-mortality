# =============================================================
# Uruguay empirical products for Section 5 and Appendix D
# =============================================================

options(error = function() {
  traceback(2)
  q("no", status = 1)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

.this_file <- tryCatch(normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = FALSE),
                       error = function(e) NA_character_)
.run_dir <- if (is.na(.this_file) || !nzchar(.this_file)) getwd() else dirname(.this_file)
.project_root <- normalizePath(file.path(.run_dir, ".."), winslash = "/", mustWork = FALSE)
if (!dir.exists(file.path(.project_root, "R"))) .project_root <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)

if (!exists("project_root") || !dir.exists(project_root)) project_root <- .project_root

source(file.path(project_root, "runs", "run_real_lung.R"))
source(file.path(project_root, "R", "09_figures_maintext.R"))
source(file.path(project_root, "R", "10_diagnostics_methodpaper.R"))

URUGUAY_OUT_BASE <- {
  val <- Sys.getenv("BAPC_URUGUAY_OUT_BASE", unset = "")
  if (!nzchar(val)) val <- getOption("BAPC_URUGUAY_OUT_BASE", file.path(BAPC_PATHS$results, "20260520_URUGUAY_CANDIDATE"))
  normalizePath(val, winslash = "/", mustWork = FALSE)
}
OUT_SECTION5 <- file.path(URUGUAY_OUT_BASE, "section5")
OUT_APPENDIXD <- file.path(URUGUAY_OUT_BASE, "appendixD")
OUT_RAW_URUGUAY <- file.path(URUGUAY_OUT_BASE, "raw")
invisible(lapply(c(URUGUAY_OUT_BASE, OUT_SECTION5, OUT_APPENDIXD, OUT_RAW_URUGUAY),
                 dir.create, recursive = TRUE, showWarnings = FALSE))

URUGUAY_SCENARIOS <- c("freeze", "up1pc", "down1pc", "quit")
URUGUAY_SCEN_LABELS <- c(
  freeze = "Frozen at 2022",
  up1pc = "Up 1% p.a.",
  down1pc = "Down 1% p.a.",
  quit = "Quit"
)
URUGUAY_SCEN_LABELS_TEX <- c(
  freeze = "Frozen at 2022",
  up1pc = "$\\uparrow$1\\%",
  down1pc = "$\\downarrow$1\\%",
  quit = "Quit"
)
URUGUAY_COLORS <- c(
  "Frozen at 2022" = "#111111",
  "Up 1% p.a." = "#C0392B",
  "Down 1% p.a." = "#1E8449",
  "Quit" = "#2874A6"
)
URUGUAY_ZONE_LINETYPES <- c(
  credible = "solid",
  caution = "longdash",
  risky = "dotted",
  beyond_max = "dotdash",
  historical = "solid"
)

fmt_int <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  ifelse(is.finite(x), format(round(x), big.mark = ",", scientific = FALSE, trim = TRUE), "")
}

fmt_num <- function(x, digits = 1) {
  x <- suppressWarnings(as.numeric(x))
  ifelse(is.finite(x), formatC(x, digits = digits, format = "f"), "")
}

latex_open <- function(colspec) {
  c(
    "\\begingroup",
    "\\renewcommand{\\arraystretch}{1.08}",
    "\\setlength{\\tabcolsep}{4pt}",
    sprintf("\\begin{tabular}{%s}", colspec),
    "\\toprule"
  )
}

latex_close <- function() c("\\bottomrule", "\\end{tabular}", "\\endgroup")

save_uruguay_plot <- function(plot, path_no_ext, width, height, dpi = 300) {
  dir.create(dirname(path_no_ext), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(paste0(path_no_ext, ".svg"), plot = plot, width = width, height = height,
                  dpi = dpi, bg = "white")
  ggplot2::ggsave(paste0(path_no_ext, ".pdf"), plot = plot, width = width, height = height,
                  dpi = dpi, bg = "white", device = grDevices::cairo_pdf)
  invisible(path_no_ext)
}

sex_public <- function(x) dplyr::recode(as.character(x), M = "Male", F = "Female", T = "Total", .default = as.character(x))

metric_public <- function(x) dplyr::recode(as.character(x), incidence = "Incidence", mortality = "Mortality", .default = as.character(x))

main_series_for_metric <- function(metric) {
  if (identical(metric, "incidence")) "I|P" else "M|I|P"
}

benchmark_series_for_metric <- function(metric) {
  if (identical(metric, "incidence")) "I" else "BAPC benchmark"
}

extract_mortality_bapc <- function(run_out, scenario_name) {
  rows <- list()
  for (sx in c("M", "F")) {
    rs <- if (sx == "M") run_out$res$resM else run_out$res$resF
    if (is.null(rs) || is.null(rs$annual_bapc) || !nrow(rs$annual_bapc)) next
    rows[[sx]] <- tibble::as_tibble(rs$annual_bapc) |>
      dplyr::transmute(
        scenario = scenario_name,
        sex = sx,
        metric = "mortality",
        series = "BAPC benchmark",
        period = as.integer(period),
        mean = as.numeric(deaths_hat),
        lwr = as.numeric(deaths_lwr),
        upr = as.numeric(deaths_upr)
      )
  }
  dplyr::bind_rows(rows)
}

combined_projection_tables <- function(run_list, scenarios = names(run_list)) {
  proj <- dplyr::bind_rows(lapply(scenarios, function(scn) {
    run_list[[scn]]$proj_tbl |> dplyr::mutate(scenario = scn)
  }))
  mort_bapc <- dplyr::bind_rows(lapply(scenarios, function(scn) extract_mortality_bapc(run_list[[scn]], scn)))
  dplyr::bind_rows(proj, mort_bapc) |>
    dplyr::mutate(
      scenario = factor(as.character(scenario), levels = URUGUAY_SCENARIOS),
      scenario_label = factor(unname(URUGUAY_SCEN_LABELS[as.character(scenario)]),
                              levels = unname(URUGUAY_SCEN_LABELS[URUGUAY_SCENARIOS])),
      sex_label = factor(sex_public(sex), levels = c("Male", "Female", "Total")),
      metric_label = factor(metric_public(metric), levels = c("Incidence", "Mortality")),
      projection_zone = factor(as.character(projection_zone), levels = c("historical", "credible", "caution", "risky", "beyond_max"))
    )
}

observed_projection_tables <- function(run_list) {
  ref <- run_list[[names(run_list)[1]]]
  rows <- list()
  for (sx in c("M", "F")) {
    rs <- if (sx == "M") ref$res$resM else ref$res$resF
    rows[[paste0(sx, "_inc")]] <- rs$inc_obs_annual |>
      dplyr::transmute(sex = sx, metric = "incidence", period = as.integer(period), obs = as.numeric(obs))
    rows[[paste0(sx, "_mort")]] <- rs$obs_annual |>
      dplyr::transmute(sex = sx, metric = "mortality", period = as.integer(period), obs = as.numeric(obs))
  }
  dplyr::bind_rows(rows) |>
    dplyr::mutate(sex_label = factor(sex_public(sex), levels = c("Male", "Female")),
                  metric_label = factor(metric_public(metric), levels = c("Incidence", "Mortality")))
}

prepare_segmented_lines <- function(df) {
  dplyr::bind_rows(lapply(split(df, interaction(df$scenario, df$metric, df$sex, df$series, drop = TRUE)), function(dd) {
    dd <- dd[order(dd$period), , drop = FALSE]
    if (!nrow(dd)) return(dd)
    zone <- as.character(dd$projection_zone)
    zone[is.na(zone)] <- "beyond_max"
    grp <- cumsum(c(TRUE, zone[-1] != zone[-length(zone)]))
    dd$zone_group <- grp
    out <- list()
    gids <- sort(unique(grp))
    for (g in gids) {
      seg <- dd[dd$zone_group == g, , drop = FALSE]
      if (g > min(gids)) {
        prev_row <- dd[max(which(dd$zone_group < g)), , drop = FALSE]
        prev_row$zone_group <- g
        prev_row$projection_zone <- seg$projection_zone[[1]]
        seg <- dplyr::bind_rows(prev_row, seg)
      }
      out[[length(out) + 1L]] <- seg
    }
    dplyr::bind_rows(out)
  }))
}

horizon_region_from_year <- function(period, frontier_row) {
  cred <- suppressWarnings(as.integer(frontier_row$end_year_credible[1]))
  caut <- suppressWarnings(as.integer(frontier_row$end_year_caution[1]))
  risk <- suppressWarnings(as.integer(frontier_row$end_year_risky[1]))
  dplyr::case_when(
    is.finite(cred) & period <= cred ~ "Credible",
    is.finite(caut) & period <= caut ~ "Caution",
    is.finite(risk) & period <= risk ~ "Risky",
    TRUE ~ NA_character_
  )
}

get_total_horizon <- function(run_list) {
  h <- dplyr::bind_rows(lapply(names(run_list), function(scn) {
    run_list[[scn]]$horizon_tbl |> dplyr::mutate(scenario = scn)
  }))
  h |> dplyr::filter(as.character(sex) == "T") |> dplyr::slice(1)
}

run_uruguay_lung_rebuilt_all_scenarios <- function(run_cfg = get("run_cfg", envir = .GlobalEnv),
                                                   scenarios = URUGUAY_SCENARIOS,
                                                   save_raw_rds = FALSE) {
  cfg_row <- run_cfg$causes_tbl[1, ]
  inputs <- build_inputs_real_cause(cfg_row)
  extra_args <- list(trend_type = "level", gammaP_method = "freeze", sd_theta_IP = 2.0)

  message(">>> Uruguay empirical run: fitting freeze benchmark")
  res_freeze <- do.call(run_pipeline_both_from_inputs, c(list(
    inputs = inputs,
    cfg_row = cfg_row,
    prev_cfg = get_prev_config(scenario = "freeze"),
    emit_prev_diag_console = FALSE
  ), extra_args))

  runs <- list()
  .pack_one <- function(res, scenario_name) {
    attr(res, "cause_id") <- cfg_row$cause_id[[1]]
    attr(res, "label") <- cfg_row$label[[1]]
    attr(res, "scenario") <- scenario_name
    params_tbl <- pack_params(res, cfg_row$cause_id[[1]], cfg_row$label[[1]]) |> dplyr::mutate(scenario = scenario_name, .before = 1)
    proj_tbl <- pack_proj(res, cfg_row$cause_id[[1]], cfg_row$label[[1]]) |> dplyr::mutate(scenario = scenario_name, .before = 1)
    horizon_tbl <- pack_horizon(res, cfg_row$cause_id[[1]], cfg_row$label[[1]]) |> dplyr::mutate(scenario = scenario_name, .before = 1)
    list(scenario = scenario_name, prev_cfg = get_prev_config(scenario = scenario_name),
         res = res, params_tbl = params_tbl, proj_tbl = proj_tbl, horizon_tbl = horizon_tbl)
  }

  runs[["freeze"]] <- .pack_one(res_freeze, "freeze")
  if (isTRUE(save_raw_rds)) saveRDS(res_freeze, file.path(OUT_RAW_URUGUAY, "res_lung_freeze.rds"))

  for (scn in setdiff(scenarios, "freeze")) {
    message(">>> Uruguay empirical run: rebuilding scenario ", scn)
    out_rebuild <- do.call(.rebuild_scenario_freeze_benchmark, c(list(
      res_base = res_freeze,
      inputs = inputs,
      cfg_row = cfg_row,
      prev_cfg_scen = get_prev_config(scenario = scn),
      overwrite_main = TRUE
    ), extra_args))
    res_freeze <- out_rebuild$res_base
    res_scen <- out_rebuild$res_scen
    runs[[scn]] <- .pack_one(res_scen, scn)
    if (isTRUE(save_raw_rds)) saveRDS(res_scen, file.path(OUT_RAW_URUGUAY, sprintf("res_lung_%s.rds", scn)))
  }
  runs[scenarios]
}

plot_uruguay_projection_panel <- function(run_list, show_ci = FALSE) {
  proj <- combined_projection_tables(run_list) |>
    dplyr::filter(series == main_series_for_metric(metric), period >= 2022) |>
    prepare_segmented_lines()
  obs <- observed_projection_tables(run_list) |>
    dplyr::filter(period >= 1998, period <= 2022)

  g <- ggplot2::ggplot() +
    ggplot2::geom_line(data = obs, ggplot2::aes(period, obs), color = "black", linewidth = 0.65) +
    ggplot2::geom_vline(xintercept = 2022.5, color = "grey45", linewidth = 0.35)
  if (isTRUE(show_ci)) {
    g <- g + ggplot2::geom_ribbon(
      data = proj,
      ggplot2::aes(period, ymin = lwr, ymax = upr, fill = scenario_label,
                   group = interaction(scenario_label, zone_group)),
      alpha = 0.08, color = NA, show.legend = FALSE
    )
  }
  g <- g +
    ggplot2::geom_line(
      data = proj,
      ggplot2::aes(period, mean, color = scenario_label, linetype = projection_zone,
                   group = interaction(scenario_label, zone_group)),
      linewidth = 0.78
    ) +
    ggplot2::facet_grid(metric_label ~ sex_label, scales = "free_y") +
    ggplot2::scale_color_manual(values = URUGUAY_COLORS, name = "Scenario") +
    ggplot2::scale_linetype_manual(
      values = URUGUAY_ZONE_LINETYPES,
      breaks = c("credible", "caution", "risky"),
      labels = c("Credible", "Caution", "Risky"),
      name = "Horizon"
    ) +
    ggplot2::scale_y_continuous(limits = c(0, NA), expand = ggplot2::expansion(mult = c(0, 0.04))) +
    ggplot2::labs(x = "Year", y = "Annual cases/deaths") +
    theme_paper_main(base_size = 10.8) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.box = "vertical",
      strip.background = ggplot2::element_rect(fill = "gray95"),
      panel.spacing = grid::unit(0.7, "lines")
    )
  if (isTRUE(show_ci)) {
    g <- g + ggplot2::scale_fill_manual(values = URUGUAY_COLORS, name = "Scenario")
  }
  g
}

extract_smoking_exposure <- function(run_list) {
  rows <- list()
  for (scn in names(run_list)) {
    for (sx in c("M", "F")) {
      rs <- if (sx == "M") run_list[[scn]]$res$resM else run_list[[scn]]$res$resF
      dat <- tryCatch(tibble::as_tibble(rs$inc_fit$rates_all_full), error = function(e) tibble::tibble())
      if (!nrow(dat)) next
      if (!"E" %in% names(dat)) dat$E <- NA_real_
      if (!"exposure" %in% names(dat)) dat$exposure <- NA_real_
      rows[[paste(scn, sx)]] <- dat |>
        dplyr::filter(period >= 2001, period <= PROJ_TO, age >= AGE_I_MIN, age <= AGE_I_MAX) |>
        dplyr::mutate(E = dplyr::coalesce(as.numeric(E), as.numeric(exposure), 1)) |>
        dplyr::group_by(period) |>
        dplyr::summarise(
          current_smoking = stats::weighted.mean(as.numeric(p_cur), E, na.rm = TRUE) * 100,
          effective_exposure = stats::weighted.mean(as.numeric(q_eff), E, na.rm = TRUE) * 100,
          .groups = "drop"
        ) |>
        tidyr::pivot_longer(c(current_smoking, effective_exposure), names_to = "measure", values_to = "value") |>
        dplyr::mutate(scenario = scn, sex = sx)
    }
  }
  dplyr::bind_rows(rows) |>
    dplyr::mutate(
      scenario_label = factor(unname(URUGUAY_SCEN_LABELS[scenario]), levels = unname(URUGUAY_SCEN_LABELS[URUGUAY_SCENARIOS])),
      sex_label = factor(sex_public(sex), levels = c("Male", "Female")),
      measure_label = factor(dplyr::recode(measure,
                                           current_smoking = "Current smoking prevalence",
                                           effective_exposure = "Effective smoking exposure"),
                             levels = c("Current smoking prevalence", "Effective smoking exposure"))
    )
}

plot_uruguay_smoking_exposure <- function(run_list) {
  df <- extract_smoking_exposure(run_list)
  ggplot2::ggplot(df, ggplot2::aes(period, value, color = scenario_label)) +
    ggplot2::geom_vline(xintercept = 2022.5, color = "grey45", linewidth = 0.35) +
    ggplot2::geom_line(linewidth = 0.78, na.rm = TRUE) +
    ggplot2::facet_grid(measure_label ~ sex_label, scales = "free_y") +
    ggplot2::scale_color_manual(values = URUGUAY_COLORS, name = "Scenario") +
    ggplot2::labs(x = "Year", y = "Percent") +
    theme_paper_main(base_size = 10.8) +
    ggplot2::theme(
      legend.position = "bottom",
      strip.background = ggplot2::element_rect(fill = "gray95"),
      panel.spacing = grid::unit(0.7, "lines")
    )
}

build_uruguay_data_overview <- function(inputs = NULL) {
  if (is.null(inputs)) inputs <- build_inputs_real_cause(run_cfg$causes_tbl[1, ])
  prev_raw <- haven::read_dta(inputs$prev_path)
  prev_names <- names(prev_raw)
  year_col <- intersect(c("period", "year", "anio", "ano", intToUtf8(c(0x61, 0xf1, 0x6f))), prev_names)[1]
  if (is.na(year_col)) stop("Could not identify year column in prevalence microdata.")
  micro <- prev_raw |>
    dplyr::transmute(
      period = as.integer(.data[[year_col]]),
      age = as.integer(edad),
      sex = ifelse(as.numeric(mujer) == 1, "F", "M"),
      fuma = as.numeric(fuma),
      w = as.numeric(expansor),
      inst_ok = as.integer(d_act) + as.integer(d_12m) + as.integer(d_30d)
    ) |>
    dplyr::filter(inst_ok == 1, period >= 1998, period <= 2022, age >= AGE_P_MIN, age <= AGE_P_MAX,
                  sex %in% c("M", "F"), is.finite(w), w > 0) |>
    dplyr::group_by(period, sex) |>
    dplyr::summarise(value = stats::weighted.mean(fuma, w, na.rm = TRUE) * 100,
                     .groups = "drop") |>
    dplyr::mutate(measure = "Current smoking prevalence")

  inc <- inputs$inc_hist_tbl |>
    dplyr::filter(period >= 1998, period <= 2022, age >= AGE_I_MIN, age <= AGE_I_MAX) |>
    dplyr::group_by(period, sex) |>
    dplyr::summarise(cases = sum(cases, na.rm = TRUE),
                     exposure = sum(exposure, na.rm = TRUE), .groups = "drop") |>
    dplyr::mutate(value = 1e5 * cases / exposure, measure = "Lung-cancer incidence rate") |>
    dplyr::select(period, sex, value, measure)

  mort <- inputs$mort_hist_tbl |>
    dplyr::filter(period >= 1998, period <= 2022, age >= AGE_M_MIN, age <= AGE_M_MAX) |>
    dplyr::group_by(period, sex) |>
    dplyr::summarise(deaths = sum(deaths, na.rm = TRUE),
                     exposure = sum(exposure, na.rm = TRUE), .groups = "drop") |>
    dplyr::mutate(value = 1e5 * deaths / exposure, measure = "Lung-cancer mortality rate") |>
    dplyr::select(period, sex, value, measure)

  dplyr::bind_rows(micro, inc, mort) |>
    dplyr::mutate(
      sex_label = factor(sex_public(sex), levels = c("Male", "Female")),
      measure_label = factor(measure, levels = c("Current smoking prevalence",
                                                 "Lung-cancer incidence rate",
                                                 "Lung-cancer mortality rate"))
    )
}

plot_uruguay_data_overview <- function(inputs = NULL) {
  df <- build_uruguay_data_overview(inputs)
  ggplot2::ggplot(df, ggplot2::aes(period, value, color = sex_label)) +
    ggplot2::geom_line(linewidth = 0.78, na.rm = TRUE) +
    ggplot2::geom_point(size = 1.25, na.rm = TRUE) +
    ggplot2::facet_wrap(~ measure_label, scales = "free_y", ncol = 1) +
    ggplot2::scale_color_manual(values = c(Male = "#1B4F72", Female = "#A93226"), name = "Sex") +
    ggplot2::labs(x = "Year", y = NULL) +
    theme_paper_main(base_size = 10.6) +
    ggplot2::theme(
      legend.position = "bottom",
      strip.background = ggplot2::element_rect(fill = "gray95")
    )
}

plot_uruguay_benchmark_comparison <- function(run_list) {
  proj <- combined_projection_tables(run_list) |>
    dplyr::filter(
      period >= 2022,
      (metric == "incidence" & series %in% c("I|P", "I")) |
        (metric == "mortality" & series %in% c("M|I|P", "BAPC benchmark"))
    ) |>
    dplyr::mutate(
      series_label = dplyr::case_when(
        series %in% c("I|P", "M|I|P") ~ "SBAPC",
        series %in% c("I", "BAPC benchmark") ~ "BAPC benchmark",
        TRUE ~ as.character(series)
      ),
      series_label = factor(series_label, levels = c("SBAPC", "BAPC benchmark"))
    ) |>
    prepare_segmented_lines()
  obs <- observed_projection_tables(run_list) |> dplyr::filter(period >= 1998, period <= 2022)
  ggplot2::ggplot() +
    ggplot2::geom_line(data = obs, ggplot2::aes(period, obs), color = "black", linewidth = 0.55) +
    ggplot2::geom_vline(xintercept = 2022.5, color = "grey45", linewidth = 0.35) +
    ggplot2::geom_line(
      data = proj,
      ggplot2::aes(period, mean, color = scenario_label, linetype = series_label,
                   group = interaction(scenario_label, series_label, zone_group)),
      linewidth = 0.68
    ) +
    ggplot2::facet_grid(metric_label ~ sex_label, scales = "free_y") +
    ggplot2::scale_color_manual(values = URUGUAY_COLORS, name = "Scenario") +
    ggplot2::scale_linetype_manual(values = c("SBAPC" = "solid", "BAPC benchmark" = "longdash"),
                                   name = "Series") +
    ggplot2::scale_y_continuous(limits = c(0, NA), expand = ggplot2::expansion(mult = c(0, 0.04))) +
    ggplot2::labs(x = "Year", y = "Annual cases/deaths") +
    theme_paper_main(base_size = 10.4) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.box = "vertical",
      strip.background = ggplot2::element_rect(fill = "gray95"),
      panel.spacing = grid::unit(0.7, "lines")
    )
}

plot_uruguay_horizon_support <- function(run_list) {
  ref <- run_list[["freeze"]] %||% run_list[[1]]
  rows <- list()
  for (sx in c("M", "F")) {
    rs <- if (sx == "M") ref$res$resM else ref$res$resF
    yd <- tryCatch(tibble::as_tibble(rs$diag$projection_horizon_year), error = function(e) tibble::tibble())
    if (!nrow(yd)) next
    rows[[sx]] <- yd |> dplyr::mutate(sex = sx)
  }
  df <- dplyr::bind_rows(rows) |>
    dplyr::filter(period > 2022) |>
    dplyr::select(dplyr::any_of(c("sex", "period", "mean_support_frac", "edge_share", "projection_zone"))) |>
    dplyr::mutate(
      sex_label = factor(sex_public(sex), levels = c("Male", "Female")),
      projection_zone = factor(as.character(projection_zone), levels = c("credible", "caution", "risky", "beyond_max"))
    )
  ggplot2::ggplot(df, ggplot2::aes(period, mean_support_frac, color = projection_zone)) +
    ggplot2::geom_hline(yintercept = c(HORIZON_SUPPORT_FLOOR_CREDIBLE, HORIZON_SUPPORT_FLOOR_CAUTION),
                        color = "grey70", linewidth = 0.35) +
    ggplot2::geom_line(linewidth = 0.78, na.rm = TRUE) +
    ggplot2::facet_wrap(~ sex_label, ncol = 1) +
    ggplot2::scale_color_manual(
      values = c(credible = "#1E8449", caution = "#D4AC0D", risky = "#C0392B", beyond_max = "#7F8C8D"),
      breaks = c("credible", "caution", "risky", "beyond_max"),
      labels = c("Credible", "Caution", "Risky", "Beyond maximum"),
      name = "Horizon"
    ) +
    ggplot2::scale_y_continuous(limits = c(0, 1), labels = scales::percent_format(accuracy = 1)) +
    ggplot2::labs(x = "Year", y = "Mean support fraction") +
    theme_paper_main(base_size = 10.6) +
    ggplot2::theme(legend.position = "bottom")
}

build_total_effects <- function(run_list) {
  proj <- combined_projection_tables(run_list) |>
    dplyr::filter(
      period > 2022,
      sex %in% c("M", "F"),
      (metric == "incidence" & series == "I|P") |
        (metric == "mortality" & series == "M|I|P")
    ) |>
    dplyr::group_by(scenario, scenario_label, metric, metric_label, period) |>
    dplyr::summarise(mean = sum(mean, na.rm = TRUE), .groups = "drop")
  freeze <- proj |>
    dplyr::filter(as.character(scenario) == "freeze") |>
    dplyr::select(metric, period, freeze_mean = mean)
  total_h <- get_total_horizon(run_list)
  proj |>
    dplyr::filter(as.character(scenario) != "freeze") |>
    dplyr::left_join(freeze, by = c("metric", "period")) |>
    dplyr::mutate(
      effect = mean - freeze_mean,
      horizon = horizon_region_from_year(period, total_h)
    ) |>
    dplyr::filter(!is.na(horizon)) |>
    dplyr::group_by(scenario, scenario_label, metric, metric_label, horizon) |>
    dplyr::summarise(
      cumulative_effect = sum(effect, na.rm = TRUE),
      cumulative_freeze = sum(freeze_mean, na.rm = TRUE),
      pct_of_freeze = 100 * cumulative_effect / pmax(cumulative_freeze, 1e-9),
      .groups = "drop"
    ) |>
    dplyr::arrange(factor(as.character(scenario), levels = c("up1pc", "down1pc", "quit")),
                   factor(horizon, levels = c("Credible", "Caution", "Risky")),
                   metric)
}

export_uruguay_cumulative_effects <- function(run_list,
                                              csv_out = file.path(OUT_SECTION5, "tab_uruguay_cumulative_effects.csv"),
                                              tex_out = file.path(OUT_SECTION5, "tab_uruguay_cumulative_effects.tex")) {
  eff <- build_total_effects(run_list)
  readr::write_csv(eff, csv_out)
  wide <- eff |>
    dplyr::select(scenario, scenario_label, horizon, metric, cumulative_effect, pct_of_freeze) |>
    tidyr::pivot_wider(names_from = metric, values_from = c(cumulative_effect, pct_of_freeze)) |>
    dplyr::arrange(factor(as.character(scenario), levels = c("up1pc", "down1pc", "quit")),
                   factor(horizon, levels = c("Credible", "Caution", "Risky")))

  lines <- c(latex_open("llrrrr"),
             "Scenario & Horizon & Inc. effect & Mort. effect & Inc. \\% & Mort. \\% \\\\",
             "\\midrule")
  last_scenario <- NULL
  for (i in seq_len(nrow(wide))) {
    row <- wide[i, ]
    scen <- if (!identical(last_scenario, as.character(row$scenario))) URUGUAY_SCEN_LABELS_TEX[as.character(row$scenario)] else ""
    if (!is.null(last_scenario) && !identical(last_scenario, as.character(row$scenario))) lines <- c(lines, "\\midrule")
    lines <- c(lines, sprintf(
      "%s & %s & %s & %s & %s & %s \\\\",
      scen, row$horizon,
      fmt_int(row$cumulative_effect_incidence),
      fmt_int(row$cumulative_effect_mortality),
      fmt_num(row$pct_of_freeze_incidence, 1),
      fmt_num(row$pct_of_freeze_mortality, 1)
    ))
    last_scenario <- as.character(row$scenario)
  }
  writeLines(c(lines, latex_close()), tex_out, useBytes = TRUE)
  invisible(eff)
}

export_uruguay_horizon_boundaries <- function(run_list,
                                              csv_out = file.path(OUT_APPENDIXD, "tab_uruguay_horizon_boundaries.csv"),
                                              tex_out = file.path(OUT_APPENDIXD, "tab_uruguay_horizon_boundaries.tex")) {
  h <- dplyr::bind_rows(lapply(names(run_list), function(scn) run_list[[scn]]$horizon_tbl |> dplyr::mutate(scenario = scn))) |>
    dplyr::filter(as.character(scenario) == "freeze") |>
    dplyr::mutate(sex_label = sex_public(sex)) |>
    dplyr::select(sex_label, last_hist_year, end_year_credible, end_year_caution, end_year_risky, max_projection_year_endogenous)
  readr::write_csv(h, csv_out)
  lines <- c(latex_open("lrrrrr"),
             "Sex & Last historical & Credible end & Caution end & Risky end & Maximum \\\\",
             "\\midrule")
  for (i in seq_len(nrow(h))) {
    row <- h[i, ]
    lines <- c(lines, sprintf("%s & %s & %s & %s & %s & %s \\\\",
                              row$sex_label,
                              fmt_int(row$last_hist_year),
                              fmt_int(row$end_year_credible),
                              fmt_int(row$end_year_caution),
                              fmt_int(row$end_year_risky),
                              fmt_int(row$max_projection_year_endogenous)))
  }
  writeLines(c(lines, latex_close()), tex_out, useBytes = TRUE)
  invisible(h)
}

export_uruguay_fit_scores <- function(run_list,
                                      csv_out = file.path(OUT_APPENDIXD, "tab_uruguay_fit_scores.csv"),
                                      tex_out = file.path(OUT_APPENDIXD, "tab_uruguay_fit_scores.tex")) {
  fs <- run_list[["freeze"]]$res$fit_scores |>
    dplyr::mutate(
      sex_label = sex_public(sex),
      model_label = dplyr::recode(
        model,
        "Prevalence BAPC" = "Prevalence APC",
        "Incidence benchmark (APC)" = "Incidence BAPC benchmark",
        "Incidence prevalence-informed" = "SBAPC incidence layer",
        "Mortality benchmark (APC)" = "Mortality BAPC benchmark",
        "Mortality anchored on I|P" = "SBAPC mortality layer",
        "Mortality anchored on I only" = "Incidence-anchored mortality layer",
        .default = model
      )
    ) |>
    dplyr::select(sex_label, model_label, model, WAIC, DIC, LCPO, dLCPO, se_dLCPO, BT_LPD, BT_RMSE)
  readr::write_csv(fs, csv_out)
  lines <- c(latex_open("llrrrr"),
             "Sex & Model & WAIC & DIC & LCPO & dLCPO \\\\",
             "\\midrule")
  last_sex <- NULL
  for (i in seq_len(nrow(fs))) {
    row <- fs[i, ]
    sx <- if (!identical(last_sex, as.character(row$sex_label))) row$sex_label else ""
    if (!is.null(last_sex) && !identical(last_sex, as.character(row$sex_label))) lines <- c(lines, "\\midrule")
    lines <- c(lines, sprintf("%s & %s & %s & %s & %s & %s \\\\",
                              sx, row$model_label,
                              fmt_num(row$WAIC, 1), fmt_num(row$DIC, 1),
                              fmt_num(row$LCPO, 1), fmt_num(row$dLCPO, 1)))
    last_sex <- as.character(row$sex_label)
  }
  writeLines(c(lines, latex_close()), tex_out, useBytes = TRUE)
  invisible(fs)
}

export_uruguay_transmission_inputs <- function(run_list,
                                               csv_out = file.path(OUT_SECTION5, "tab_uruguay_transmission_inputs.csv"),
                                               tex_out = file.path(OUT_SECTION5, "tab_uruguay_transmission_inputs.tex")) {
  params <- run_list[["freeze"]]$params_tbl |>
    dplyr::mutate(sex_label = sex_public(sex)) |>
    dplyr::select(sex_label, rr_inc, mort_kernel_max_lag, mort_kernel_total_mass)
  readr::write_csv(params, csv_out)
  lines <- c(latex_open("lrrr"),
             "Sex & Incidence RR & Max lag & Mortality mass \\\\",
             "\\midrule")
  for (i in seq_len(nrow(params))) {
    row <- params[i, ]
    lines <- c(lines, sprintf("%s & %s & %s & %s \\\\",
                              row$sex_label, fmt_num(row$rr_inc, 2),
                              fmt_int(row$mort_kernel_max_lag),
                              fmt_num(row$mort_kernel_total_mass, 3)))
  }
  writeLines(c(lines, latex_close()), tex_out, useBytes = TRUE)
  invisible(params)
}

write_uruguay_notes_and_inventories <- function() {
  section5_fig <- c(
    "# Figure Titles and Notes: Section 5",
    "",
    "## fig_uruguay_smoking_exposure",
    "Files: fig_uruguay_smoking_exposure.svg, fig_uruguay_smoking_exposure.pdf",
    "Title: Estimated smoking exposure under prevalence scenarios",
    "Note: Panels report model-implied current smoking prevalence and effective smoking exposure for the lung-cancer incidence layer, aggregated over the incidence age range using exposure weights. Source: Own elaboration.",
    "",
    "## fig_uruguay_projection_panel",
    "Files: fig_uruguay_projection_panel.svg, fig_uruguay_projection_panel.pdf",
    "Title: Lung-cancer incidence and mortality projections in Uruguay",
    "Note: Black lines show historical preprocessed inputs. Colored lines show SBAPC projections under smoking-prevalence scenarios; line type indicates the endogenous horizon region. Source: Own elaboration."
  )
  writeLines(section5_fig, file.path(OUT_SECTION5, "figure_titles_notes.md"), useBytes = TRUE)

  section5_tab <- c(
    "# Table Titles and Notes: Section 5",
    "",
    "## tab_uruguay_cumulative_effects",
    "Title: Cumulative scenario effects in the Uruguay application",
    "Note: Effects are cumulative projected incident cases or deaths relative to the frozen-prevalence baseline, aggregated across sexes within each endogenous horizon region. Source: Own elaboration.",
    "",
    "## tab_uruguay_transmission_inputs",
    "Title: Lung-cancer transmission inputs used in the Uruguay application",
    "Note: Incidence relative risks and post-diagnosis mortality kernel summaries are fixed external inputs used by the empirical pipeline. Source: Own elaboration."
  )
  writeLines(section5_tab, file.path(OUT_SECTION5, "table_titles_notes.md"), useBytes = TRUE)

  appendix_fig <- c(
    "# Figure Titles and Notes: Appendix D",
    "",
    "## fig_uruguay_data_overview",
    "Files: fig_uruguay_data_overview.svg, fig_uruguay_data_overview.pdf",
    "Title: Historical empirical inputs for Uruguay",
    "Note: The figure summarizes weighted smoking prevalence and preprocessed lung-cancer incidence and mortality inputs before model-based projection. Source: Own elaboration.",
    "",
    "## fig_uruguay_benchmark_comparison",
    "Files: fig_uruguay_benchmark_comparison.svg, fig_uruguay_benchmark_comparison.pdf",
    "Title: SBAPC and BAPC benchmark projections",
    "Note: The display contrasts scenario-responsive SBAPC projections with the scenario-blind BAPC benchmark. Source: Own elaboration.",
    "",
    "## fig_uruguay_horizon_support",
    "Files: fig_uruguay_horizon_support.svg, fig_uruguay_horizon_support.pdf",
    "Title: Endogenous support diagnostics for Uruguay",
    "Note: Lines report the exposure-weighted mean support fraction by projection year and sex. Source: Own elaboration."
  )
  writeLines(appendix_fig, file.path(OUT_APPENDIXD, "figure_titles_notes.md"), useBytes = TRUE)

  appendix_tab <- c(
    "# Table Titles and Notes: Appendix D",
    "",
    "## tab_uruguay_horizon_boundaries",
    "Title: Endogenous horizon boundaries in the Uruguay application",
    "Note: Boundaries are computed from the Lexis-support diagnostics used to classify projection years. Source: Own elaboration.",
    "",
    "## tab_uruguay_fit_scores",
    "Title: Historical fit diagnostics for the Uruguay freeze baseline",
    "Note: Fit statistics are reported for the freeze-baseline empirical run and should be interpreted as diagnostics, not as the main validation target. Source: Own elaboration."
  )
  writeLines(appendix_tab, file.path(OUT_APPENDIXD, "table_titles_notes.md"), useBytes = TRUE)

  section5_inv <- c(
    "# Section 5 Float Inventory",
    "",
    "| Filename | Document | Priority | Aggregation | SVG+PDF | Source note | Purpose |",
    "|---|---|---|---|---|---|---|",
    "| `fig_uruguay_smoking_exposure` | Main text | Useful | Sex-specific annual aggregates | Yes | Yes | Shows the upstream smoking signal and effective exposure entering the incidence layer. |",
    "| `fig_uruguay_projection_panel` | Main text | Essential | Sex-specific annual counts | Yes | Yes | Main empirical projection figure for incidence and mortality under the four smoking scenarios. |",
    "| `tab_uruguay_cumulative_effects` | Main text | Useful | Both-sex cumulative effects | Not applicable | Yes | Compact numerical summary of scenario contrasts relative to frozen prevalence. |",
    "| `tab_uruguay_transmission_inputs` | Main text or Appendix D | Optional | Sex-specific inputs | Not applicable | Yes | Documents the fixed external inputs used by the lung-cancer transmission block. |"
  )
  writeLines(section5_inv, file.path(OUT_SECTION5, "section5_float_inventory.md"), useBytes = TRUE)

  appendix_inv <- c(
    "# Appendix D Float Inventory",
    "",
    "| Filename | Document | Priority | Aggregation | SVG+PDF | Source note | Purpose |",
    "|---|---|---|---|---|---|---|",
    "| `fig_uruguay_data_overview` | Appendix D | Essential | Historical sex-specific annual aggregates | Yes | Yes | Describes the empirical inputs entering the Uruguay application. |",
    "| `fig_uruguay_benchmark_comparison` | Appendix D | Useful | Sex-specific annual counts | Yes | Yes | Shows how the sequential scenario-responsive projection differs from the BAPC benchmark. |",
    "| `fig_uruguay_horizon_support` | Appendix D | Useful | Sex-specific horizon diagnostics | Yes | Yes | Displays the support deterioration underlying the horizon categories. |",
    "| `tab_uruguay_horizon_boundaries` | Appendix D | Useful | Sex-specific and total frontier rows | Not applicable | Yes | Reports projection horizon boundary years. |",
    "| `tab_uruguay_fit_scores` | Appendix D | Optional | Freeze-baseline fit statistics | Not applicable | Yes | Historical fit/backtesting diagnostic. |"
  )
  writeLines(appendix_inv, file.path(OUT_APPENDIXD, "appendixD_float_inventory.md"), useBytes = TRUE)

  invisible(TRUE)
}

export_uruguay_products <- function(run_list, inputs = NULL) {
  if (is.null(inputs)) inputs <- build_inputs_real_cause(run_cfg$causes_tbl[1, ])

  readr::write_csv(combined_projection_tables(run_list), file.path(URUGUAY_OUT_BASE, "uruguay_projection_long.csv"))
  readr::write_csv(observed_projection_tables(run_list), file.path(URUGUAY_OUT_BASE, "uruguay_observed_long.csv"))
  readr::write_csv(extract_smoking_exposure(run_list), file.path(URUGUAY_OUT_BASE, "uruguay_smoking_exposure_long.csv"))

  save_uruguay_plot(plot_uruguay_smoking_exposure(run_list),
                    file.path(OUT_SECTION5, "fig_uruguay_smoking_exposure"),
                    width = 8.9, height = 6.0)
  save_uruguay_plot(plot_uruguay_projection_panel(run_list, show_ci = FALSE),
                    file.path(OUT_SECTION5, "fig_uruguay_projection_panel"),
                    width = 9.2, height = 6.8)

  export_uruguay_cumulative_effects(run_list)
  export_uruguay_transmission_inputs(run_list)

  save_uruguay_plot(plot_uruguay_data_overview(inputs),
                    file.path(OUT_APPENDIXD, "fig_uruguay_data_overview"),
                    width = 7.6, height = 7.6)
  save_uruguay_plot(plot_uruguay_benchmark_comparison(run_list),
                    file.path(OUT_APPENDIXD, "fig_uruguay_benchmark_comparison"),
                    width = 9.4, height = 6.9)
  save_uruguay_plot(plot_uruguay_horizon_support(run_list),
                    file.path(OUT_APPENDIXD, "fig_uruguay_horizon_support"),
                    width = 7.4, height = 5.4)

  export_uruguay_horizon_boundaries(run_list)
  export_uruguay_fit_scores(run_list)
  write_uruguay_notes_and_inventories()

  invisible(list(section5 = OUT_SECTION5, appendixD = OUT_APPENDIXD))
}

replicate_uruguay_empirical <- function(scenarios = URUGUAY_SCENARIOS,
                                        save_raw_rds = FALSE) {
  inputs <- build_inputs_real_cause(run_cfg$causes_tbl[1, ])
  runs <- run_uruguay_lung_rebuilt_all_scenarios(run_cfg = run_cfg, scenarios = scenarios, save_raw_rds = save_raw_rds)
  export_uruguay_products(runs, inputs = inputs)
  invisible(list(runs = runs, out_base = URUGUAY_OUT_BASE))
}

if (sys.nframe() == 0L) {
  replicate_uruguay_empirical()
}
