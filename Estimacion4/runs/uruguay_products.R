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
  "Frozen at 2022" = "#FDAE61",
  "Up 1% p.a." = "#D73027",
  "Down 1% p.a." = "#00796B",
  "Quit" = "#512DA8"
)
URUGUAY_ZONE_LINETYPES <- c(
  credible = "solid",
  caution = "longdash",
  risky = "dotted",
  beyond_max = "dotdash",
  historical = "solid"
)
URUGUAY_SITE_COLORS <- c(
  "Lung" = "#4E79A7",
  "Stomach" = "#F28E2B",
  "Pancreas" = "#E15759",
  "Bladder" = "#76B7B2",
  "Oral cavity and pharynx" = "#59A14F",
  "Esophagus" = "#EDC948",
  "Kidney" = "#B07AA1",
  "Larynx" = "#FF9DA7",
  "Cervix" = "#9C755F"
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

site_public <- function(cause_id, fallback = NULL) {
  cause_id <- as.character(cause_id)
  fallback_chr <- if (is.null(fallback)) cause_id else rep_len(as.character(fallback), length(cause_id))
  vapply(seq_along(cause_id), function(i) get_cause_label_en(cause_id[[i]], fallback = fallback_chr[[i]] %||% cause_id[[i]]), character(1))
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
  out <- dplyr::bind_rows(rows)
  horizon_year <- tryCatch(
    projection_common_max_year_from_res_both(run_out$res, policy = "endogenous_max"),
    error = function(e) NA_integer_
  )
  if (is.finite(horizon_year)) {
    out <- clip_to_year(out, max_year = horizon_year, year_var = "period")
  }

  horizon_zone <- dplyr::bind_rows(
    tryCatch(run_out$res$resM$diag$projection_horizon_year, error = function(e) NULL),
    tryCatch(run_out$res$resF$diag$projection_horizon_year, error = function(e) NULL)
  )
  if (is.data.frame(horizon_zone) && nrow(horizon_zone)) {
    horizon_zone <- horizon_zone |>
      dplyr::transmute(
        sex = as.character(sex),
        period = suppressWarnings(as.integer(period)),
        projection_zone = as.character(projection_zone)
      )
    out <- out |> dplyr::left_join(horizon_zone, by = c("sex", "period"))
  } else {
    out$projection_zone <- NA_character_
  }

  hist_cutoff <- tryCatch(suppressWarnings(as.integer(run_out$res$combined$last_hist_year)[1]), error = function(e) NA_integer_)
  if (is.finite(hist_cutoff)) {
    out <- out |>
      dplyr::mutate(
        projection_zone = dplyr::case_when(
          period <= hist_cutoff ~ "historical",
          TRUE ~ dplyr::coalesce(projection_zone, "beyond_max")
        ),
        projection_zone = factor(projection_zone, levels = c("historical", "credible", "caution", "risky", "beyond_max"))
      )
  }
  out
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

get_common_endogenous_year <- function(run_list) {
  h <- dplyr::bind_rows(lapply(names(run_list), function(scn) {
    run_list[[scn]]$horizon_tbl |> dplyr::mutate(scenario = scn)
  }))
  vals <- h |>
    dplyr::filter(as.character(sex) == "T") |>
    dplyr::pull(max_projection_year_endogenous)
  vals <- suppressWarnings(as.integer(vals))
  vals <- vals[is.finite(vals)]
  if (!length(vals)) return(NA_integer_)
  min(vals, na.rm = TRUE)
}

selected_endogenous_years <- function(run_list, candidates = c(2022L, 2030L, 2040L)) {
  max_year <- get_common_endogenous_year(run_list)
  yrs <- unique(as.integer(c(candidates, max_year)))
  yrs[is.finite(yrs) & (is.na(max_year) | yrs <= max_year)]
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

build_uruguay_projection_panel_data <- function(run_list) {
  combined_projection_tables(run_list) |>
    dplyr::filter(
      period >= 2022,
      (metric == "incidence" & series == "I|P") |
        (metric == "mortality" & series == "M|I|P")
    )
}

plot_uruguay_projection_panel <- function(run_list, show_ci = FALSE) {
  proj <- build_uruguay_projection_panel_data(run_list) |>
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

build_uruguay_lung_mortality_by_sex_data <- function(run_list) {
  proj <- combined_projection_tables(run_list) |>
    dplyr::filter(period >= 2022, sex %in% c("M", "F"), metric == "mortality", series == "M|I|P") |>
    dplyr::transmute(
      source = "Projected",
      scenario = as.character(scenario),
      scenario_label,
      sex = as.character(sex),
      sex_label,
      metric = "mortality",
      series = "SBAPC",
      period = as.integer(period),
      deaths = as.numeric(mean),
      lwr = as.numeric(lwr),
      upr = as.numeric(upr),
      projection_zone
    )
  obs <- observed_projection_tables(run_list) |>
    dplyr::filter(period >= 1998, period <= 2022, sex %in% c("M", "F"), metric == "mortality") |>
    dplyr::transmute(
      source = "Observed",
      scenario = NA_character_,
      scenario_label = factor(NA_character_, levels = levels(proj$scenario_label)),
      sex = as.character(sex),
      sex_label,
      metric = "mortality",
      series = "Observed",
      period = as.integer(period),
      deaths = as.numeric(obs),
      lwr = NA_real_,
      upr = NA_real_,
      projection_zone = factor("historical", levels = levels(proj$projection_zone))
    )
  dplyr::bind_rows(obs, proj)
}

plot_uruguay_lung_mortality_by_sex <- function(run_list) {
  dat <- build_uruguay_lung_mortality_by_sex_data(run_list)
  proj <- dat |>
    dplyr::filter(source == "Projected") |>
    prepare_segmented_lines()
  obs <- dat |> dplyr::filter(source == "Observed")

  ggplot2::ggplot() +
    ggplot2::geom_line(data = obs, ggplot2::aes(period, deaths), color = "black", linewidth = 0.65) +
    ggplot2::geom_vline(xintercept = 2022.5, color = "grey45", linewidth = 0.35) +
    ggplot2::geom_ribbon(
      data = proj,
      ggplot2::aes(period, ymin = lwr, ymax = upr, fill = scenario_label,
                   group = interaction(scenario_label, zone_group)),
      alpha = 0.10,
      color = NA,
      show.legend = FALSE
    ) +
    ggplot2::geom_line(
      data = proj,
      ggplot2::aes(period, deaths, color = scenario_label, linetype = projection_zone,
                   group = interaction(scenario_label, zone_group)),
      linewidth = 0.82
    ) +
    ggplot2::facet_wrap(~ sex_label, ncol = 2, scales = "free_y") +
    ggplot2::scale_color_manual(values = URUGUAY_COLORS, name = "Scenario") +
    ggplot2::scale_fill_manual(values = URUGUAY_COLORS, guide = "none") +
    ggplot2::scale_linetype_manual(
      values = URUGUAY_ZONE_LINETYPES,
      breaks = c("credible", "caution", "risky"),
      labels = c("Credible", "Caution", "Risky"),
      name = "Horizon"
    ) +
    ggplot2::scale_y_continuous(limits = c(0, NA), expand = ggplot2::expansion(mult = c(0, 0.04))) +
    ggplot2::labs(x = "Year", y = "Annual lung-cancer deaths") +
    theme_paper_main(base_size = 11.0) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.box = "vertical",
      strip.background = ggplot2::element_rect(fill = "gray95")
    )
}

export_uruguay_lung_mortality_selected_years <- function(
    run_list,
    csv_out = file.path(OUT_SECTION5, "tab_uruguay_lung_mortality_selected_years.csv"),
    tex_out = file.path(OUT_SECTION5, "tab_uruguay_lung_mortality_selected_years.tex")) {
  years <- selected_endogenous_years(run_list)
  proj <- combined_projection_tables(run_list) |>
    dplyr::filter(sex %in% c("M", "F"), metric == "mortality", series == "M|I|P", period %in% years) |>
    dplyr::transmute(
      sex = as.character(sex),
      sex_label = sex_public(sex),
      scenario = as.character(scenario),
      scenario_label = unname(URUGUAY_SCEN_LABELS[scenario]),
      period = as.integer(period),
      deaths = as.numeric(mean)
    )
  obs_2022 <- observed_projection_tables(run_list) |>
    dplyr::filter(sex %in% c("M", "F"), metric == "mortality", period == 2022L) |>
    dplyr::transmute(sex = as.character(sex), period = as.integer(period), obs_deaths = as.numeric(obs))
  tab <- proj |>
    dplyr::left_join(obs_2022, by = c("sex", "period")) |>
    dplyr::mutate(deaths = dplyr::if_else(period == 2022L & is.finite(obs_deaths), obs_deaths, deaths)) |>
    dplyr::select(sex, sex_label, scenario, scenario_label, period, deaths) |>
    dplyr::arrange(factor(sex, levels = c("M", "F")), factor(scenario, levels = URUGUAY_SCENARIOS), period)
  readr::write_csv(tab, csv_out)

  wide <- tab |>
    tidyr::pivot_wider(names_from = period, values_from = deaths, names_prefix = "y_") |>
    dplyr::arrange(factor(sex, levels = c("M", "F")), factor(scenario, levels = URUGUAY_SCENARIOS))
  year_cols <- paste0("y_", years)
  lines <- c(
    latex_open(paste0("ll", paste(rep("r", length(years)), collapse = ""))),
    paste(c("Sex", "Scenario", as.character(years)), collapse = " & ") |> paste0(" \\\\"),
    "\\midrule"
  )
  last_sex <- NULL
  for (i in seq_len(nrow(wide))) {
    row <- wide[i, ]
    sx <- as.character(row$sex)
    if (!is.null(last_sex) && !identical(last_sex, sx)) lines <- c(lines, "\\midrule")
    sex_cell <- if (!identical(last_sex, sx)) row$sex_label else ""
    vals <- vapply(year_cols, function(cc) fmt_int(row[[cc]]), character(1))
    lines <- c(lines, paste(c(sex_cell, row$scenario_label, vals), collapse = " & ") |> paste0(" \\\\"))
    last_sex <- sx
  }
  writeLines(c(lines, latex_close()), tex_out, useBytes = TRUE)
  invisible(tab)
}

build_uruguay_mortality_effects <- function(run_list, segmented = FALSE) {
  proj <- combined_projection_tables(run_list) |>
    dplyr::filter(period > 2022, sex %in% c("M", "F"), metric == "mortality", series == "M|I|P") |>
    dplyr::group_by(scenario, scenario_label, period) |>
    dplyr::summarise(mean = sum(mean, na.rm = TRUE), .groups = "drop")
  freeze <- proj |>
    dplyr::filter(as.character(scenario) == "freeze") |>
    dplyr::select(period, freeze_mean = mean)
  total_h <- get_total_horizon(run_list)
  out <- proj |>
    dplyr::filter(as.character(scenario) != "freeze") |>
    dplyr::left_join(freeze, by = "period") |>
    dplyr::mutate(
      effect = mean - freeze_mean,
      horizon = horizon_region_from_year(period, total_h),
      projection_zone = factor(tolower(horizon), levels = c("credible", "caution", "risky", "beyond_max")),
      metric = "mortality_effect",
      sex = "T",
      series = "SBAPC"
    ) |>
    dplyr::filter(!is.na(horizon))
  if (isTRUE(segmented)) out <- prepare_segmented_lines(out)
  out
}

plot_uruguay_mortality_effects <- function(run_list) {
  df <- build_uruguay_mortality_effects(run_list, segmented = TRUE)
  ggplot2::ggplot(df, ggplot2::aes(period, effect, color = scenario_label, linetype = projection_zone,
                                   group = interaction(scenario_label, zone_group))) +
    ggplot2::geom_hline(yintercept = 0, color = "grey50", linewidth = 0.35) +
    ggplot2::geom_line(linewidth = 0.82, na.rm = TRUE) +
    ggplot2::scale_color_manual(values = URUGUAY_COLORS, name = "Scenario") +
    ggplot2::scale_linetype_manual(
      values = URUGUAY_ZONE_LINETYPES,
      breaks = c("credible", "caution", "risky"),
      labels = c("Credible", "Caution", "Risky"),
      name = "Horizon"
    ) +
    ggplot2::labs(x = "Year", y = "Annual deaths relative to frozen prevalence") +
    theme_paper_main(base_size = 11.2) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.box = "vertical"
    )
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

build_uruguay_benchmark_comparison_data <- function(run_list) {
  combined_projection_tables(run_list) |>
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
    )
}

plot_uruguay_benchmark_comparison <- function(run_list) {
  proj <- build_uruguay_benchmark_comparison_data(run_list) |>
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

build_uruguay_horizon_support_data <- function(run_list) {
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
  df
}

plot_uruguay_horizon_support <- function(run_list) {
  df <- build_uruguay_horizon_support_data(run_list)
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

build_uruguay_lung_mortality_uncertainty_data <- function(run_list, scenarios = c("freeze", "quit")) {
  combined_projection_tables(run_list) |>
    dplyr::filter(
      period >= 2022,
      sex %in% c("M", "F"),
      metric == "mortality",
      series == "M|I|P",
      as.character(scenario) %in% scenarios
    ) |>
    dplyr::transmute(
      scenario = as.character(scenario),
      scenario_label,
      sex = as.character(sex),
      sex_label,
      metric = "mortality",
      series = "SBAPC",
      period = as.integer(period),
      mean = as.numeric(mean),
      lwr = as.numeric(lwr),
      upr = as.numeric(upr),
      projection_zone
    )
}

plot_uruguay_lung_mortality_uncertainty <- function(run_list) {
  df <- build_uruguay_lung_mortality_uncertainty_data(run_list) |>
    prepare_segmented_lines()
  ggplot2::ggplot(df, ggplot2::aes(period, mean, color = scenario_label, fill = scenario_label,
                                   group = interaction(scenario_label, zone_group))) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lwr, ymax = upr), alpha = 0.12, color = NA) +
    ggplot2::geom_line(ggplot2::aes(linetype = projection_zone), linewidth = 0.78, na.rm = TRUE) +
    ggplot2::facet_wrap(~ sex_label, ncol = 2, scales = "free_y") +
    ggplot2::scale_color_manual(values = URUGUAY_COLORS, name = "Scenario") +
    ggplot2::scale_fill_manual(values = URUGUAY_COLORS, name = "Scenario") +
    ggplot2::scale_linetype_manual(
      values = URUGUAY_ZONE_LINETYPES,
      breaks = c("credible", "caution", "risky"),
      labels = c("Credible", "Caution", "Risky"),
      name = "Horizon"
    ) +
    ggplot2::scale_y_continuous(limits = c(0, NA), expand = ggplot2::expansion(mult = c(0, 0.04))) +
    ggplot2::labs(x = "Year", y = "Annual lung-cancer deaths") +
    theme_paper_main(base_size = 10.7) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.box = "vertical",
      strip.background = ggplot2::element_rect(fill = "gray95")
    )
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
                                              csv_out = file.path(OUT_APPENDIXD, "tab_uruguay_cumulative_effects.csv"),
                                              tex_out = file.path(OUT_APPENDIXD, "tab_uruguay_cumulative_effects.tex")) {
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
    dplyr::select(sex_label, end_year_credible, end_year_caution, end_year_risky, max_projection_year_endogenous)
  readr::write_csv(h, csv_out)
  lines <- c(latex_open("lrrrr"),
             "Sex & Credible end & Caution end & Risky end & Maximum endogenous year \\\\",
             "\\midrule")
  for (i in seq_len(nrow(h))) {
    row <- h[i, ]
    lines <- c(lines, sprintf("%s & %s & %s & %s & %s \\\\",
                              row$sex_label,
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

build_uruguay_fit_diagnostics_compact <- function(run_list) {
  fs <- run_list[["freeze"]]$res$fit_scores |>
    dplyr::mutate(
      sex_label = sex_public(sex),
      model_family = dplyr::case_when(
        grepl("Prevalence", model, ignore.case = TRUE) ~ "Prevalence",
        grepl("Incidence", model, ignore.case = TRUE) ~ "Incidence",
        grepl("Mortality", model, ignore.case = TRUE) ~ "Mortality",
        TRUE ~ "Other"
      ),
      estimator = dplyr::case_when(
        grepl("benchmark", model, ignore.case = TRUE) ~ "BAPC benchmark",
        grepl("SBAPC|anchor", model, ignore.case = TRUE) ~ "SBAPC",
        grepl("Prevalence", model, ignore.case = TRUE) ~ "APC",
        TRUE ~ model
      )
    )
  base <- fs |>
    dplyr::filter(estimator == "BAPC benchmark") |>
    dplyr::select(sex, model_family, WAIC_benchmark = WAIC, DIC_benchmark = DIC, LCPO_benchmark = LCPO)
  fs |>
    dplyr::left_join(base, by = c("sex", "model_family")) |>
    dplyr::mutate(
      delta_WAIC = WAIC - WAIC_benchmark,
      delta_DIC = DIC - DIC_benchmark,
      dLCPO = dplyr::coalesce(dLCPO, LCPO - LCPO_benchmark),
      BT_RMSE = suppressWarnings(as.numeric(BT_RMSE))
    ) |>
    dplyr::filter(model_family %in% c("Prevalence", "Incidence", "Mortality")) |>
    dplyr::select(sex, sex_label, model_family, estimator, delta_WAIC, delta_DIC, dLCPO, BT_RMSE) |>
    dplyr::arrange(factor(sex, levels = c("M", "F")),
                   factor(model_family, levels = c("Prevalence", "Incidence", "Mortality")),
                   factor(estimator, levels = c("APC", "BAPC benchmark", "SBAPC")))
}

export_uruguay_fit_diagnostics_compact <- function(
    run_list,
    csv_out = file.path(OUT_APPENDIXD, "tab_uruguay_fit_diagnostics_compact.csv"),
    tex_out = file.path(OUT_APPENDIXD, "tab_uruguay_fit_diagnostics_compact.tex")) {
  fs <- build_uruguay_fit_diagnostics_compact(run_list)
  readr::write_csv(fs, csv_out)
  lines <- c(
    latex_open("lllrrrr"),
    "Sex & Layer & Estimator & $\\Delta$WAIC & $\\Delta$DIC & dLCPO & RMSE \\\\",
    "\\midrule"
  )
  last_sex <- NULL
  for (i in seq_len(nrow(fs))) {
    row <- fs[i, ]
    sx <- as.character(row$sex)
    if (!is.null(last_sex) && !identical(last_sex, sx)) lines <- c(lines, "\\midrule")
    sex_cell <- if (!identical(last_sex, sx)) row$sex_label else ""
    lines <- c(lines, sprintf(
      "%s & %s & %s & %s & %s & %s & %s \\\\",
      sex_cell, row$model_family, row$estimator,
      fmt_num(row$delta_WAIC, 1), fmt_num(row$delta_DIC, 1),
      fmt_num(row$dLCPO, 1), fmt_num(row$BT_RMSE, 1)
    ))
    last_sex <- sx
  }
  writeLines(c(lines, latex_close()), tex_out, useBytes = TRUE)
  invisible(fs)
}

assemble_single_sex_res_both <- function(res_sex, sex_sel) {
  resM <- if (identical(sex_sel, "M")) res_sex else NULL
  resF <- if (identical(sex_sel, "F")) res_sex else NULL
  list(
    resM = resM,
    resF = resF,
    combined = list(
      annual_bapc = res_sex$annual_bapc %||% tibble::tibble(),
      annual_anchor = res_sex$annual_anchor %||% tibble::tibble(),
      annual_anchor_noP = res_sex$annual_anchor_noP %||% tibble::tibble(),
      obs_annual = res_sex$obs_annual %||% tibble::tibble(),
      last_hist_year = tryCatch(res_sex$diag$last_hist_year, error = function(e) NA_integer_),
      projection_horizon_frontier = tryCatch(res_sex$diag$projection_horizon_frontier, error = function(e) tibble::tibble()),
      max_projection_year_endogenous = tryCatch(
        projection_max_year_from_res_sex(res_sex, policy = "endogenous_max"),
        error = function(e) NA_integer_
      )
    ),
    fit_scores = res_sex$fit_scores %||% tibble::tibble()
  )
}

run_pipeline_single_sex_from_inputs <- function(inputs,
                                                cfg_row,
                                                sex_sel = c("M", "F"),
                                                prev_cfg = NULL,
                                                gammaP_method = GAMMAP_METHOD,
                                                trend_type = TREND_TYPE,
                                                emit_prev_diag_console = EMIT_PREV_DIAG_CONSOLE,
                                                emit_prev_diag_write = TRUE,
                                                ...) {
  validate_bapc_inputs(inputs)
  sex_sel <- match.arg(sex_sel)
  cfg_row <- tibble::as_tibble(cfg_row)
  stopifnot(nrow(cfg_row) == 1)
  if (is.null(prev_cfg)) prev_cfg <- get_prev_config(scenario = "freeze")
  res <- run_pipeline_sex(
    sex_sel = sex_sel,
    period_min_m = PERIOD_M_MIN,
    period_max_m = PERIOD_M_MAX,
    period_min_p = if ("PERIOD_P_MIN" %in% names(cfg_row)) cfg_row$PERIOD_P_MIN[[1]] else PERIOD_M_MIN,
    period_max_p = if ("PERIOD_P_MAX" %in% names(cfg_row)) cfg_row$PERIOD_P_MAX[[1]] else PERIOD_M_MAX,
    age_min_m = cfg_row$AGE_M_MIN[[1]],
    age_max_m = cfg_row$AGE_M_MAX[[1]],
    age_min_p = cfg_row$AGE_P_MIN[[1]],
    age_max_p = cfg_row$AGE_P_MAX[[1]],
    age_min_i = cfg_row$AGE_I_MIN[[1]],
    age_max_i = cfg_row$AGE_I_MAX[[1]],
    L_I = L_I_DEFAULT,
    Da_I = DA_I,
    bridge_inc_years = BRIDGE_INC_YEARS,
    prev_cfg = prev_cfg,
    L_I_max_years = if (exists(".extract_scalar", inherits = TRUE)) .extract_scalar(cfg_row$L_I_MAX_YEARS, L_I_MAX_YEARS) else cfg_row$L_I_MAX_YEARS[[1]],
    mort_period_shock_years = if (exists(".extract_intvec", inherits = TRUE)) .extract_intvec(cfg_row$MORT_SHOCK_YEARS) else integer(0),
    mort_downweight_years = if (identical(sex_sel, "F")) {
      if (exists(".extract_intvec", inherits = TRUE)) .extract_intvec(cfg_row$DOWNWEIGHT_F) else integer(0)
    } else {
      integer(0)
    },
    mort_downweight_weight = if (identical(sex_sel, "F")) MORT_DOWNWEIGHT_WEIGHT_F else 1,
    mort_hist_tbl = inputs$mort_hist_tbl,
    pop_all_tbl = inputs$pop_all_tbl,
    inc_hist_tbl = inputs$inc_hist_tbl,
    path_prev_dta = if (!is.null(inputs$prev_path)) inputs$prev_path else PATH_PREV_DTA,
    prev_micro_df = if (is.data.frame(inputs$prev_data)) inputs$prev_data else NULL,
    gammaP_method = gammaP_method,
    trend_type = trend_type,
    use_age_slope = FALSE,
    tech_scenario = MORT_TREND_SCENARIO,
    delta_tech = DELTA_TECH,
    cause_id_override = if ("cause_id" %in% names(cfg_row)) as.character(cfg_row$cause_id[[1]]) else NA_character_,
    ...
  )
  assemble_single_sex_res_both(res, sex_sel = sex_sel)
}

run_uruguay_multisite_rebuilt_all_scenarios <- function(causes_tbl = get("causes", envir = .GlobalEnv),
                                                        scenarios = URUGUAY_SCENARIOS,
                                                        save_raw_rds = FALSE) {
  extra_args <- list(trend_type = "level", gammaP_method = "freeze", sd_theta_IP = 2.0)
  out <- list()
  failures <- list()

  .pack_one <- function(res, cfg_row, scenario_name) {
    cause_id <- as.character(cfg_row$cause_id[[1]])
    label <- as.character(cfg_row$label[[1]])
    attr(res, "cause_id") <- cause_id
    attr(res, "label") <- label
    attr(res, "scenario") <- scenario_name
    params_tbl <- pack_params(res, cause_id, label) |> dplyr::mutate(scenario = scenario_name, .before = 1)
    proj_tbl <- pack_proj(res, cause_id, label) |> dplyr::mutate(scenario = scenario_name, .before = 1)
    horizon_tbl <- pack_horizon(res, cause_id, label) |> dplyr::mutate(scenario = scenario_name, .before = 1)
    list(scenario = scenario_name, cause_id = cause_id, label = label, res = res,
         params_tbl = params_tbl, proj_tbl = proj_tbl, horizon_tbl = horizon_tbl)
  }

  for (i in seq_len(nrow(causes_tbl))) {
    cfg_row <- tibble::as_tibble(causes_tbl[i, ])
    cause_id <- as.character(cfg_row$cause_id[[1]])
    label <- as.character(cfg_row$label[[1]])
    message(">>> Uruguay multisite empirical run: fitting freeze for ", cause_id, " - ", label)
    cause_runs <- tryCatch({
      inputs <- build_inputs_real_cause(cfg_row)
      inc_sexes <- sort(unique(as.character(inputs$inc_hist_tbl$sex)))
      inc_sexes <- inc_sexes[inc_sexes %in% c("M", "F")]
      if (length(inc_sexes) == 1L) {
        message(">>> Uruguay multisite empirical run: using single-sex fit for ", cause_id, " (", inc_sexes[[1]], ")")
        res_freeze <- do.call(run_pipeline_single_sex_from_inputs, c(list(
          inputs = inputs,
          cfg_row = cfg_row,
          sex_sel = inc_sexes[[1]],
          prev_cfg = get_prev_config(scenario = "freeze"),
          emit_prev_diag_console = FALSE
        ), extra_args))
      } else {
        res_freeze <- do.call(run_pipeline_both_from_inputs, c(list(
          inputs = inputs,
          cfg_row = cfg_row,
          prev_cfg = get_prev_config(scenario = "freeze"),
          emit_prev_diag_console = FALSE
        ), extra_args))
      }
      runs <- list(freeze = .pack_one(res_freeze, cfg_row, "freeze"))
      if (isTRUE(save_raw_rds)) {
        saveRDS(res_freeze, file.path(OUT_RAW_URUGUAY, sprintf("res_%s_freeze.rds", cause_id)))
      }
      for (scn in setdiff(scenarios, "freeze")) {
        message(">>> Uruguay multisite empirical run: rebuilding ", cause_id, " scenario ", scn)
        out_rebuild <- do.call(.rebuild_scenario_freeze_benchmark, c(list(
          res_base = res_freeze,
          inputs = inputs,
          cfg_row = cfg_row,
          prev_cfg_scen = get_prev_config(scenario = scn),
          overwrite_main = TRUE
        ), extra_args))
        res_freeze <- out_rebuild$res_base
        res_scen <- out_rebuild$res_scen
        runs[[scn]] <- .pack_one(res_scen, cfg_row, scn)
        if (isTRUE(save_raw_rds)) {
          saveRDS(res_scen, file.path(OUT_RAW_URUGUAY, sprintf("res_%s_%s.rds", cause_id, scn)))
        }
      }
      runs[scenarios]
    }, error = function(e) {
      failures[[length(failures) + 1L]] <<- tibble::tibble(
        cause_id = cause_id,
        label = label,
        error = conditionMessage(e)
      )
      NULL
    })
    if (!is.null(cause_runs)) out[[cause_id]] <- cause_runs
  }

  list(runs = out, failures = dplyr::bind_rows(failures))
}

multisite_projection_long <- function(multisite_runs) {
  dplyr::bind_rows(lapply(names(multisite_runs), function(cause_id) {
    dplyr::bind_rows(lapply(names(multisite_runs[[cause_id]]), function(scn) {
      multisite_runs[[cause_id]][[scn]]$proj_tbl
    }))
  })) |>
    dplyr::mutate(
      scenario = factor(as.character(scenario), levels = URUGUAY_SCENARIOS),
      scenario_label = factor(unname(URUGUAY_SCEN_LABELS[as.character(scenario)]),
                              levels = unname(URUGUAY_SCEN_LABELS[URUGUAY_SCENARIOS])),
      site_label = factor(site_public(cause_id, label), levels = names(URUGUAY_SITE_COLORS)),
      sex_label = factor(sex_public(sex), levels = c("Male", "Female", "Total")),
      projection_zone = factor(as.character(projection_zone), levels = c("historical", "credible", "caution", "risky", "beyond_max"))
    )
}

multisite_horizon_long <- function(multisite_runs) {
  dplyr::bind_rows(lapply(names(multisite_runs), function(cause_id) {
    dplyr::bind_rows(lapply(names(multisite_runs[[cause_id]]), function(scn) {
      multisite_runs[[cause_id]][[scn]]$horizon_tbl
    }))
  }))
}

get_multisite_common_endogenous_year <- function(multisite_runs) {
  h <- multisite_horizon_long(multisite_runs)
  vals <- h |>
    dplyr::filter(as.character(scenario) == "freeze", as.character(sex) == "T") |>
    dplyr::pull(max_projection_year_endogenous)
  vals <- suppressWarnings(as.integer(vals))
  vals <- vals[is.finite(vals)]
  if (!length(vals)) return(NA_integer_)
  min(vals, na.rm = TRUE)
}

build_multisite_observed_stack_data <- function(multisite_runs) {
  dplyr::bind_rows(lapply(names(multisite_runs), function(cause_id) {
    run_freeze <- multisite_runs[[cause_id]][["freeze"]]
    obs <- tryCatch(tibble::as_tibble(run_freeze$res$combined$obs_annual), error = function(e) tibble::tibble())
    if (!nrow(obs)) return(NULL)
    obs |>
      dplyr::transmute(
        cause_id = cause_id,
        label = run_freeze$label,
        site_label = site_public(cause_id, run_freeze$label),
        period = as.integer(period),
        deaths = as.numeric(obs)
      )
  }))
}

build_multisite_stack_data <- function(multisite_runs, scenarios = URUGUAY_SCENARIOS) {
  common_max <- get_multisite_common_endogenous_year(multisite_runs)
  proj <- multisite_projection_long(multisite_runs) |>
    dplyr::filter(
      as.character(scenario) %in% scenarios,
      metric == "mortality",
      series == "M|I|P",
      sex %in% c("M", "F"),
      period > 2022,
      is.na(common_max) | period <= common_max
    ) |>
    dplyr::group_by(scenario, scenario_label, cause_id, label, site_label, period) |>
    dplyr::summarise(deaths = sum(mean, na.rm = TRUE), .groups = "drop") |>
    dplyr::mutate(source = "Projected")
  obs <- build_multisite_observed_stack_data(multisite_runs) |>
    dplyr::filter(period >= 1998, period <= 2022) |>
    tidyr::crossing(scenario = scenarios) |>
    dplyr::mutate(
      scenario_label = factor(unname(URUGUAY_SCEN_LABELS[scenario]),
                              levels = unname(URUGUAY_SCEN_LABELS[URUGUAY_SCENARIOS])),
      source = "Observed"
    )
  dplyr::bind_rows(obs, proj) |>
    dplyr::mutate(
      site_label = factor(as.character(site_label), levels = names(URUGUAY_SITE_COLORS)),
      scenario = factor(as.character(scenario), levels = URUGUAY_SCENARIOS),
      scenario_label = factor(as.character(scenario_label), levels = unname(URUGUAY_SCEN_LABELS[URUGUAY_SCENARIOS])),
      common_endogenous_year = common_max
    ) |>
    dplyr::arrange(scenario, site_label, period)
}

plot_multisite_mortality_stack <- function(multisite_runs, scenarios = URUGUAY_SCENARIOS) {
  df <- build_multisite_stack_data(multisite_runs, scenarios = scenarios)
  ggplot2::ggplot(df, ggplot2::aes(period, deaths, fill = site_label, group = site_label)) +
    ggplot2::geom_area(alpha = 0.92, color = "white", linewidth = 0.08) +
    ggplot2::geom_vline(xintercept = 2022.5, color = "grey45", linewidth = 0.35) +
    ggplot2::facet_wrap(~ scenario_label, ncol = if (length(scenarios) > 2) 2 else length(scenarios)) +
    ggplot2::scale_fill_manual(values = URUGUAY_SITE_COLORS, name = "Cancer site", drop = TRUE) +
    ggplot2::scale_y_continuous(limits = c(0, NA), expand = ggplot2::expansion(mult = c(0, 0.03))) +
    ggplot2::labs(x = "Year", y = "Annual cancer deaths") +
    theme_paper_main(base_size = if (length(scenarios) > 2) 9.7 else 10.7) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.box = "vertical",
      strip.background = ggplot2::element_rect(fill = "gray95"),
      panel.spacing = grid::unit(0.8, "lines")
    )
}

export_multisite_stack_products <- function(multisite_result) {
  runs <- multisite_result$runs
  failures <- multisite_result$failures
  if (!is.data.frame(failures) || !ncol(failures)) {
    failures <- tibble::tibble(cause_id = character(), label = character(), error = character())
  }
  readr::write_csv(failures, file.path(OUT_SECTION5, "multisite_stack_failures.csv"))
  unlink(file.path(OUT_SECTION5, paste0("fig_uruguay_multisite_mortality_stack.", c("svg", "pdf"))))
  unlink(file.path(OUT_SECTION5, "fig_uruguay_multisite_mortality_stack_data.csv"))

  data_fq <- build_multisite_stack_data(runs, scenarios = c("freeze", "quit"))
  readr::write_csv(data_fq, file.path(OUT_SECTION5, "fig_uruguay_multisite_mortality_stack_freeze_quit_data.csv"))
  save_uruguay_plot(plot_multisite_mortality_stack(runs, scenarios = c("freeze", "quit")),
                    file.path(OUT_SECTION5, "fig_uruguay_multisite_mortality_stack_freeze_quit"),
                    width = 8.4, height = 4.8)

  common_max <- get_multisite_common_endogenous_year(runs)
  included <- sort(names(runs))
  recommendation <- c(
    "# Multisite Stack Recommendation",
    "",
    "Recommended main-text candidate: `fig_uruguay_multisite_mortality_stack_freeze_quit`.",
    "",
    "The two-panel version gives a cleaner baseline-versus-cessation contrast. The moderate-change scenarios are visually close to Frozen at 2022 in the stacked composition, so the four-panel version is not recommended for the main text.",
    "",
    sprintf("Both variants are clipped to the common endogenous projection year across included cancer sites: %s.", fmt_int(common_max)),
    "",
    sprintf("Included sites: %s.", paste(site_public(included, included), collapse = ", ")),
    "",
    if (nrow(failures)) {
      paste0("Failed or excluded sites: ", paste(sprintf("%s (%s)", failures$cause_id, failures$error), collapse = "; "), ".")
    } else {
      "No cancer site failed estimation in the multisite exercise."
    },
    "",
    "Conceptual role: the figure shows that the same SBAPC architecture can be used to summarize a broader smoking-attributable cancer burden, while retaining site composition rather than collapsing immediately to a single total.",
    "",
    "Source: Own elaboration."
  )
  writeLines(recommendation, file.path(OUT_SECTION5, "multisite_stack_recommendation.md"), useBytes = TRUE)
  invisible(multisite_result)
}

export_uruguay_transmission_inputs <- function(run_list,
                                               csv_out = file.path(OUT_APPENDIXD, "tab_uruguay_transmission_inputs.csv"),
                                               tex_out = file.path(OUT_APPENDIXD, "tab_uruguay_transmission_inputs.tex")) {
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

write_transmission_inputs_recommendation <- function() {
  lines <- c(
    "# Transmission Inputs Table Recommendation",
    "",
    "Recommendation: do not use the current lung-cancer transmission-input table as a central Section 5 float.",
    "",
    "Reason: the empirical section should foreground substantive projected mortality levels and the multisite extension. The compact input table is useful for documentation, but it does not carry a main empirical result.",
    "",
    "Placement: Appendix D, or Appendix B if the paper centralizes all external-link parameters in a methods appendix.",
    "",
    "Expansion required for a central methods table: sex-specific incidence relative risks, the risk-reversion schedule after cessation, post-diagnosis mortality interval probabilities, the annualized kernel total mass, and the maximum lag.",
    "",
    "Source: Own elaboration."
  )
  writeLines(lines, file.path(OUT_APPENDIXD, "transmission_inputs_table_recommendation.md"), useBytes = TRUE)
  invisible(TRUE)
}

write_uruguay_notes_and_inventories <- function() {
  section5_fig <- c(
    "# Figure Titles and Notes: Section 5",
    "",
    "## fig_uruguay_lung_mortality_by_sex",
    "Files: fig_uruguay_lung_mortality_by_sex.svg, fig_uruguay_lung_mortality_by_sex.pdf",
    "Title: Projected lung-cancer mortality in Uruguay",
    "Note: Black lines show observed annual lung-cancer deaths through 2022. Colored lines and ribbons show SBAPC projected annual deaths and interval summaries under smoking-prevalence scenarios, separately for Male and Female. Line type indicates the endogenous horizon region; projections are clipped to the common endogenous horizon. Source: Own elaboration.",
    "",
    "## fig_uruguay_multisite_mortality_stack_freeze_quit",
    "Files: fig_uruguay_multisite_mortality_stack_freeze_quit.svg, fig_uruguay_multisite_mortality_stack_freeze_quit.pdf",
    "Title: Projected smoking-attributable cancer mortality under baseline and cessation",
    "Note: Stacked areas report annual deaths summed over sexes for the nine smoking-attributable cancer sites under Frozen at 2022 and Quit. Site colors are constant across panels; projections are clipped to the common endogenous horizon across included sites. Source: Own elaboration."
  )
  writeLines(section5_fig, file.path(OUT_SECTION5, "figure_titles_notes.md"), useBytes = TRUE)

  section5_tab <- c(
    "# Table Titles and Notes: Section 5",
    "",
    "## tab_uruguay_lung_mortality_selected_years",
    "Title: Projected annual lung-cancer deaths in selected years",
    "Note: Values are annual deaths by sex and smoking-prevalence scenario. The 2022 column is the observed historical endpoint; future columns are SBAPC projections and are restricted to the common endogenous horizon. Source: Own elaboration."
  )
  writeLines(section5_tab, file.path(OUT_SECTION5, "table_titles_notes.md"), useBytes = TRUE)

  appendix_fig <- c(
    "# Figure Titles and Notes: Appendix D",
    "",
    "## fig_uruguay_mortality_effects",
    "Files: fig_uruguay_mortality_effects.svg, fig_uruguay_mortality_effects.pdf",
    "Title: Annual mortality effects of smoking-prevalence scenarios",
    "Note: Lines report projected annual lung-cancer deaths relative to Frozen at 2022, summed over sexes. Line type indicates the endogenous horizon region. Source: Own elaboration.",
    "",
    "## fig_uruguay_data_overview",
    "Files: fig_uruguay_data_overview.svg, fig_uruguay_data_overview.pdf",
    "Title: Historical empirical inputs for Uruguay",
    "Note: Smoking prevalence comes from weighted harmonized survey microdata over ages 20-65. Incidence and mortality are preprocessed historical inputs over ages 35-89 and are shown as crude annual rates per 100,000 population. Source: Own elaboration.",
    "",
    "## fig_uruguay_benchmark_comparison",
    "Files: fig_uruguay_benchmark_comparison.svg, fig_uruguay_benchmark_comparison.pdf",
    "Title: SBAPC and BAPC benchmark projections",
    "Note: The display contrasts scenario-responsive SBAPC projections with the scenario-blind BAPC benchmark over the common endogenous projection horizon. Source: Own elaboration.",
    "",
    "## fig_uruguay_horizon_support",
    "Files: fig_uruguay_horizon_support.svg, fig_uruguay_horizon_support.pdf",
    "Title: Endogenous support diagnostics for Uruguay",
    "Note: Lines report the exposure-weighted mean support fraction by projection year and sex. Thresholds mark the support levels used to classify credible, caution, risky, and beyond-maximum projection regions. Source: Own elaboration."
  )
  writeLines(appendix_fig, file.path(OUT_APPENDIXD, "figure_titles_notes.md"), useBytes = TRUE)

  appendix_tab <- c(
    "# Table Titles and Notes: Appendix D",
    "",
    "## tab_uruguay_horizon_boundaries",
    "Title: Endogenous horizon boundaries in the Uruguay application",
    "Note: Boundaries are computed from the Lexis-support diagnostics used to classify projection years. The last historical year is 2022 for all rows. Source: Own elaboration.",
    "",
    "## tab_uruguay_fit_scores",
    "Title: Historical fit diagnostics for the Uruguay freeze baseline",
    "Note: Fit statistics are reported for the freeze-baseline empirical run and should be interpreted as diagnostics, not as the main validation target. Source: Own elaboration.",
    "",
    "## tab_uruguay_fit_diagnostics_compact",
    "Title: Compact historical fit diagnostics",
    "Note: Delta information criteria compare SBAPC layers with the corresponding BAPC benchmark where such a benchmark is meaningful. Prevalence APC is reported as a standalone diagnostic. Source: Own elaboration.",
    "",
    "## tab_uruguay_cumulative_effects",
    "Title: Cumulative scenario effects in the Uruguay application",
    "Note: Effects are cumulative projected incident cases or deaths relative to Frozen at 2022, aggregated across sexes within each endogenous horizon region. Source: Own elaboration.",
    "",
    "## tab_uruguay_transmission_inputs",
    "Title: Lung-cancer transmission inputs used in the Uruguay application",
    "Note: Incidence relative risks and post-diagnosis mortality kernel summaries are fixed external inputs used by the empirical pipeline. This table is recommended as appendix documentation rather than a central empirical result. Source: Own elaboration."
  )
  writeLines(appendix_tab, file.path(OUT_APPENDIXD, "table_titles_notes.md"), useBytes = TRUE)

  section5_inv <- c(
    "# Section 5 Float Inventory",
    "",
    "| Filename | Location | Priority | Type | Aggregation | SVG/PDF | Data CSV | Source note | Purpose |",
    "|---|---|---|---|---|---|---|---|---|",
    "| `fig_uruguay_lung_mortality_by_sex` | Main text | Essential | Substantive | Sex-specific annual mortality counts | Yes | Yes | Yes | Main empirical lung-cancer mortality projection by sex and scenario. |",
    "| `tab_uruguay_lung_mortality_selected_years` | Main text | Useful | Substantive | Sex-specific selected-year annual deaths | Not applicable | Yes | Yes | Compact numerical companion to the main lung-cancer mortality figure. |",
    "| `fig_uruguay_multisite_mortality_stack_freeze_quit` | Main text | Essential | Substantive | Both-sex annual deaths by cancer site | Yes | Yes | Yes | Shows the broader nine-site potential of the framework and total burden composition under baseline and cessation. |"
  )
  writeLines(section5_inv, file.path(OUT_SECTION5, "section5_float_inventory.md"), useBytes = TRUE)

  appendix_inv <- c(
    "# Appendix D Float Inventory",
    "",
    "| Filename | Location | Priority | Type | Aggregation | SVG/PDF | Data CSV | Source note | Purpose |",
    "|---|---|---|---|---|---|---|---|---|",
    "| `fig_uruguay_mortality_effects` | Appendix D | Useful | Diagnostic | Both-sex annual effects | Yes | Yes | Yes | Shows mortality effects relative to Frozen at 2022. |",
    "| `fig_uruguay_data_overview` | Appendix D | Essential | Descriptive | Historical sex-specific annual aggregates | Yes | Yes | Yes | Describes the empirical inputs entering the Uruguay application. |",
    "| `fig_uruguay_horizon_support` | Appendix D | Useful | Diagnostic | Sex-specific support diagnostics | Yes | Yes | Yes | Displays support deterioration underlying the horizon categories. |",
    "| `fig_uruguay_benchmark_comparison` | Appendix D | Useful | Diagnostic | Sex-specific annual counts | Yes | Yes | Yes | Contrasts scenario-responsive SBAPC with the scenario-blind BAPC benchmark. |",
    "| `tab_uruguay_horizon_boundaries` | Appendix D | Useful | Diagnostic | Sex-specific and total frontier rows | Not applicable | Yes | Yes | Reports projection horizon boundary years. |",
    "| `tab_uruguay_fit_scores` | Appendix D | Optional | Diagnostic | Freeze-baseline fit statistics | Not applicable | Yes | Yes | Full historical fit/backtesting diagnostic. |",
    "| `tab_uruguay_fit_diagnostics_compact` | Appendix D | Useful | Diagnostic | Freeze-baseline fit statistics | Not applicable | Yes | Yes | Compact diagnostic comparison of SBAPC layers and BAPC benchmarks. |",
    "| `tab_uruguay_cumulative_effects` | Appendix D | Useful | Diagnostic | Both-sex cumulative effects | Not applicable | Yes | Yes | Scenario effects relative to Frozen at 2022 by endogenous horizon region. |",
    "| `tab_uruguay_transmission_inputs` | Appendix D | Optional | Descriptive | Sex-specific external inputs | Not applicable | Yes | Yes | Documents fixed external lung-cancer transmission inputs. |"
  )
  writeLines(appendix_inv, file.path(OUT_APPENDIXD, "appendixD_float_inventory.md"), useBytes = TRUE)

  invisible(TRUE)
}

export_uruguay_products <- function(run_list, inputs = NULL) {
  if (is.null(inputs)) inputs <- build_inputs_real_cause(run_cfg$causes_tbl[1, ])

  readr::write_csv(combined_projection_tables(run_list), file.path(URUGUAY_OUT_BASE, "uruguay_projection_long.csv"))
  readr::write_csv(observed_projection_tables(run_list), file.path(URUGUAY_OUT_BASE, "uruguay_observed_long.csv"))
  readr::write_csv(extract_smoking_exposure(run_list), file.path(URUGUAY_OUT_BASE, "uruguay_smoking_exposure_long.csv"))
  readr::write_csv(build_uruguay_mortality_effects(run_list), file.path(URUGUAY_OUT_BASE, "uruguay_mortality_effects_long.csv"))
  readr::write_csv(build_uruguay_lung_mortality_by_sex_data(run_list), file.path(OUT_SECTION5, "fig_uruguay_lung_mortality_by_sex_data.csv"))
  readr::write_csv(build_uruguay_mortality_effects(run_list), file.path(OUT_APPENDIXD, "fig_uruguay_mortality_effects_data.csv"))
  readr::write_csv(build_uruguay_data_overview(inputs), file.path(OUT_APPENDIXD, "fig_uruguay_data_overview_data.csv"))
  readr::write_csv(build_uruguay_benchmark_comparison_data(run_list), file.path(OUT_APPENDIXD, "fig_uruguay_benchmark_comparison_data.csv"))
  readr::write_csv(build_uruguay_horizon_support_data(run_list), file.path(OUT_APPENDIXD, "fig_uruguay_horizon_support_data.csv"))

  unlink(file.path(OUT_SECTION5, paste0(c(
    "fig_uruguay_smoking_exposure",
    "fig_uruguay_projection_panel",
    "fig_uruguay_mortality_effects"
  ), rep(c(".svg", ".pdf"), each = 3))))
  unlink(file.path(OUT_SECTION5, c(
    "tab_uruguay_cumulative_effects.csv",
    "tab_uruguay_cumulative_effects.tex",
    "tab_uruguay_transmission_inputs.csv",
    "tab_uruguay_transmission_inputs.tex"
  )))
  unlink(file.path(OUT_APPENDIXD, c(
    paste0("fig_uruguay_projection_panel.", c("svg", "pdf")),
    "fig_uruguay_projection_panel_data.csv",
    paste0("fig_uruguay_lung_mortality_uncertainty.", c("svg", "pdf")),
    "fig_uruguay_lung_mortality_uncertainty_data.csv"
  )))

  save_uruguay_plot(plot_uruguay_lung_mortality_by_sex(run_list),
                    file.path(OUT_SECTION5, "fig_uruguay_lung_mortality_by_sex"),
                    width = 7.4, height = 4.6)
  export_uruguay_lung_mortality_selected_years(run_list)

  save_uruguay_plot(plot_uruguay_mortality_effects(run_list),
                    file.path(OUT_APPENDIXD, "fig_uruguay_mortality_effects"),
                    width = 7.4, height = 4.6)

  save_uruguay_plot(plot_uruguay_data_overview(inputs),
                    file.path(OUT_APPENDIXD, "fig_uruguay_data_overview"),
                    width = 7.6, height = 7.6)
  save_uruguay_plot(plot_uruguay_benchmark_comparison(run_list),
                    file.path(OUT_APPENDIXD, "fig_uruguay_benchmark_comparison"),
                    width = 9.4, height = 6.9)
  save_uruguay_plot(plot_uruguay_horizon_support(run_list),
                    file.path(OUT_APPENDIXD, "fig_uruguay_horizon_support"),
                    width = 7.4, height = 5.4)

  export_uruguay_cumulative_effects(run_list)
  export_uruguay_transmission_inputs(run_list)
  export_uruguay_horizon_boundaries(run_list)
  export_uruguay_fit_scores(run_list)
  export_uruguay_fit_diagnostics_compact(run_list)
  write_transmission_inputs_recommendation()
  write_uruguay_notes_and_inventories()

  invisible(list(section5 = OUT_SECTION5, appendixD = OUT_APPENDIXD))
}

replicate_uruguay_empirical <- function(scenarios = URUGUAY_SCENARIOS,
                                        save_raw_rds = FALSE,
                                        run_multisite = TRUE) {
  inputs <- build_inputs_real_cause(run_cfg$causes_tbl[1, ])
  runs <- run_uruguay_lung_rebuilt_all_scenarios(run_cfg = run_cfg, scenarios = scenarios, save_raw_rds = save_raw_rds)
  export_uruguay_products(runs, inputs = inputs)
  multisite <- NULL
  if (isTRUE(run_multisite)) {
    multisite <- run_uruguay_multisite_rebuilt_all_scenarios(scenarios = scenarios, save_raw_rds = save_raw_rds)
    export_multisite_stack_products(multisite)
    write_uruguay_notes_and_inventories()
  }
  invisible(list(runs = runs, multisite = multisite, out_base = URUGUAY_OUT_BASE))
}

if (sys.nframe() == 0L) {
  replicate_uruguay_empirical()
}
