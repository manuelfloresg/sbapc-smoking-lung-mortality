
# =============================================================
# Main-text paper figure helpers
# =============================================================

paper_base_family <- function(preferred = c("CMU Serif", "Computer Modern Serif", "Latin Modern Roman", "serif")) {
  opt_family <- getOption("BAPC_PAPER_BASE_FAMILY", NULL)
  if (is.character(opt_family) && length(opt_family) && nzchar(opt_family[[1]])) {
    return(opt_family[[1]])
  }
  # Use a robust default for PDF export; users can opt in to a specific CM family
  # via options(BAPC_PAPER_BASE_FAMILY = "CMU Serif") once their device supports it.
  "serif"
}

theme_paper_main <- function(base_size = 10.5,
                             base_family = paper_base_family(),
                             legend_position = "bottom") {
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "plain", size = base_size * 1.20, hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = base_size * 0.90, hjust = 0.5),
      plot.caption = ggplot2::element_text(size = base_size * 0.78, hjust = 0),
      axis.title = ggplot2::element_text(size = base_size * 0.95),
      axis.text = ggplot2::element_text(size = base_size * 0.85),
      strip.text = ggplot2::element_text(size = base_size * 0.94),
      legend.position = legend_position,
      legend.box = "vertical",
      legend.direction = "horizontal",
      legend.box.just = "center",
      legend.title = ggplot2::element_text(size = base_size * 0.88),
      legend.text = ggplot2::element_text(size = base_size * 0.82),
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      legend.box.margin = ggplot2::margin(0, 0, 0, 0),
      plot.margin = ggplot2::margin(6, 7, 5, 6)
    )
}

paper_scenario_palette <- function(scenarios = names(scenario_labels_en), palette = c("auto", "method", "uruguay9")) {
  palette <- match.arg(palette)
  scenarios_chr <- as.character(scenarios)
  if (identical(palette, "auto")) {
    palette <- if ("down3pc" %in% scenarios_chr) "uruguay9" else "method"
  }
  pal_method <- c(
    "freeze"  = "#111111",
    "up1pc"   = "#C0392B",
    "down1pc" = "#1E8449",
    "quit"    = "#2874A6"
  )
  pal_uy9 <- c(
    "freeze"  = "#C0392B",
    "down1pc" = "#D4AC0D",
    "down3pc" = "#1E8449",
    "quit"    = "#2874A6"
  )
  pal <- if (identical(palette, "method")) pal_method else pal_uy9
  pal[intersect(names(pal), scenarios_chr)]
}

paper_zone_linetypes <- function() {
  c(
    "credible"   = "solid",
    "caution"    = "longdash",
    "risky"      = "dotted",
    "beyond_max" = "dotdash"
  )
}

register_figure_output <- function(file,
                                   registry_file,
                                   figure_group,
                                   panel_id = NA_character_,
                                   output_mode = c("qc", "paper_panel", "paper_figure"),
                                   note_template = NA_character_,
                                   extra = NULL) {
  output_mode <- match.arg(output_mode)
  dir.create(dirname(registry_file), recursive = TRUE, showWarnings = FALSE)
  row <- tibble::tibble(
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    figure_group = figure_group,
    panel_id = panel_id,
    output_mode = output_mode,
    file_path = normalizePath(file, winslash = "/", mustWork = FALSE),
    note_template = note_template
  )
  if (!is.null(extra) && length(extra)) {
    row <- dplyr::bind_cols(row, tibble::as_tibble(extra))
  }
  row <- dplyr::mutate(row, dplyr::across(dplyr::everything(), as.character))
  old <- if (file.exists(registry_file)) {
    tryCatch(
      readr::read_csv(
        registry_file,
        col_types = readr::cols(.default = readr::col_character()),
        show_col_types = FALSE
      ),
      error = function(e) NULL
    )
  } else NULL
  if (!is.null(old)) {
    old <- dplyr::mutate(old, dplyr::across(dplyr::everything(), as.character))
  }
  out <- dplyr::bind_rows(old, row)
  readr::write_csv(out, registry_file)
  invisible(row)
}


.paper_save_plot <- function(plot, file, width, height, dpi = 300) {
  ext <- tolower(tools::file_ext(file))
  if (identical(ext, "pdf") && capabilities("cairo")) {
    ggplot2::ggsave(
      filename = file, plot = plot, width = width, height = height,
      dpi = dpi, bg = "white", device = grDevices::cairo_pdf
    )
  } else {
    ggplot2::ggsave(
      filename = file, plot = plot, width = width, height = height,
      dpi = dpi, bg = "white"
    )
  }
}

save_plot_qc <- function(plot, file, width = 8.5, height = 5.2, dpi = 300,
                         registry_file = NULL, figure_group = NA_character_, panel_id = NA_character_,
                         note_template = NA_character_, extra = NULL) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  .paper_save_plot(plot = plot, file = file, width = width, height = height, dpi = dpi)
  if (!is.null(registry_file)) {
    register_figure_output(file, registry_file = registry_file, figure_group = figure_group,
                           panel_id = panel_id, output_mode = "qc",
                           note_template = note_template, extra = extra)
  }
  invisible(file)
}

save_plot_panel <- function(plot, file, width = 4.6, height = 3.6, dpi = 300,
                            registry_file = NULL, figure_group = NA_character_, panel_id = NA_character_,
                            note_template = NA_character_, extra = NULL) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  .paper_save_plot(plot = plot, file = file, width = width, height = height, dpi = dpi)
  if (!is.null(registry_file)) {
    register_figure_output(file, registry_file = registry_file, figure_group = figure_group,
                           panel_id = panel_id, output_mode = "paper_panel",
                           note_template = note_template, extra = extra)
  }
  invisible(file)
}

save_plot_figure <- function(plot, file, width = 9.2, height = 6.8, dpi = 300,
                             registry_file = NULL, figure_group = NA_character_, panel_id = NA_character_,
                             note_template = NA_character_, extra = NULL) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  .paper_save_plot(plot = plot, file = file, width = width, height = height, dpi = dpi)
  if (!is.null(registry_file)) {
    register_figure_output(file, registry_file = registry_file, figure_group = figure_group,
                           panel_id = panel_id, output_mode = "paper_figure",
                           note_template = note_template, extra = extra)
  }
  invisible(file)
}

.prepare_scenario_projection_segments <- function(df) {
  stopifnot(all(c("scenario", "period", "mean", "lwr", "upr", "projection_zone") %in% names(df)))
  dplyr::bind_rows(lapply(split(df, df$scenario), function(dd) {
    dd <- dd[order(dd$period), , drop = FALSE]
    if (!nrow(dd)) return(dd)
    zone <- as.character(dd$projection_zone)
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

.empirical_obs_df <- function(run_out, metric = c("incidence", "mortality"), sex = c("M", "F")) {
  metric <- match.arg(metric)
  sex <- match.arg(sex)
  rs <- if (sex == "M") run_out$res$resM else run_out$res$resF
  if (metric == "incidence") {
    rs$inc_obs_annual %>% dplyr::transmute(period = period, obs = obs)
  } else {
    rs$obs_annual %>% dplyr::transmute(period = period, obs = obs)
  }
}

.build_empirical_panel_data <- function(run_list,
                                       metric = c("incidence", "mortality"),
                                       sex = c("M", "F"),
                                       scenarios = NULL) {
  metric <- match.arg(metric)
  sex <- match.arg(sex)
  if (is.null(scenarios)) scenarios <- names(run_list)
  scenarios <- intersect(as.character(scenarios), names(run_list))
  if (!length(scenarios)) stop("No matching scenarios found in run_list.")
  series_target <- if (metric == "incidence") "I|P" else "M|I|P"
  ref <- run_list[[scenarios[[1]]]]
  obs <- .empirical_obs_df(ref, metric = metric, sex = sex)
  last_hist_year <- suppressWarnings(as.integer(ref$horizon_tbl$last_hist_year[ref$horizon_tbl$sex == sex][1]))
  proj <- dplyr::bind_rows(lapply(scenarios, function(scn) {
    x <- run_list[[scn]]$proj_tbl %>%
      dplyr::filter(sex == !!sex, metric == !!metric, series == !!series_target, period > !!last_hist_year) %>%
      dplyr::mutate(scenario = scn)
    y0 <- obs %>% dplyr::filter(period == last_hist_year) %>% dplyr::slice_tail(n = 1)
    if (nrow(y0)) {
      x <- dplyr::bind_rows(
        tibble::tibble(
          cause_id = unique(x$cause_id)[1],
          label = unique(x$label)[1],
          sex = sex,
          metric = metric,
          series = series_target,
          period = y0$period,
          mean = y0$obs,
          lwr = NA_real_,
          upr = NA_real_,
          projection_zone = factor("credible", levels = c("historical", "credible", "caution", "risky", "beyond_max")),
          scenario = scn
        ),
        x
      )
    }
    x
  }))
  proj <- proj %>%
    dplyr::mutate(
      scenario = factor(scenario, levels = scenarios, labels = unname(scenario_labels_en[scenarios])),
      projection_zone = factor(as.character(projection_zone), levels = c("credible", "caution", "risky", "beyond_max"))
    )
  proj <- .prepare_scenario_projection_segments(proj)
  list(obs = obs, proj = proj, last_hist_year = last_hist_year)
}

plot_empirical_scenario_panel <- function(run_list,
                                          metric = c("incidence", "mortality"),
                                          sex = c("M", "F"),
                                          scenarios = NULL,
                                          mode = c("paper_panel", "qc"),
                                          show_ci = FALSE,
                                          panel_title = NULL,
                                          base_size = 10.5,
                                          base_family = paper_base_family()) {
  metric <- match.arg(metric)
  sex <- match.arg(sex)
  mode <- match.arg(mode)
  dat <- .build_empirical_panel_data(run_list, metric = metric, sex = sex, scenarios = scenarios)
  obs <- dat$obs
  proj <- dat$proj
  last_hist_year <- dat$last_hist_year
  sex_lab <- dplyr::recode(sex, "M" = "Men", "F" = "Women")
  metric_lab <- dplyr::recode(metric, "incidence" = "Incidence", "mortality" = "Mortality")
  if (is.null(panel_title)) panel_title <- paste(metric_lab, "—", sex_lab)
  pal_raw <- paper_scenario_palette(scenarios = scenarios, palette = "auto")
  pal <- setNames(unname(pal_raw[scenarios]), unname(scenario_labels_en[scenarios]))
  pal <- pal[levels(proj$scenario)]
  g <- ggplot2::ggplot() +
    ggplot2::geom_line(data = obs, ggplot2::aes(x = period, y = obs), linewidth = 0.8, color = "black") +
    ggplot2::geom_vline(xintercept = last_hist_year + 0.5, linetype = 3, linewidth = 0.45, color = "grey35")
  if (isTRUE(show_ci)) {
    g <- g + ggplot2::geom_ribbon(
      data = proj,
      ggplot2::aes(x = period, ymin = lwr, ymax = upr, fill = scenario, group = interaction(scenario, zone_group)),
      alpha = 0.10, colour = NA, show.legend = FALSE
    )
  }
  g <- g +
    ggplot2::geom_line(
      data = proj,
      ggplot2::aes(x = period, y = mean, color = scenario, linetype = projection_zone,
                   group = interaction(scenario, zone_group)),
      linewidth = 0.90
    ) +
    ggplot2::scale_color_manual(values = pal, name = "Smoking prevalence scenario") +
    ggplot2::scale_linetype_manual(
      values = paper_zone_linetypes(),
      breaks = c("credible", "caution", "risky"),
      labels = c("Credible", "Caution", "Risky"),
      name = "Projection reliability"
    ) +
    ggplot2::guides(
      linetype = ggplot2::guide_legend(order = 1, nrow = 1, byrow = TRUE),
      color = ggplot2::guide_legend(order = 2, nrow = 1, byrow = TRUE),
      fill = ggplot2::guide_legend(order = 2, nrow = 1, byrow = TRUE)
    ) +
    ggplot2::labs(x = NULL, y = if (metric == "incidence") "Cases" else "Deaths") +
    ggplot2::scale_y_continuous(limits = c(0, NA), expand = c(0, 0)) +
    theme_paper_main(base_size = base_size, base_family = base_family)
  if (isTRUE(show_ci)) {
    g <- g + ggplot2::scale_fill_manual(values = pal, name = "Smoking prevalence scenario")
  }
  if (mode == "qc") {
    max_year <- suppressWarnings(max(proj$period, na.rm = TRUE))
    cap <- if (isTRUE(show_ci)) "Shaded bands show 95% credible intervals." else NULL
    g <- g + ggplot2::labs(
      title = paste0(panel_title, ": observed and projected"),
      subtitle = paste0("Observed through ", last_hist_year, " | Projected through ", max_year),
      caption = cap
    )
  } else {
    g <- g + ggplot2::labs(title = panel_title) +
      ggplot2::theme(
        plot.subtitle = ggplot2::element_blank(),
        plot.caption = ggplot2::element_blank()
      )
  }
  g
}

assemble_empirical_lung_main_figure <- function(run_list,
                                                scenarios = NULL,
                                                mode = c("paper_figure", "qc"),
                                                show_ci = FALSE,
                                                base_size = 10.5,
                                                base_family = paper_base_family()) {
  mode <- match.arg(mode)
  if (is.null(scenarios)) scenarios <- names(run_list)
  panel_mode <- if (mode == "paper_figure") "paper_panel" else "qc"
  p1 <- plot_empirical_scenario_panel(run_list, metric = "incidence", sex = "M", scenarios = scenarios,
                                      mode = panel_mode, show_ci = show_ci, panel_title = "Incidence — Men",
                                      base_size = base_size, base_family = base_family)
  p2 <- plot_empirical_scenario_panel(run_list, metric = "incidence", sex = "F", scenarios = scenarios,
                                      mode = panel_mode, show_ci = show_ci, panel_title = "Incidence — Women",
                                      base_size = base_size, base_family = base_family)
  p3 <- plot_empirical_scenario_panel(run_list, metric = "mortality", sex = "M", scenarios = scenarios,
                                      mode = panel_mode, show_ci = show_ci, panel_title = "Mortality — Men",
                                      base_size = base_size, base_family = base_family)
  p4 <- plot_empirical_scenario_panel(run_list, metric = "mortality", sex = "F", scenarios = scenarios,
                                      mode = panel_mode, show_ci = show_ci, panel_title = "Mortality — Women",
                                      base_size = base_size, base_family = base_family)
  g <- (p1 + p2) / (p3 + p4) +
    patchwork::plot_layout(guides = "collect") &
    ggplot2::theme(legend.position = "bottom", legend.box = "vertical")
  if (mode == "qc") {
    g <- g + patchwork::plot_annotation(
      title = "Lung cancer projections under smoking prevalence scenarios",
      subtitle = if (isTRUE(show_ci)) "Observed segment in black; line type indicates projection reliability zone; bands show 95% credible intervals." else "Observed segment in black; line type indicates projection reliability zone."
    )
  }
  g
}


export_empirical_lung_panel_set <- function(run_list,
                                            out_dir,
                                            scenarios = NULL,
                                            show_ci = FALSE,
                                            registry_file = file.path(out_dir, "figure_registry.csv"),
                                            stem = "panel_empirical_lung") {
  if (is.null(scenarios)) scenarios <- names(run_list)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  panels <- list(
    A = plot_empirical_scenario_panel(run_list, metric = "incidence", sex = "M", scenarios = scenarios, mode = "paper_panel", show_ci = show_ci, panel_title = "Incidence — Men"),
    B = plot_empirical_scenario_panel(run_list, metric = "incidence", sex = "F", scenarios = scenarios, mode = "paper_panel", show_ci = show_ci, panel_title = "Incidence — Women"),
    C = plot_empirical_scenario_panel(run_list, metric = "mortality", sex = "M", scenarios = scenarios, mode = "paper_panel", show_ci = show_ci, panel_title = "Mortality — Men"),
    D = plot_empirical_scenario_panel(run_list, metric = "mortality", sex = "F", scenarios = scenarios, mode = "paper_panel", show_ci = show_ci, panel_title = "Mortality — Women")
  )
  note <- paper_figure_note_template(show_ci = show_ci)
  files <- lapply(names(panels), function(tag) {
    file <- file.path(out_dir, paste0(stem, "_", tag, ".pdf"))
    save_plot_panel(panels[[tag]], file,
                    registry_file = registry_file,
                    figure_group = "main_empirical_lung",
                    panel_id = tag,
                    note_template = note,
                    extra = list(scenarios = paste(scenarios, collapse = ";"), show_ci = show_ci))
    file
  })
  names(files) <- names(panels)
  invisible(list(panels = panels, files = unlist(files), registry_file = registry_file))
}

paper_figure_note_template <- function(show_ci = FALSE) {
  if (isTRUE(show_ci)) {
    "Observed segments are shown in black. Line type indicates the endogenous projection-reliability zone (credible, caution, risky). Shaded bands show 95% credible intervals."
  } else {
    "Observed segments are shown in black. Line type indicates the endogenous projection-reliability zone (credible, caution, risky)."
  }
}

export_empirical_lung_main_figure <- function(run_list,
                                              out_dir,
                                              scenarios = NULL,
                                              show_ci = FALSE,
                                              registry_file = file.path(out_dir, "figure_registry.csv"),
                                              stem = "fig_main_empirical_lung") {
  if (is.null(scenarios)) scenarios <- names(run_list)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  fig <- assemble_empirical_lung_main_figure(run_list, scenarios = scenarios, mode = "paper_figure", show_ci = show_ci)
  note <- paper_figure_note_template(show_ci = show_ci)
  fig_file <- file.path(out_dir, paste0(stem, ".pdf"))
  save_plot_figure(fig, fig_file, registry_file = registry_file,
                   figure_group = "main_empirical_lung", panel_id = "A-D",
                   note_template = note,
                   extra = list(scenarios = paste(scenarios, collapse = ";"), show_ci = show_ci))
  invisible(list(plot = fig, file = fig_file, registry_file = registry_file))
}
