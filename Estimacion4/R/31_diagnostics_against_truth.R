# 31_diagnostics_against_truth.R
# Diagnostic tools to compare pipeline results against simulated ground truth (DGP)

compare_pipeline_to_truth <- function(res_both, sim, out_dir = NULL, prefix = "") {
  if (nzchar(prefix) && !grepl("_$", prefix)) prefix <- paste0(prefix, "_")
  if (is.null(out_dir)) {
    out_dir <- file.path(BAPC_PATHS$results, "truth_diagnostics")
  }
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  last_hist <- res_both$meta$last_hist %||% 2022
  
  # 1) Incidencia: Comparar tasas proyectadas
  # ------------------------------------------
  diag_inc <- list()
  for (sx in c("M", "F")) {
    res_sex <- if (sx == "M") res_both$resM else res_both$resF
    if (is.null(res_sex)) {
      message("compare_pipeline_to_truth: res_sex is NULL for ", sx)
      next
    }
    
    # Intentar obtener rates de inc_fit o inc_fit_bapc
    rates_df <- res_sex$inc_fit$rates_all %||% res_sex$inc_fit_bapc$rates_all
    if (is.null(rates_df)) {
      message("compare_pipeline_to_truth: rates_all is NULL for ", sx)
      next
    }
    
    # Estimación (Ponderada por exposición)
    hat_raw <- res_sex$inc_fit$rates_all_full
    if (is.null(hat_raw)) {
      message("compare_pipeline_to_truth: rates_all_full is NULL for Informed ", sx)
      next
    }
    
    hat <- hat_raw %>%
      dplyr::filter(age >= 35, age <= 89) %>%
      dplyr::group_by(period) %>%
      dplyr::summarise(
        cases_hat = sum(rate_hat * E, na.rm = TRUE),
        exposure_hat = sum(E, na.rm = TRUE),
        support_frac = stats::weighted.mean(support_frac, w = E, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::mutate(rate_hat = cases_hat / exposure_hat)
    
    # Verdad (DGP, filtrada por las mismas edades)
    true <- sim$inc_truth_grid %>%
      dplyr::filter(as.character(sex) == sx, age >= 35, age <= 89)
    
    if (!("exposure" %in% names(true))) {
      true <- true %>% dplyr::left_join(sim$pop_all, by = c("sex", "age", "period"))
    }

    true <- true %>%
      dplyr::group_by(period) %>%
      dplyr::summarise(
        cases_true = sum(rateI_scen_true * exposure, na.rm = TRUE),
        exposure_true = sum(exposure, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::mutate(rate_true = cases_true / exposure_true)

    # Benchmark (BAPC puro, sin prevalencia)
    bapc_raw <- res_sex$inc_fit_bapc$rates_all_full
    if (is.null(bapc_raw)) {
      message("compare_pipeline_to_truth: rates_all_full is NULL for BAPC ", sx)
      next
    }
    
    bapc_df <- bapc_raw %>%
      dplyr::filter(age >= 35, age <= 89) %>%
      dplyr::group_by(period) %>%
      dplyr::summarise(
        cases_bapc = sum(rate_hat * E, na.rm = TRUE),
        exposure_bapc = sum(E, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::mutate(rate_bapc = cases_bapc / exposure_bapc)
    
    if (nrow(hat) > 0 && nrow(true) > 0) {
      comp <- dplyr::inner_join(hat, true, by = "period") %>%
        dplyr::left_join(bapc_df %>% dplyr::select(period, rate_bapc), by = "period") %>%
        dplyr::mutate(sex = sx, error = rate_hat - rate_true, rel_error = error / pmax(rate_true, 1e-12))
      diag_inc[[sx]] <- comp
    } else {
      message("compare_pipeline_to_truth: empty hat or true for incidence ", sx)
    }
  }
  
  df_inc <- dplyr::bind_rows(diag_inc)
  if (nrow(df_inc) > 0) {
    df_inc$sex <- factor(df_inc$sex, levels = c("M", "F"))
    # Escalar para el gráfico (100k)
    df_plot_inc <- df_inc %>%
      dplyr::mutate(dplyr::across(c(rate_true, rate_hat, rate_bapc), ~ .x * 100000))
      
    p_inc <- ggplot(df_plot_inc, aes(x = period)) +
      geom_line(aes(y = rate_true, color = "Truth"), linetype = "dashed", linewidth = 1.1) +
      geom_line(aes(y = rate_hat, color = "Informed (SBAPC)"), linewidth = 1.3) +
      geom_line(aes(y = rate_bapc, color = "Benchmark (BAPC)"), linetype = "dotted", linewidth = 1.0) +
      geom_vline(xintercept = last_hist, linetype = "dotted", color = "gray40") +
      facet_wrap(~sex, scales = "free_y") +
      scale_color_manual(values = c("Truth" = "black", "Informed (SBAPC)" = "#CD5C5C", "Benchmark (BAPC)" = "#4682B4")) +
      labs(title = "Incidence Rate: Estimate vs Truth",
           subtitle = paste("Cause:", res_both$meta$cause_id %||% "simulated", "| Scale: per 100,000"),
           y = "Weighted Rate", color = "Source", x = "Year") +
      theme_minimal(base_family = "sans") +
      theme(legend.position = "bottom",
            panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = "white", color = NA),
            plot.background = element_rect(fill = "white", color = NA),
            strip.text = element_text(face = "bold", size = 11))
            
    ggsave(file.path(out_dir, paste0(prefix, "comparison_incidence_truth.png")), p_inc, 
           width = 10, height = 6, bg = "white", dpi = 300)
  }
  
  # 2) Mortalidad: Comparar muertes proyectadas
  # ------------------------------------------
  diag_mort <- list()
  for (sx in c("M", "F")) {
    res_sex <- if (sx == "M") res_both$resM else res_both$resF
    if (is.null(res_sex)) next
    
    # Estimación (Unir historia y futuro)
    obs_m <- res_sex$obs_annual %>% dplyr::mutate(deaths_hat = obs) %>% dplyr::select(period, deaths_hat)
    proj_m <- res_sex$annual_anchor %>% dplyr::filter(period > last_hist) %>% dplyr::select(period, deaths_hat)
    hat_m <- dplyr::bind_rows(obs_m, proj_m)
    
    true_m <- sim$mort_truth_grid %>%
      dplyr::filter(as.character(sex) == sx, age >= 35, age <= 89) %>%
      dplyr::group_by(period) %>%
      dplyr::summarise(deaths_true = sum(mort_deaths_scen_true, na.rm = TRUE), .groups = "drop")

    # Benchmark (BAPC puro)
    bapc_m <- res_sex$annual_bapc %>% dplyr::select(period, deaths_bapc = deaths_hat)
    
    # Informed-noP (Inc -> Mort, no Prev)
    nop_m <- res_sex$annual_anchor_noP %>% dplyr::select(period, deaths_noP = deaths_hat)
    
    if (nrow(hat_m) > 0 && nrow(true_m) > 0) {
      comp_m <- dplyr::inner_join(hat_m, true_m, by = "period") %>%
        dplyr::left_join(bapc_m, by = "period") %>%
        dplyr::left_join(nop_m, by = "period") %>%
        dplyr::mutate(sex = sx, error = deaths_hat - deaths_true)
      diag_mort[[sx]] <- comp_m
    }
  }
  
  df_mort <- dplyr::bind_rows(diag_mort)
  if (nrow(df_mort) > 0) {
    df_mort$sex <- factor(df_mort$sex, levels = c("M", "F"))
    p_mort <- ggplot(df_mort, aes(x = period)) +
      geom_line(aes(y = deaths_true, color = "Truth"), linetype = "dashed", linewidth = 1.1) +
      geom_line(aes(y = deaths_hat, color = "Informed (SBAPC)"), linewidth = 1.3) +
      geom_line(aes(y = deaths_noP, color = "Informed (No Prev)"), linetype = "dotdash", linewidth = 1.0) +
      geom_line(aes(y = deaths_bapc, color = "Benchmark (BAPC)"), linetype = "dotted", linewidth = 1.0) +
      geom_vline(xintercept = last_hist, linetype = "dotted", color = "gray40") +
      facet_wrap(~sex, scales = "free_y") +
      scale_color_manual(values = c("Truth" = "black", "Informed (SBAPC)" = "#CD5C5C", "Informed (No Prev)" = "#ff7f0e", "Benchmark (BAPC)" = "#4682B4")) +
      labs(title = "Mortality: Estimate vs Truth (Full Period)",
           subtitle = paste("Cause:", res_both$meta$cause_id %||% "simulated", "| Vertical line at", last_hist),
           y = "Total Deaths", color = "Source", x = "Year") +
      theme_minimal(base_family = "sans") +
      theme(legend.position = "bottom",
            panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = "white", color = NA),
            plot.background = element_rect(fill = "white", color = NA),
            strip.text = element_text(face = "bold", size = 11))
            
    ggsave(file.path(out_dir, paste0(prefix, "comparison_mortality_truth.png")), p_mort, 
           width = 10, height = 6, bg = "white", dpi = 300)
    
    metrics <- df_mort %>%
      dplyr::group_by(sex) %>%
      dplyr::summarise(
        hist_bias = 100 * mean(error[period <= last_hist], na.rm = TRUE) / pmax(mean(deaths_true[period <= last_hist], na.rm = TRUE), 1e-12),
        proj_bias = 100 * mean(error[period > last_hist], na.rm = TRUE) / pmax(mean(deaths_true[period > last_hist], na.rm = TRUE), 1e-12),
        .groups = "drop"
      )
    
    # Nuevo: Métricas desglosadas por confiabilidad (Reliability)
    # Necesitamos el support_frac que está en la incidencia
    # Mapeamos cada año de proyección a su categoría de confiabilidad predominante
    rel_map <- df_inc %>%
      dplyr::group_by(sex, period) %>%
      dplyr::summarise(
        avg_support = mean(support_frac, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        reliability = dplyr::case_when(
          period <= last_hist ~ "Historical",
          avg_support >= 0.50  ~ "Credible",
          avg_support >= 0.33  ~ "Caution",
          TRUE                 ~ "Risky"
        )
      )
    
    metrics_rel <- df_mort %>%
      dplyr::left_join(rel_map, by = c("sex", "period")) %>%
      dplyr::group_by(sex, reliability) %>%
      dplyr::summarise(
        bias = 100 * mean(deaths_hat - deaths_true, na.rm = TRUE) / pmax(mean(deaths_true, na.rm = TRUE), 1e-12),
        .groups = "drop"
      )
    
    write.csv(metrics, file.path(out_dir, paste0(prefix, "truth_comparison_metrics.csv")), row.names = FALSE)
    write.csv(metrics_rel, file.path(out_dir, paste0(prefix, "truth_comparison_reliability.csv")), row.names = FALSE)
    write.csv(df_inc, file.path(out_dir, paste0(prefix, "truth_comparison_incidence.csv")), row.names = FALSE)
    write.csv(df_mort, file.path(out_dir, paste0(prefix, "truth_comparison_mortality.csv")), row.names = FALSE)
  }
  
  # 3) Horizonte y Soporte: Extraer para el manuscrito
  # --------------------------------------------------
  diag_support <- list()
  for (sx in c("M", "F")) {
    res_sex <- if (sx == "M") res_both$resM else res_both$resF
    if (is.null(res_sex)) next
    
    # El soporte está en rates_all_full
    supp <- tryCatch(res_sex$inc_fit$rates_all_full, error = function(e) NULL)
    if (is.null(supp)) next
    
    supp_sum <- supp %>%
      dplyr::mutate(
        reliability = dplyr::case_when(
          period <= last_hist ~ "Historical",
          support_frac >= 0.50  ~ "Credible",
          support_frac >= 0.33  ~ "Caution",
          TRUE                 ~ "Risky"
        )
      ) %>%
      dplyr::group_by(reliability) %>%
      dplyr::summarise(
        n_cells = n(),
        avg_support_frac = mean(support_frac, na.rm = TRUE),
        min_support_frac = if(any(is.finite(support_frac))) min(support_frac, na.rm = TRUE) else NA_real_,
        n_clamped_period = sum(period_is_clamped %||% 0, na.rm = TRUE),
        n_edge_cohort = sum(cohort_is_edge %||% 0, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::mutate(sex = sx)
    
    diag_support[[sx]] <- supp_sum
  }
  
  df_supp <- dplyr::bind_rows(diag_support)
  if (nrow(df_supp) > 0) {
    write.csv(df_supp, file.path(out_dir, paste0(prefix, "truth_comparison_support.csv")), row.names = FALSE)
  }
  
  message("Diagnostics saved to: ", out_dir)
  return(list(inc = df_inc, mort = df_mort, metrics = if (exists("metrics")) metrics else NULL, support = df_supp))
}

export_latex_bias_summary <- function(metrics_df, file_out) {
  # Simple LaTeX tabular generator for paper bias summary
  if (!requireNamespace("tidyr", quietly = TRUE)) stop("tidyr required")
  
  # Agregamos por DGP y Escenario para el main text
  summary_tab <- metrics_df %>%
    dplyr::group_by(dgp, scenario, sex) %>%
    dplyr::summarise(
      Mean_Hist_Bias = mean(hist_bias, na.rm = TRUE),
      Mean_Proj_Bias = mean(proj_bias, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Formatear como tabla LaTeX básica
  lines <- c(
    "\\begin{tabular}{lllcc}",
    "\\hline",
    "DGP & Scenario & Sex & Hist Bias (%) & Proj Bias (%) \\\\",
    "\\hline"
  )
  
  for (i in 1:nrow(summary_tab)) {
    row <- summary_tab[i, ]
    lines <- c(lines, sprintf("%s & %s & %s & %.2f & %.2f \\\\", 
                              row$dgp, row$scenario, row$sex, row$Mean_Hist_Bias, row$Mean_Proj_Bias))
  }
  
  lines <- c(lines, "\\hline", "\\end{tabular}")
  
  writeLines(lines, file_out)
  message("LaTeX table saved to: ", file_out)
}
