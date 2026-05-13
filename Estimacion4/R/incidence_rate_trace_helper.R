trace_incidence_rate_chain <- function(run_obj,
                                       out_dir,
                                       prefix,
                                       last_hist = 2022L) {
  stopifnot(is.list(run_obj), is.character(out_dir), length(out_dir) == 1, is.character(prefix), length(prefix) == 1)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  `%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

  safe_num <- function(x) suppressWarnings(as.numeric(x))

  pick_truth <- function(sim, sex_lab) {
    if (is.null(sim) || !is.data.frame(sim$inc_truth_grid) || !is.data.frame(sim$pop_all)) {
      return(list(cells = data.frame(), annual = data.frame()))
    }
    cells <- sim$inc_truth_grid |>
      dplyr::filter(as.character(sex) == sex_lab, period > last_hist) |>
      dplyr::left_join(
        sim$pop_all |>
          dplyr::filter(as.character(sex) == sex_lab) |>
          dplyr::select(sex, age, period, exposure),
        by = c("sex", "age", "period")
      ) |>
      dplyr::mutate(
        exposure = safe_num(exposure),
        rate_true = safe_num(rateI_base_true),
        cases_true = exposure * rate_true
      ) |>
      dplyr::select(sex, age, period, rate_true, exposure, cases_true)

    annual <- cells |>
      dplyr::group_by(period) |>
      dplyr::summarise(
        rate_true_ann = stats::weighted.mean(rate_true, w = exposure, na.rm = TRUE),
        cases_true = sum(cases_true, na.rm = TRUE),
        exposure_total = sum(exposure, na.rm = TRUE),
        .groups = "drop"
      )

    list(cells = cells, annual = annual)
  }

  process_sex <- function(res_sex, sim, sex_lab) {
    full <- res_sex$inc_fit$rates_all_full %||% data.frame()
    annual_out <- res_sex$inc_annual_cond %||% data.frame()
    if (!is.data.frame(full)) full <- data.frame()
    if (!is.data.frame(annual_out)) annual_out <- data.frame()

    truth <- pick_truth(sim, sex_lab)

    cells <- full |>
      dplyr::filter(period > last_hist) |>
      dplyr::mutate(
        sex = as.character(sex),
        E = safe_num(E),
        logE = safe_num(logE),
        z_prev = safe_num(z_prev),
        inc_tech_offset = safe_num(inc_tech_offset),
        coef_fc_offset_I = safe_num(coef_fc_offset_I),
        coef_fc_offset_I_effective = safe_num(coef_fc_offset_I_effective),
        coef_fc_posthoc_adj = safe_num(coef_fc_posthoc_adj),
        rate_hat = safe_num(rate_hat),
        rate_lwr = safe_num(rate_lwr),
        rate_upr = safe_num(rate_upr),
        offset_mult = exp(dplyr::coalesce(inc_tech_offset, 0) + dplyr::coalesce(coef_fc_offset_I, 0)),
        rate_hat_plus_offset = rate_hat * offset_mult,
        log_rate_hat = log(pmax(rate_hat, 1e-12)),
        log_rate_hat_plus_offset = log(pmax(rate_hat_plus_offset, 1e-12))
      ) |>
      dplyr::left_join(truth$cells, by = c("sex", "age", "period")) |>
      dplyr::mutate(
        log_rate_true = log(pmax(rate_true, 1e-12)),
        log_gap = log_rate_hat - log_rate_true,
        log_gap_plus_offset = log_rate_hat_plus_offset - log_rate_true,
        rate_ratio = rate_hat / rate_true,
        rate_ratio_plus_offset = rate_hat_plus_offset / rate_true,
        cases_hat_from_cells = rate_hat * E,
        cases_hat_from_cells_plus_offset = rate_hat_plus_offset * E,
        cases_ratio_from_cells = cases_hat_from_cells / cases_true,
        cases_ratio_from_cells_plus_offset = cases_hat_from_cells_plus_offset / cases_true,
        hbin = dplyr::case_when(
          horizon <= 5 ~ "1_5",
          horizon <= 10 ~ "6_10",
          horizon <= 20 ~ "11_20",
          TRUE ~ "21p"
        )
      )

    annual_from_cells <- cells |>
      dplyr::group_by(period) |>
      dplyr::summarise(
        rate_hat_ann = stats::weighted.mean(rate_hat, w = E, na.rm = TRUE),
        rate_hat_plus_offset_ann = stats::weighted.mean(rate_hat_plus_offset, w = E, na.rm = TRUE),
        cases_hat_from_cells = sum(cases_hat_from_cells, na.rm = TRUE),
        cases_hat_from_cells_plus_offset = sum(cases_hat_from_cells_plus_offset, na.rm = TRUE),
        .groups = "drop"
      ) |>
      dplyr::left_join(annual_out |>
                         dplyr::filter(period > last_hist) |>
                         dplyr::transmute(period, cases_hat_output = safe_num(cases_hat)),
                       by = "period") |>
      dplyr::left_join(truth$annual, by = "period") |>
      dplyr::mutate(
        rate_ratio = rate_hat_ann / rate_true_ann,
        rate_ratio_plus_offset = rate_hat_plus_offset_ann / rate_true_ann,
        cases_ratio_from_cells = cases_hat_from_cells / cases_true,
        cases_ratio_from_cells_plus_offset = cases_hat_from_cells_plus_offset / cases_true,
        cases_ratio_output = cases_hat_output / cases_true,
        consistency_output = cases_hat_output / cases_hat_from_cells
      )

    summary_tbl <- dplyr::bind_rows(
      cells |>
        dplyr::summarise(
          scope = "cells",
          mean_rate_ratio = mean(rate_ratio, na.rm = TRUE),
          median_rate_ratio = stats::median(rate_ratio, na.rm = TRUE),
          p10_rate_ratio = as.numeric(stats::quantile(rate_ratio, 0.10, na.rm = TRUE)),
          p90_rate_ratio = as.numeric(stats::quantile(rate_ratio, 0.90, na.rm = TRUE)),
          mean_rate_ratio_plus_offset = mean(rate_ratio_plus_offset, na.rm = TRUE),
          mean_log_gap = mean(log_gap, na.rm = TRUE),
          mean_log_gap_plus_offset = mean(log_gap_plus_offset, na.rm = TRUE),
          corr_log = stats::cor(log_rate_hat, log_rate_true, use = "complete.obs"),
          corr_log_plus_offset = stats::cor(log_rate_hat_plus_offset, log_rate_true, use = "complete.obs")
        ),
      annual_from_cells |>
        dplyr::summarise(
          scope = "annual",
          mean_rate_ratio = mean(rate_ratio, na.rm = TRUE),
          median_rate_ratio = stats::median(rate_ratio, na.rm = TRUE),
          p10_rate_ratio = as.numeric(stats::quantile(rate_ratio, 0.10, na.rm = TRUE)),
          p90_rate_ratio = as.numeric(stats::quantile(rate_ratio, 0.90, na.rm = TRUE)),
          mean_rate_ratio_plus_offset = mean(rate_ratio_plus_offset, na.rm = TRUE),
          mean_log_gap = mean(log(pmax(rate_hat_ann,1e-12)) - log(pmax(rate_true_ann,1e-12)), na.rm = TRUE),
          mean_log_gap_plus_offset = mean(log(pmax(rate_hat_plus_offset_ann,1e-12)) - log(pmax(rate_true_ann,1e-12)), na.rm = TRUE),
          corr_log = stats::cor(log(pmax(rate_hat_ann,1e-12)), log(pmax(rate_true_ann,1e-12)), use = "complete.obs"),
          corr_log_plus_offset = stats::cor(log(pmax(rate_hat_plus_offset_ann,1e-12)), log(pmax(rate_true_ann,1e-12)), use = "complete.obs")
        )
    ) |>
      dplyr::mutate(sex = sex_lab, .before = 1)

    horizon_tbl <- cells |>
      dplyr::group_by(hbin) |>
      dplyr::summarise(
        mean_rate_ratio = mean(rate_ratio, na.rm = TRUE),
        mean_rate_ratio_plus_offset = mean(rate_ratio_plus_offset, na.rm = TRUE),
        .groups = "drop"
      ) |>
      dplyr::mutate(sex = sex_lab, .before = 1)

    list(cells = cells, annual = annual_from_cells, summary = summary_tbl, horizon = horizon_tbl)
  }

  sxM <- process_sex(run_obj$res$resM, run_obj$sim, "M")
  sxF <- process_sex(run_obj$res$resF, run_obj$sim, "F")

  write.csv(dplyr::bind_rows(sxM$cells, sxF$cells), file.path(out_dir, paste0(prefix, "__inc_rate_trace_cells.csv")), row.names = FALSE)
  write.csv(dplyr::bind_rows(sxM$annual, sxF$annual), file.path(out_dir, paste0(prefix, "__inc_rate_trace_annual.csv")), row.names = FALSE)
  write.csv(dplyr::bind_rows(sxM$summary, sxF$summary), file.path(out_dir, paste0(prefix, "__inc_rate_trace_summary.csv")), row.names = FALSE)
  write.csv(dplyr::bind_rows(sxM$horizon, sxF$horizon), file.path(out_dir, paste0(prefix, "__inc_rate_trace_horizon.csv")), row.names = FALSE)

  invisible(list(M = sxM, F = sxF))
}
