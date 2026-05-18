# =============================================================
# QC helpers (function-only; no side effects on source)
# =============================================================

to_num <- function(x) suppressWarnings(as.numeric(x))

# --- Parameter QC flags ---
qc_params_flags <- function(params_tbl, LI_top = L_I_MAX_YEARS) {
  params_tbl %>%
    dplyr::mutate(
      flag_rr_inc_invalid = is.na(rr_inc) | (rr_inc <= 1),
      flag_LI_en_tope     = !is.na(L_I) & (L_I >= LI_top),
      flag_betaI_distante = is.na(beta_I_mean) | abs(beta_I_mean - 1) > 0.30
    )
}

# --- Projection QC on main series ---
qc_proj_flags <- function(proj_tbl, growth_ratio_thr = 2, end_width_rel_thr = 1.0) {
  if (is.null(proj_tbl) || !nrow(proj_tbl)) return(tibble::tibble())

  proj_use <- proj_tbl %>%
    dplyr::filter(
      projection_zone %in% c("credible", "caution", "risky", "beyond_max"),
      (metric == "incidence" & series %in% c("I", "I|P")) |
        (metric == "mortality" & series == "M|I|P")
    ) %>%
    dplyr::arrange(cause_id, sex, metric, series, period)

  if (!nrow(proj_use)) return(tibble::tibble())

  by_series <- proj_use %>%
    dplyr::group_by(cause_id, label, sex, metric, series) %>%
    dplyr::summarise(
      flag_CI_mal_ordenado = any((lwr > mean) | (mean > upr), na.rm = TRUE),
      flag_CI_muy_ancho = {
        width_rel <- (upr - lwr) / pmax(abs(mean), 1)
        stats::median(width_rel[is.finite(width_rel)], na.rm = TRUE) > 0.8
      },
      flag_CI_muy_ancho_end = {
        width_rel <- (upr - lwr) / pmax(abs(mean), 1)
        ww <- width_rel[projection_zone %in% c("credible", "caution")]
        if (!length(ww) || all(!is.finite(ww))) ww <- width_rel
        ww <- ww[is.finite(ww)]
        if (!length(ww)) FALSE else dplyr::last(ww) > end_width_rel_thr
      },
      flag_salto_inicial = {
        mm <- mean
        if (length(mm) < 3) {
          FALSE
        } else {
          deltas <- diff(mm)
          first_delta <- deltas[1]
          rest_delta <- deltas[-1]
          thr <- if (!length(rest_delta) || all(!is.finite(rest_delta))) NA_real_ else 3 * stats::IQR(abs(rest_delta), na.rm = TRUE)
          rel_thr <- 0.50 * pmax(abs(mm[1]), 1)
          isTRUE(is.finite(first_delta) & ((is.finite(thr) & abs(first_delta) > thr) | abs(first_delta) > rel_thr))
        }
      },
      growth_ratio_credible = {
        mm <- mean[projection_zone == "credible"]
        mm <- mm[is.finite(mm)]
        if (length(mm) >= 2 && abs(mm[1]) > 1e-12) dplyr::last(mm) / mm[1] else NA_real_
      },
      flag_growth_extreme = dplyr::first(metric) == "incidence" & is.finite(growth_ratio_credible) & growth_ratio_credible > growth_ratio_thr,
      .groups = "drop"
    ) %>%
    dplyr::mutate(flag_key = dplyr::case_when(
      metric == "incidence" & series == "I" ~ "inc_apc",
      metric == "incidence" & series == "I|P" ~ "inc_main",
      metric == "mortality" & series == "M|I|P" ~ "mort_main",
      TRUE ~ NA_character_
    )) %>%
    dplyr::filter(!is.na(flag_key)) %>%
    dplyr::select(-metric, -series)

  if (!nrow(by_series)) return(tibble::tibble())

  by_series %>%
    tidyr::pivot_wider(
      id_cols = c(cause_id, label, sex),
      names_from = flag_key,
      values_from = c(flag_CI_mal_ordenado, flag_CI_muy_ancho, flag_CI_muy_ancho_end, flag_salto_inicial, flag_growth_extreme, growth_ratio_credible),
      names_glue = "{flag_key}_{.value}"
    )
}

# --- Consolidar flags por causa/sexo (resumen amigable) ---
make_qc_summary <- function(params_tbl, proj_tbl) {
  if (is.null(params_tbl) || nrow(params_tbl) == 0) return(tibble::tibble())

  need_cols <- c("beta_mode", "rr_inc", "prev_sign", "L_I", "bridge_years")
  for (nm in need_cols) if (!nm %in% names(params_tbl)) params_tbl[[nm]] <- NA_real_

  params0 <- params_tbl %>%
    dplyr::mutate(
      rr_inc       = to_num(rr_inc),
      prev_sign    = to_num(prev_sign),
      L_I          = to_num(L_I),
      bridge_years = suppressWarnings(as.integer(bridge_years))
    )

  qc <- params0 %>%
    dplyr::mutate(
      flag_rr_inc_missing = is.na(rr_inc),
      flag_rr_inc_invalid = is.na(rr_inc) | (rr_inc <= 1),
      flag_prev_sign_bad = is.na(prev_sign) | !(prev_sign %in% c(-1, 1)),
      flag_bridge_neg    = is.na(bridge_years) | (bridge_years < 0)
    )

  opt_flags <- c("flag_LI_en_tope", "flag_betaI_distante")
  for (nm in opt_flags) if (!nm %in% names(qc)) qc[[nm]] <- FALSE

  proj_flags <- tryCatch(qc_proj_flags(proj_tbl), error = function(e) tibble::tibble())

  qc %>%
    dplyr::left_join(proj_flags, by = c("cause_id", "label", "sex")) %>%
    dplyr::mutate(across(dplyr::starts_with("flag_"), ~ dplyr::coalesce(as.logical(.x), FALSE))) %>%
    dplyr::mutate(qc_score = rowSums(dplyr::across(dplyr::starts_with("flag_")), na.rm = TRUE)) %>%
    dplyr::select(cause_id, label, sex,
                  beta_mode, rr_inc, prev_sign, L_I, bridge_years,
                  dplyr::starts_with("flag_"), dplyr::contains("growth_ratio_credible"), qc_score) %>%
    dplyr::arrange(dplyr::desc(qc_score), cause_id, sex)
}

write_qc_outputs <- function(params_tbl,
                             proj_tbl,
                             out_file = file.path(DIAGNOSTICS_RESULTS_DIR, "qc_flags_by_cause.csv"),
                             print_top = TRUE,
                             top_n = 12) {
  qc_tbl <- make_qc_summary(params_tbl, proj_tbl)
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(qc_tbl, out_file)

  if (isTRUE(print_top) && nrow(qc_tbl) > 0) {
    cat("\n=== QC: top causes with most flags ===\n")
    print(
      qc_tbl %>%
        dplyr::select(
          dplyr::any_of(c("cause_id", "label", "sex", "qc_score",
                          "flag_rr_inc_invalid", "flag_LI_en_tope", "flag_betaI_distante")),
          dplyr::starts_with("flag_inc_main"),
          dplyr::starts_with("flag_mort_main"),
          dplyr::starts_with("flag_inc_apc")
        ) %>%
        dplyr::slice_head(n = top_n)
    )
  }

  invisible(qc_tbl)
}
