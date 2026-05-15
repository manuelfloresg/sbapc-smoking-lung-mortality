# 6) Wrapper para ambos sexos y suma total
# =============================================================
source("R/04b_rebuilder_helpers.R")

run_pipeline_both <- function(
    anchor_pseudo_w = ANCHOR_PSEUDO_W, sd_cohfix = SD_COHORT_RESID, sd_beta = SD_BETA_FIXED,
    age_min_m = AGE_M_MIN, age_max_m = AGE_M_MAX,
    age_min_p = AGE_P_MIN, age_max_p = AGE_P_MAX,
    age_min_i = AGE_I_MIN, age_max_i = AGE_I_MAX,
    L_I = L_I_DEFAULT, Da_I = DA_I,
    mort_trend_scenario = MORT_TREND_SCENARIO, delta_tech = DELTA_TECH,
    inc_include_trend = (INC_TREND_DEGREE > 0),
    inc_trend_scenario = INC_TREND_SCENARIO,
    delta_inc = DELTA_INC,
    sd_beta_I = SD_BETA_I,
    use_weighted_cohort = USE_WEIGHTED_COHORT,
    A_star = NA_real_, beta_mode = c("estimate","prior_ols","offset","fixed_rr_offset"), beta_force = NULL,
    gammaP_method = GAMMAP_METHOD, trend_type = TREND_TYPE, use_age_slope = FALSE,
    path_prev_dta = PATH_PREV_DTA,
    prev_micro_df = NULL,
    prev_cfg = NULL,
    bridge_inc_years = BRIDGE_INC_YEARS,
    L_I_max_years = L_I_MAX_YEARS,
    mort_period_shock_years = integer(0),
    mort_downweight_years_M = integer(0),
    mort_downweight_years_F = MORT_DOWNWEIGHT_YEARS_F,
    mort_downweight_weight_M = 1,
    mort_downweight_weight_F = MORT_DOWNWEIGHT_WEIGHT_F,
    mort_hist_tbl = mort_hist, pop_all_tbl = pop_all, inc_hist_tbl = NULL,
    emit_prev_diag_console = EMIT_PREV_DIAG_CONSOLE,
    emit_prev_diag_write = TRUE,
    cause_id_override = NA_character_,
    prev_inc_channel_mode = PREV_INC_CHANNEL_MODE,
    rr_inc = NA_real_,
    rr_mort = NA_real_
){
  if (is.null(prev_cfg)) prev_cfg <- make_prev_config()
  beta_mode <- match.arg(beta_mode, c("estimate","prior_ols","offset","fixed_rr_offset"))
  gammaP_method <- match.arg(gammaP_method)
  mort_trend_scenario <- match.arg(mort_trend_scenario)
  prev_inc_channel_mode <- "stock_former"
  
  # ---------- Estimación por sexo ----------
  # ---------- Estimación por sexo ----------
  resM <- run_pipeline_sex(
      sex_sel = "M",
      period_min_m = PERIOD_M_MIN, period_max_m = PERIOD_M_MAX,
      age_min_m = age_min_m, age_max_m = age_max_m,
      age_min_p = age_min_p, age_max_p = age_max_p,
      age_min_i = age_min_i, age_max_i = age_max_i,
      L_I = L_I, Da_I = Da_I, bridge_inc_years = bridge_inc_years,
      tech_scenario = mort_trend_scenario, delta_tech = delta_tech,
      inc_include_trend = inc_include_trend,
      inc_trend_scenario = inc_trend_scenario,
      delta_inc = delta_inc, sd_beta_I = sd_beta_I,
      use_weighted_cohort = use_weighted_cohort,
      anchor_pseudo_w = anchor_pseudo_w,
      sd_cohort_resid = sd_cohfix, sd_beta_fixed = sd_beta,
      use_age_slope = use_age_slope, A_star = A_star,
      beta_mode = beta_mode, beta_force = beta_force,
      gammaP_method = gammaP_method, trend_type = trend_type,
      path_prev_dta = path_prev_dta,
      prev_micro_df = prev_micro_df,
      prev_cfg = prev_cfg,
      L_I_max_years = L_I_max_years,
      mort_period_shock_years = mort_period_shock_years,
      mort_downweight_years = mort_downweight_years_M,
      mort_downweight_weight = mort_downweight_weight_M,
      mort_hist_tbl = mort_hist_tbl, pop_all_tbl = pop_all_tbl, inc_hist_tbl = inc_hist_tbl,
      cause_id_override = cause_id_override,
      rr_inc = rr_inc,
      rr_mort = rr_mort
    )
  
  resF <- run_pipeline_sex(
      sex_sel = "F",
      period_min_m = PERIOD_M_MIN, period_max_m = PERIOD_M_MAX,
      age_min_m = age_min_m, age_max_m = age_max_m,
      age_min_p = age_min_p, age_max_p = age_max_p,
      age_min_i = age_min_i, age_max_i = age_max_i,
      L_I = L_I, Da_I = Da_I, bridge_inc_years = bridge_inc_years,
      tech_scenario = mort_trend_scenario, delta_tech = delta_tech,
      inc_include_trend = inc_include_trend,
      inc_trend_scenario = inc_trend_scenario,
      delta_inc = delta_inc, sd_beta_I = sd_beta_I,
      use_weighted_cohort = use_weighted_cohort,
      anchor_pseudo_w = anchor_pseudo_w,
      sd_cohort_resid = sd_cohfix, sd_beta_fixed = sd_beta,
      use_age_slope = use_age_slope, A_star = A_star,
      beta_mode = beta_mode, beta_force = beta_force,
      gammaP_method = gammaP_method, trend_type = trend_type,
      path_prev_dta = path_prev_dta,
      prev_micro_df = prev_micro_df,
      prev_cfg = prev_cfg,
      L_I_max_years = L_I_max_years,
      mort_period_shock_years = mort_period_shock_years,
      mort_downweight_years = mort_downweight_years_F,
      mort_downweight_weight = mort_downweight_weight_F,
      mort_hist_tbl = mort_hist_tbl, pop_all_tbl = pop_all_tbl, inc_hist_tbl = inc_hist_tbl,
      cause_id_override = cause_id_override,
      rr_inc = rr_inc,
      rr_mort = rr_mort
    )
  
  if (is.null(resM) && is.null(resF)) stop("run_pipeline_both: no se pudo estimar ni M ni F.")
  # ---------- Helpers para joins seguros ----------
  add_missing_cols <- function(df, cols) {
    if (is.null(df)) return(NULL)
    for (nm in cols) if (!nm %in% names(df)) df[[nm]] <- NA_real_
    df
  }
  
  # ---------- annual_bapc (comb) ----------
  ab_m <- if (!is.null(resM) && !is.null(resM$annual_bapc) && nrow(resM$annual_bapc) > 0)
    dplyr::rename(resM$annual_bapc,
                  deaths_hat_M = deaths_hat, deaths_lwr_M = deaths_lwr, deaths_upr_M = deaths_upr) else NULL
  ab_f <- if (!is.null(resF) && !is.null(resF$annual_bapc) && nrow(resF$annual_bapc) > 0)
    dplyr::rename(resF$annual_bapc,
                  deaths_hat_F = deaths_hat, deaths_lwr_F = deaths_lwr, deaths_upr_F = deaths_upr) else NULL
  
  if (is.null(ab_m) && is.null(ab_f)) {
    comb <- tibble::tibble(period = integer(), deaths_hat = double(), deaths_lwr = double(), deaths_upr = double())
  } else if (is.null(ab_m)) {
    comb <- add_missing_cols(ab_f, c("deaths_hat_M","deaths_lwr_M","deaths_upr_M"))
  } else if (is.null(ab_f)) {
    comb <- add_missing_cols(ab_m, c("deaths_hat_F","deaths_lwr_F","deaths_upr_F"))
  } else {
    comb <- dplyr::full_join(ab_m, ab_f, by = "period")
  }
  
  comb <- comb %>%
    dplyr::mutate(
      deaths_hat = dplyr::coalesce(.data$deaths_hat_M, 0) + dplyr::coalesce(.data$deaths_hat_F, 0),
      deaths_lwr = dplyr::coalesce(.data$deaths_lwr_M, 0) + dplyr::coalesce(.data$deaths_lwr_F, 0),
      deaths_upr = dplyr::coalesce(.data$deaths_upr_M, 0) + dplyr::coalesce(.data$deaths_upr_F, 0)
    ) %>%
    dplyr::select(period, deaths_hat, deaths_lwr, deaths_upr) %>%
    dplyr::arrange(period)
  
  # ---------- annual_anchor (comb_anchor) ----------
  aa_m <- if (!is.null(resM) && !is.null(resM$annual_anchor) && nrow(resM$annual_anchor) > 0)
    dplyr::rename(resM$annual_anchor,
                  deaths_hat_M = deaths_hat, deaths_lwr_M = deaths_lwr, deaths_upr_M = deaths_upr) else NULL
  aa_f <- if (!is.null(resF) && !is.null(resF$annual_anchor) && nrow(resF$annual_anchor) > 0)
    dplyr::rename(resF$annual_anchor,
                  deaths_hat_F = deaths_hat, deaths_lwr_F = deaths_lwr, deaths_upr_F = deaths_upr) else NULL
  
  if (is.null(aa_m) && is.null(aa_f)) {
    comb_anchor <- tibble::tibble(period = integer(), deaths_hat = double(), deaths_lwr = double(), deaths_upr = double())
  } else if (is.null(aa_m)) {
    comb_anchor <- add_missing_cols(aa_f, c("deaths_hat_M","deaths_lwr_M","deaths_upr_M"))
  } else if (is.null(aa_f)) {
    comb_anchor <- add_missing_cols(aa_m, c("deaths_hat_F","deaths_lwr_F","deaths_upr_F"))
  } else {
    comb_anchor <- dplyr::full_join(aa_m, aa_f, by = "period")
  }
  
  comb_anchor <- comb_anchor %>%
    dplyr::mutate(
      deaths_hat = dplyr::coalesce(.data$deaths_hat_M, 0) + dplyr::coalesce(.data$deaths_hat_F, 0),
      deaths_lwr = dplyr::coalesce(.data$deaths_lwr_M, 0) + dplyr::coalesce(.data$deaths_lwr_F, 0),
      deaths_upr = dplyr::coalesce(.data$deaths_upr_M, 0) + dplyr::coalesce(.data$deaths_upr_F, 0)
    ) %>%
    dplyr::select(period, deaths_hat, deaths_lwr, deaths_upr) %>%
    dplyr::arrange(period)
  
  # ---------- annual_anchor_noP (comb_noP) ----------
  an_m <- if (!is.null(resM) && !is.null(resM$annual_anchor_noP) && nrow(resM$annual_anchor_noP) > 0)
    dplyr::rename(resM$annual_anchor_noP,
                  deaths_hat_M = deaths_hat, deaths_lwr_M = deaths_lwr, deaths_upr_M = deaths_upr) else NULL
  an_f <- if (!is.null(resF) && !is.null(resF$annual_anchor_noP) && nrow(resF$annual_anchor_noP) > 0)
    dplyr::rename(resF$annual_anchor_noP,
                  deaths_hat_F = deaths_hat, deaths_lwr_F = deaths_lwr, deaths_upr_F = deaths_upr) else NULL
  
  if (is.null(an_m) && is.null(an_f)) {
    comb_noP <- tibble::tibble(period = integer(), deaths_hat = double(), deaths_lwr = double(), deaths_upr = double())
  } else if (is.null(an_m)) {
    comb_noP <- add_missing_cols(an_f, c("deaths_hat_M","deaths_lwr_M","deaths_upr_M"))
  } else if (is.null(an_f)) {
    comb_noP <- add_missing_cols(an_m, c("deaths_hat_F","deaths_lwr_F","deaths_upr_F"))
  } else {
    comb_noP <- dplyr::full_join(an_m, an_f, by = "period")
  }
  
  comb_noP <- comb_noP %>%
    dplyr::mutate(
      deaths_hat = dplyr::coalesce(.data$deaths_hat_M, 0) + dplyr::coalesce(.data$deaths_hat_F, 0),
      deaths_lwr = dplyr::coalesce(.data$deaths_lwr_M, 0) + dplyr::coalesce(.data$deaths_lwr_F, 0),
      deaths_upr = dplyr::coalesce(.data$deaths_upr_M, 0) + dplyr::coalesce(.data$deaths_upr_F, 0)
    ) %>%
    dplyr::select(period, deaths_hat, deaths_lwr, deaths_upr) %>%
    dplyr::arrange(period)

  # ---------- obs_annual (obs_tot) ----------
  obs_m <- if (!is.null(resM) && !is.null(resM$obs_annual) && nrow(resM$obs_annual) > 0) resM$obs_annual else NULL
  obs_f <- if (!is.null(resF) && !is.null(resF$obs_annual) && nrow(resF$obs_annual) > 0) resF$obs_annual else NULL
  
  if (is.null(obs_m) && is.null(obs_f)) {
    obs_tot <- tibble::tibble(period = integer(), obs = double())
  } else if (is.null(obs_m)) {
    obs_tot <- obs_f %>% dplyr::select(period, obs) %>% dplyr::arrange(period)
  } else if (is.null(obs_f)) {
    obs_tot <- obs_m %>% dplyr::select(period, obs) %>% dplyr::arrange(period)
  } else {
    obs_tot <- dplyr::full_join(obs_m, obs_f, by = "period") %>%
      dplyr::mutate(obs = dplyr::coalesce(.data$obs.x, 0) + dplyr::coalesce(.data$obs.y, 0)) %>%
      dplyr::select(period, obs) %>% dplyr::arrange(period)
  }
  
  last_hist_year <- max(
    c(
      if (!is.null(resM) && !is.null(resM$diag$last_hist_year)) resM$diag$last_hist_year else -Inf,
      if (!is.null(resF) && !is.null(resF$diag$last_hist_year)) resF$diag$last_hist_year else -Inf
    ),
    na.rm = TRUE
  )

  horizon_tbl <- dplyr::bind_rows(
    tryCatch(resM$diag$projection_horizon_frontier, error = function(e) NULL),
    tryCatch(resF$diag$projection_horizon_frontier, error = function(e) NULL)
  )
  horizon_common_year <- projection_common_max_year_from_res_both(list(resM = resM, resF = resF), policy = "endogenous_max", default = NA_integer_)
  
  res_both <- list(
    params = list(
      anchor_pseudo_w = anchor_pseudo_w, sd_cohfix = sd_cohfix, sd_beta = sd_beta,
      age_min_m = age_min_m, age_max_m = age_max_m, age_min_p = age_min_p, age_max_p = age_max_p,
      age_min_i = age_min_i, age_max_i = age_max_i, L_I = L_I, Da_I = Da_I,
      mort_trend_scenario = mort_trend_scenario, delta_tech = delta_tech,
      inc_include_trend = inc_include_trend, inc_trend_scenario = inc_trend_scenario,
      delta_inc = delta_inc, sd_beta_I = sd_beta_I, use_weighted_cohort = use_weighted_cohort,
      A_star = NA_real_, beta_mode = beta_mode, gammaP_method = gammaP_method, trend_type = trend_type,
      method_policy_by_sex = list(
        M = tryCatch(resM$params$method_policy, error = function(e) NULL),
        F = tryCatch(resF$params$method_policy, error = function(e) NULL)
      ),
      method_policy_table_by_sex = dplyr::bind_rows(
        tryCatch(dplyr::mutate(resM$params$method_policy_table, sex = "M"), error = function(e) NULL),
        tryCatch(dplyr::mutate(resF$params$method_policy_table, sex = "F"), error = function(e) NULL)
      )
    ),
    resM = resM, resF = resF,
    combined = list(
      annual_bapc      = comb,
      annual_anchor    = comb_anchor,
      annual_anchor_noP = comb_noP,
      obs_annual       = obs_tot,
      last_hist_year   = last_hist_year,
      projection_horizon_frontier = horizon_tbl,
      max_projection_year_endogenous = horizon_common_year
    )
  )
  
  # ---------- Scores combinados (M + F) ----------
  scores_M <- if (!inherits(resM, "try-error") && !is.null(resM$fit_scores)) resM$fit_scores else NULL
  scores_F <- if (!inherits(resF, "try-error") && !is.null(resF$fit_scores)) resF$fit_scores else NULL
  res_both$fit_scores <- dplyr::bind_rows(scores_M, scores_F)
  
  # === exportar diagnóstico de PREV para este único escenario ===
  cause_label_rb <- tryCatch(unique(mort_hist_tbl$cause)[1], error = function(e) "unknown_cause")
  try(emit_diag_prev(res_both, cause_label = cause_label_rb, also_print = emit_prev_diag_console, write_csv = emit_prev_diag_write), silent = TRUE)
  
  return(res_both)
}



# =============================================================


run_pipeline_both_from_inputs <- function(inputs,
                                          cfg_row,
                                          prev_cfg = NULL,
                                          emit_prev_diag_console = EMIT_PREV_DIAG_CONSOLE,
                                          emit_prev_diag_write = TRUE,
                                          ...) {
  validate_bapc_inputs(inputs)
  cfg_row <- tibble::as_tibble(cfg_row)
  stopifnot(nrow(cfg_row) == 1)
  .bapc_verbose("inputs$mort_hist_tbl is NULL? ", is.null(inputs$mort_hist_tbl))
  run_pipeline_both(
    age_min_m = cfg_row$AGE_M_MIN[[1]], age_max_m = cfg_row$AGE_M_MAX[[1]],
    age_min_p = cfg_row$AGE_P_MIN[[1]], age_max_p = cfg_row$AGE_P_MAX[[1]],
    age_min_i = cfg_row$AGE_I_MIN[[1]], age_max_i = cfg_row$AGE_I_MAX[[1]],
    L_I = L_I_DEFAULT, Da_I = DA_I,
    bridge_inc_years = BRIDGE_INC_YEARS,
    prev_cfg = prev_cfg,
    L_I_max_years = if (exists(".extract_scalar", inherits = TRUE)) .extract_scalar(cfg_row$L_I_MAX_YEARS, L_I_MAX_YEARS) else cfg_row$L_I_MAX_YEARS[[1]],
    mort_period_shock_years = if (exists(".extract_intvec", inherits = TRUE)) .extract_intvec(cfg_row$MORT_SHOCK_YEARS) else integer(0),
    mort_downweight_years_F = if (exists(".extract_intvec", inherits = TRUE)) .extract_intvec(cfg_row$DOWNWEIGHT_F) else integer(0),
    mort_downweight_weight_F = MORT_DOWNWEIGHT_WEIGHT_F,
    mort_hist_tbl = inputs$mort_hist_tbl,
    pop_all_tbl   = inputs$pop_all_tbl,
    inc_hist_tbl  = inputs$inc_hist_tbl,
    path_prev_dta = if (!is.null(inputs$prev_path)) inputs$prev_path else PATH_PREV_DTA,
    prev_micro_df = if (is.data.frame(inputs$prev_data)) inputs$prev_data else NULL,
    A_star = NA_real_, beta_mode = BETA_MODE,
    gammaP_method = GAMMAP_METHOD, trend_type = TREND_TYPE, use_age_slope = FALSE,
    mort_trend_scenario = MORT_TREND_SCENARIO, delta_tech = DELTA_TECH,
    prev_inc_channel_mode = PREV_INC_CHANNEL_MODE,
    emit_prev_diag_console = emit_prev_diag_console,
    emit_prev_diag_write = emit_prev_diag_write,
    cause_id_override = if ("cause_id" %in% names(cfg_row)) as.character(cfg_row$cause_id[[1]]) else NA_character_,
    ...
  )
}

# freeze-benchmark rebuild helper for mortality anchor
# Extracted from the current simulation wrapper and formalized as a reusable helper.
# The benchmark residual comes from the base/freeze branch, and the scenario-specific
# mortality anchor is rebuilt by combining that residual with the scenario-specific offset.

.rebuild_anchor_freeze_benchmark <- function(res_base, res_scen, overwrite_main = TRUE) {
  `%||%` <- function(x, y) if (is.null(x)) y else x

  .rebuild_one_sex <- function(base_sex, scen_sex) {
    if (is.null(base_sex) || is.null(scen_sex)) return(list(base = base_sex, scen = scen_sex))
    req <- c("mort_anchor_pred_detail", "mort_anchor_data_cond", "annual_anchor")
    if (!all(req %in% names(base_sex)) || !all(req %in% names(scen_sex))) {
      return(list(base = base_sex, scen = scen_sex))
    }
    dat_b  <- tibble::as_tibble(base_sex$mort_anchor_data_cond)
    dat_s  <- tibble::as_tibble(scen_sex$mort_anchor_data_cond)
    pred_b <- tibble::as_tibble(base_sex$mort_anchor_pred_detail)
    pred_s <- tibble::as_tibble(scen_sex$mort_anchor_pred_detail)
    if (!nrow(dat_b) || !nrow(dat_s) || nrow(dat_b) != nrow(pred_b) || nrow(dat_s) != nrow(pred_s)) {
      return(list(base = base_sex, scen = scen_sex))
    }
    keys <- c("age","period","cohort")
    base_cells <- dplyr::bind_cols(
      dat_b |> dplyr::select(dplyr::all_of(keys), E = dplyr::any_of("E"), mort_ext_deaths = dplyr::any_of("mort_ext_deaths")),
      pred_b |> dplyr::select(hist_flag, offset_total, eta_resid_hat, mu_hat, mu_lwr = dplyr::any_of("mu_lwr"), mu_upr = dplyr::any_of("mu_upr"))
    ) |> dplyr::rename(base_offset_total = offset_total, base_eta_resid_hat = eta_resid_hat,
                       base_mu_hat = mu_hat, base_mu_lwr = mu_lwr, base_mu_upr = mu_upr)
    scen_cells <- dplyr::bind_cols(
      dat_s |> dplyr::select(dplyr::all_of(keys), E = dplyr::any_of("E"), mort_ext_deaths = dplyr::any_of("mort_ext_deaths")),
      pred_s |> dplyr::select(hist_flag, offset_total, eta_resid_hat, mu_hat, mu_lwr = dplyr::any_of("mu_lwr"), mu_upr = dplyr::any_of("mu_upr"))
    ) |> dplyr::rename(scen_offset_total = offset_total, scen_eta_resid_hat = eta_resid_hat,
                       scen_mu_hat = mu_hat, scen_mu_lwr = mu_lwr, scen_mu_upr = mu_upr)
    jj <- dplyr::full_join(
      base_cells |> dplyr::select(dplyr::all_of(keys), hist_flag, base_offset_total, base_eta_resid_hat, base_mu_hat, base_mu_lwr, base_mu_upr),
      scen_cells |> dplyr::select(dplyr::all_of(keys), hist_flag, scen_offset_total, scen_eta_resid_hat, scen_mu_hat, scen_mu_lwr, scen_mu_upr),
      by = c(keys, "hist_flag")
    )
    if (!nrow(jj)) return(list(base = base_sex, scen = scen_sex))
    rel_lwr <- .safe_num(jj$base_mu_lwr) / pmax(.safe_num(jj$base_mu_hat), 1e-12)
    rel_upr <- .safe_num(jj$base_mu_upr) / pmax(.safe_num(jj$base_mu_hat), 1e-12)
    rel_lwr[!is.finite(rel_lwr) | rel_lwr < 0] <- 1
    rel_upr[!is.finite(rel_upr) | rel_upr < 0] <- 1
    cand_mu_hat <- exp(.safe_num(jj$base_eta_resid_hat) + .safe_num(jj$scen_offset_total))
    cand_mu_lwr <- pmax(0, cand_mu_hat * rel_lwr)
    cand_mu_upr <- pmax(0, cand_mu_hat * rel_upr)
    pred_s_common <- pred_s |>
      dplyr::mutate(
        mu_hat = cand_mu_hat,
        mu_lwr = cand_mu_lwr,
        mu_upr = cand_mu_upr,
        anchor_variant = "commonresid_baseeta"
      )
    ann_s_common <- pred_s_common |>
      dplyr::group_by(period) |>
      dplyr::summarise(deaths_hat = sum(mu_hat, na.rm = TRUE),
                       deaths_lwr = sum(mu_lwr, na.rm = TRUE),
                       deaths_upr = sum(mu_upr, na.rm = TRUE), .groups = "drop")
    if (is.null(base_sex$annual_anchor_inla)) base_sex$annual_anchor_inla <- base_sex$annual_anchor
    if (is.null(base_sex$annual_anchor_raw_inla)) base_sex$annual_anchor_raw_inla <- base_sex$annual_anchor_raw %||% base_sex$annual_anchor
    if (is.null(base_sex$mort_anchor_pred_detail_inla)) base_sex$mort_anchor_pred_detail_inla <- base_sex$mort_anchor_pred_detail
    if (is.null(scen_sex$annual_anchor_inla)) scen_sex$annual_anchor_inla <- scen_sex$annual_anchor
    if (is.null(scen_sex$annual_anchor_raw_inla)) scen_sex$annual_anchor_raw_inla <- scen_sex$annual_anchor_raw %||% scen_sex$annual_anchor
    if (is.null(scen_sex$mort_anchor_pred_detail_inla)) scen_sex$mort_anchor_pred_detail_inla <- scen_sex$mort_anchor_pred_detail
    scen_sex$annual_anchor_commonresid_raw <- ann_s_common
    scen_sex$annual_anchor_commonresid <- ann_s_common
    scen_sex$mort_anchor_pred_detail_commonresid <- pred_s_common
    if (isTRUE(overwrite_main)) {
      scen_sex$annual_anchor_raw <- ann_s_common
      scen_sex$annual_anchor <- ann_s_common
      scen_sex$mort_anchor_pred_detail <- pred_s_common
    }
    list(base = base_sex, scen = scen_sex)
  }
  .rebuild_combined <- function(res_both) {
    if (is.null(res_both) || is.null(res_both$resM) || is.null(res_both$resF) || is.null(res_both$combined)) return(res_both)
    
    add_missing_cols <- function(df, cols) {
      if (is.null(df)) return(NULL)
      for (nm in cols) if (!nm %in% names(df)) df[[nm]] <- NA_real_
      df
    }

    # Helper for aggregating sex-specific tables
    .agg_sexes <- function(m_tbl, f_tbl) {
      if (is.null(m_tbl) && is.null(f_tbl)) return(NULL)
      m_rn <- if (!is.null(m_tbl)) dplyr::rename(m_tbl, deaths_hat_M = deaths_hat, deaths_lwr_M = deaths_lwr, deaths_upr_M = deaths_upr) else NULL
      f_rn <- if (!is.null(f_tbl)) dplyr::rename(f_tbl, deaths_hat_F = deaths_hat, deaths_lwr_F = deaths_lwr, deaths_upr_F = deaths_upr) else NULL
      
      if (is.null(m_rn)) {
        res <- add_missing_cols(f_rn, c("deaths_hat_M","deaths_lwr_M","deaths_upr_M"))
      } else if (is.null(f_rn)) {
        res <- add_missing_cols(m_rn, c("deaths_hat_F","deaths_lwr_F","deaths_upr_F"))
      } else {
        res <- dplyr::full_join(m_rn, f_rn, by = "period")
      }
      res %>%
        dplyr::mutate(
          deaths_hat = dplyr::coalesce(.data$deaths_hat_M, 0) + dplyr::coalesce(.data$deaths_hat_F, 0),
          deaths_lwr = dplyr::coalesce(.data$deaths_lwr_M, 0) + dplyr::coalesce(.data$deaths_lwr_F, 0),
          deaths_upr = dplyr::coalesce(.data$deaths_upr_M, 0) + dplyr::coalesce(.data$deaths_upr_F, 0)
        ) %>%
        dplyr::select(period, deaths_hat, deaths_lwr, deaths_upr) %>%
        dplyr::arrange(period)
    }

    # Update Informed SBAPC (always changes)
    res_both$combined$annual_anchor <- .agg_sexes(res_both$resM$annual_anchor, res_both$resF$annual_anchor)
    
    # Update Benchmarks (should be identical if sex-specifics are identical)
    res_both$combined$annual_bapc <- .agg_sexes(res_both$resM$annual_bapc, res_both$resF$annual_bapc)
    res_both$combined$annual_anchor_noP <- .agg_sexes(res_both$resM$annual_anchor_noP, res_both$resF$annual_anchor_noP)
    
    res_both
  }
  mm <- .rebuild_one_sex(res_base$resM, res_scen$resM)
  ff <- .rebuild_one_sex(res_base$resF, res_scen$resF)
  res_base$resM <- mm$base; res_scen$resM <- mm$scen
  res_base$resF <- ff$base; res_scen$resF <- ff$scen
  res_base <- .rebuild_combined(res_base)
  res_scen <- .rebuild_combined(res_scen)
  list(res_base = res_base, res_scen = res_scen)
}

# freeze-benchmark rebuild helper for scenario branches in simulation
# Reconstructs scenario-specific incidence from the freeze/base branch,
# then rebuilds mortality from that scenario-specific incidence and the
# freeze residual benchmark.

.rebuild_incidence_freeze_benchmark <- function(res_base,
                                                inputs,
                                                cfg_row,
                                                prev_cfg_scen,
                                                overwrite_main = TRUE) {
  `%||%` <- function(x, y) if (is.null(x)) y else x


  .norm_sex <- function(x) toupper(substr(as.character(x)[1], 1, 1))

  .annualise_inc <- function(rates_all, pop_tbl, sex_sel, last_hist) {
    if (is.null(rates_all) || !nrow(rates_all)) {
      return(tibble::tibble(period = integer(), cases_hat = numeric(), rate_hat = numeric(), rate_lwr = numeric(), rate_upr = numeric()))
    }
    pop_sex <- pop_tbl %>%
      dplyr::mutate(sex = as.character(sex)) %>%
      dplyr::filter(sex == sex_sel) %>%
      dplyr::select(sex, age, period, exposure)
    rates_all %>%
      dplyr::mutate(sex = as.character(sex), age = as.integer(age), period = as.integer(period)) %>%
      dplyr::left_join(pop_sex, by = c("sex", "age", "period")) %>%
      dplyr::filter(period > last_hist) %>%
      dplyr::mutate(exposure = pmax(.safe_num(exposure), 1e-12),
                    rate_hat = pmax(.safe_num(rate_hat), 1e-12),
                    rate_lwr = pmax(.safe_num(rate_lwr), 1e-12),
                    rate_upr = pmax(.safe_num(rate_upr), 1e-12),
                    cases_hat = exposure * rate_hat,
                    cases_lwr = exposure * rate_lwr,
                    cases_upr = exposure * rate_upr) %>%
      dplyr::group_by(period) %>%
      dplyr::summarise(cases_hat = sum(cases_hat, na.rm = TRUE),
                       cases_lwr = sum(cases_lwr, na.rm = TRUE),
                       cases_upr = sum(cases_upr, na.rm = TRUE), .groups = "drop")
  }

  .prep_mort_template_attach <- function(df) {
    tibble::as_tibble(df) %>%
      dplyr::select(-dplyr::any_of(c(
        "E", "logE",
        "mort_ext_deaths", "mort_ext_rate",
        "log_mort_ext", "log_mort_ext_rate"
      )))
  }

  cause_id_cur <- tryCatch(as.character(cfg_row$cause_id[[1]]), error = function(e) NA_character_)
  pop_all_tbl <- inputs$pop_all %||% inputs$pop_all_tbl
  if (is.null(pop_all_tbl) || !is.data.frame(pop_all_tbl) || !nrow(pop_all_tbl)) {
    stop("No encontré pop_all en inputs para reconstruir incidencia freeze-benchmark.")
  }

  # Usamos asignación simple; R maneja el copy-on-write eficientemente sin desbordar el buffer de serialización
  res_scen <- res_base

  .rebuild_one_sex <- function(base_sex, scen_sex, sex_sel) {
    if (is.null(base_sex) || is.null(scen_sex)) return(list(base = base_sex, scen = scen_sex))

    sex_sel <- .norm_sex(sex_sel)
    if (is.null(base_sex$inc_fit$rates_all_full) || !nrow(base_sex$inc_fit$rates_all_full)) {
      return(list(base = base_sex, scen = scen_sex))
    }

    base_full <- tibble::as_tibble(base_sex$inc_fit$rates_all_full) %>%
      dplyr::mutate(sex = as.character(sex), age = as.integer(age), period = as.integer(period))
    last_hist <- as.integer(base_sex$inc_fit$last_year_inc %||% base_sex$diag$last_hist_year %||% 2022L)

    if (is.null(base_sex$inc_annual_cond_inla)) base_sex$inc_annual_cond_inla <- base_sex$inc_annual_cond
    if (is.null(scen_sex$inc_annual_cond_inla)) scen_sex$inc_annual_cond_inla <- scen_sex$inc_annual_cond
    if (is.null(base_sex$inc_fit_inla)) base_sex$inc_fit_inla <- base_sex$inc_fit
    if (is.null(scen_sex$inc_fit_inla)) scen_sex$inc_fit_inla <- scen_sex$inc_fit

    if (identical(normalize_prev_scenario_name(prev_cfg_scen$scenario %||% "freeze"), "freeze")) {
      scen_sex$diag$scenario_build_mode <- "freeze_copy"
      scen_sex$diag$scenario_build_benchmark <- "freeze"
      return(list(base = base_sex, scen = scen_sex))
    }

    rr_use <- base_sex$params$rr_inc
    if (is.null(rr_use) || !is.finite(rr_use)) {
      rr_use <- tryCatch(get_inc_rr_by_cause_sex(cause_id_cur, sex_sel), error = function(e) SIM_RR_I_DEFAULT)
    }
    rr_use <- .safe_num(rr_use)[1]
    if (!is.finite(rr_use) || rr_use <= 1) rr_use <- SIM_RR_I_DEFAULT

    fut <- base_full %>% dplyr::filter(period > last_hist)
    if (!nrow(fut)) {
      scen_sex$diag$scenario_build_mode <- "freeze_copy_no_future"
      scen_sex$diag$scenario_build_benchmark <- "freeze"
      return(list(base = base_sex, scen = scen_sex))
    }

    # Rebuild the smoking-exposure channel with the same stock-former helper used by
    # the main PREV->INC engine. Do not apply the legacy q_eff correction shortcut here.
    # CRITICAL: We MUST pass the full grid (history + future) to the stock builder 
    # to preserve the memory of former smokers across the 2022/2023 boundary.
    full_stock_grid <- base_full %>%
      dplyr::select(dplyr::any_of(c(
        "sex", "age", "period", "cohort",
        "age_id", "period_id", "cohort_id",
        "E", "logE", "exposure"
      )))

    full_stock_surface <- build_prev_rr_offset_stock_for_inc(
      df_inc_grid = full_stock_grid,
      fit_prev = base_sex$fit_prev,
      cause_id = cause_id_cur,
      rr_inc = rr_use,
      prev_inla = tryCatch(base_sex$fit_prev$.args$data, error = function(e) NULL),
      sex_sel = sex_sel,
      gammaP_method = GAMMAP_METHOD,
      trend_type = TREND_TYPE,
      prev_cfg = prev_cfg_scen,
      age_min_p = AGE_P_MIN,
      age_max_p = AGE_P_MAX,
      backcast_period_mode = PREV_BACKCAST_MODE,
      backcast_cohort_mode = PREV_BACKCAST_COHORT_MODE,
      post65_mode = PREV_POST65_MODE,
      quit_horizon_years = PREV_INC_MAX_QUIT_YEARS,
      return_internal = FALSE
    )

    stock_surface <- full_stock_surface %>% dplyr::filter(period > last_hist)
    join_keys_scen <- c("sex", "age", "period", "cohort")
    stock_cols <- c(
      "q_eff", "z_prev", "p_cur", "delta_p_cur", "quit_flow",
      "p_never", "p_former_total", "noncurrent_rescale",
      "offset_prev_rr", "rr_inc", "quit_horizon_years",
      "prev_source", "within_prev_age_support", "within_prev_period_support",
      "within_prev_observed_support", "within_prev_support",
      "prev_scenario_name", "prev_scenario_applied",
      "offset_prev_rr", "rr_inc", "quit_horizon_years"
    )

    # Scoped Stock Surface
    scen_stock_to_join <- stock_surface[, intersect(names(stock_surface), c(join_keys_scen, stock_cols)), drop = FALSE]
    
    # Ensure join keys are same type (Character)
    fut$sex <- as.character(fut$sex)
    scen_stock_to_join$sex <- as.character(scen_stock_to_join$sex)
    
    # Universal Merge (Base R is safer here for complex joins)
    fut_merged <- merge(
      fut, 
      scen_stock_to_join, 
      by = join_keys_scen, 
      all.x = TRUE, 
      suffixes = c("", "_scen")
    )
    
    fut <- tibble::as_tibble(fut_merged)

    # Diagnostic check for join integrity
    n_miss <- sum(is.na(fut$offset_prev_rr_scen))
    if (n_miss > 0 && n_miss < nrow(fut)) {
       .bapc_verbose(sprintf("[REBUILD] Warning: %d/%d rows in %s scenario didn't match rebuilt stock. Coalescing.", n_miss, nrow(fut), scen))
    }

    # Extract sensitivity coefficient (theta). 
    bz_hat <- get_incidence_sensitivity_coef(
      beta_mode = base_sex$meta$beta_mode, 
      beta_P_eff = base_sex$inc_fit$beta_P_eff
    )
    
    fut <- fut %>%
      dplyr::mutate(
        off_epi_base = .safe_num(dplyr::coalesce(offset_prev_rr, 0)),
        off_epi_scen = .safe_num(dplyr::coalesce(offset_prev_rr_scen, offset_prev_rr, off_epi_base)),
        
        # Scenario Incidence Shift: Apply Delta-Offset logic (M-2.2)
        rate_scen = apply_incidence_scenario_shift(rate_hat, off_epi_scen, off_epi_base, bz_hat),
        
        # Preserve uncertainty relatives
        rel_lwr = pmax(rate_lwr / pmax(rate_hat, 1e-12), 1e-12),
        rel_upr = pmax(rate_upr / pmax(rate_hat, 1e-12), 1e-12),
        rate_lwr_scen = pmax(rate_scen * rel_lwr, 1e-12),
        rate_upr_scen = pmax(rate_scen * rel_upr, 1e-12),
        
        # Update Case Counts
        mu_scen = rate_scen * E,
        mu_lwr_scen = rate_lwr_scen * E,
        mu_upr_scen = rate_upr_scen * E
      )

    # Metadata updates
    fut$offset_prev_rr <- fut$off_epi_scen
    fut$z_prev <- fut$off_epi_scen # Synchronize z_prev for diagnostics
    
    # Update other stock columns from the joined surface
    stock_meta_cols <- setdiff(stock_cols, c("offset_prev_rr", "z_prev"))
    for (col in stock_meta_cols) {
      col_scen <- paste0(col, "_scen")
      if (col_scen %in% names(fut)) {
        fut[[col]] <- fut[[col_scen]]
      }
    }

    rate_scen <- fut$rate_scen
    rate_lwr_scen <- fut$rate_lwr_scen
    rate_upr_scen <- fut$rate_upr_scen
    mu_scen <- fut$mu_scen
    mu_lwr_scen <- fut$mu_lwr_scen
    mu_upr_scen <- fut$mu_upr_scen
    
    off_total_scen <- fut$off_epi_scen + (.safe_num(fut$coef_fc_offset_I) - fut$off_epi_base)
    
    eta_offset_scen <- dplyr::coalesce(.safe_num(fut$inc_tech_offset), 0) + off_total_scen
    eta_apc_base <- .safe_num(fut$eta_apc_manual)
    eta_total_scen <- ifelse(is.finite(eta_apc_base), eta_apc_base + eta_offset_scen, log(pmax(rate_scen, 1e-12)))

    fut_rebuilt <- fut %>%
      dplyr::mutate(
        coef_fc_offset_I_epi_new = off_epi_scen,
        coef_fc_offset_I_new = (off_epi_scen + (.safe_num(coef_fc_offset_I) - off_epi_base)),
        coef_fc_offset_I_raw_new = (off_epi_scen + (.safe_num(coef_fc_offset_I) - off_epi_base)),
        coef_fc_offset_I_effective_new = (off_epi_scen + (.safe_num(coef_fc_offset_I) - off_epi_base)),
        coef_fc_posthoc_adj_new = 1,
        coef_fc_posthoc_lock_mode_new = "freeze_benchmark",
        rate_hat_new = rate_scen,
        rate_lwr_new = rate_lwr_scen,
        rate_upr_new = rate_upr_scen,
        mu_hat_new = mu_scen,
        mu_lwr_new = mu_lwr_scen,
        mu_upr_new = mu_upr_scen,
        z_prev_new = off_epi_scen,
        offset_prev_rr_new = off_epi_scen
      ) %>%
      dplyr::select(dplyr::any_of(join_keys_scen), dplyr::ends_with("_new"))

    # Diagnostic columns to preserve/update
    stock_diag_cols <- c("p_cur", "q_eff", "z_prev", "delta_p_cur", "quit_flow", 
                         "p_never", "p_former_total", "noncurrent_rescale", 
                         "offset_prev_rr", "rr_inc", "quit_horizon_years")

    join_keys <- intersect(c("sex", "age", "period", "cohort"), names(base_full))
    scen_full <- base_full %>%
      dplyr::left_join(
        fut_rebuilt,
        by = join_keys
      )
    
    # Apply new values where they exist (future)
    update_cols <- intersect(names(scen_full), paste0(stock_diag_cols, "_new"))
    for (uc in update_cols) {
      target <- gsub("_new$", "", uc)
      scen_full[[target]] <- dplyr::coalesce(.safe_num(scen_full[[uc]]), .safe_num(scen_full[[target]]))
    }
    
    # Update rates and case counts
    # Update rates and case counts safely
    for (nm in c("rate_hat", "rate_lwr", "rate_upr", "mu_hat", "mu_lwr", "mu_upr")) {
      new_nm <- paste0(nm, "_new")
      if (new_nm %in% names(scen_full)) {
        scen_full[[nm]] <- dplyr::coalesce(.safe_num(scen_full[[new_nm]]), .safe_num(scen_full[[nm]]))
        scen_full[[new_nm]] <- NULL
      }
    }

    # If q_eff was absent historically, keep it absent; otherwise future gets rebuilt values.
    if ("q_eff" %in% names(scen_full)) {
      n_sf <- nrow(scen_full)
      q_eff_new_vec <- if ("q_eff_new" %in% names(scen_full)) .safe_num(scen_full$q_eff_new) else rep(NA_real_, n_sf)
      q_eff_old_vec <- if ("q_eff" %in% names(scen_full)) .safe_num(scen_full$q_eff) else rep(NA_real_, n_sf)
      scen_full$q_eff <- dplyr::coalesce(q_eff_new_vec, q_eff_old_vec)
    }
    
    # Rebuild logic:
    # 1. We take the Pure BAPC Incidence from the base object (inertial trend, no scenarios)
    # 2. We take the Scenario Incidence from the current fit (includes prevalence impact)
    # 3. We rebuild mortality using these two distinct info levels.
    
    inc_bapc_rates_all <- base_sex$inc_fit_bapc$rates_all_full %>%
      dplyr::select(sex, age, period, rate_hat)
    
    # Validation: Ensure BAPC incidence is available
    if (is.null(inc_bapc_rates_all) || nrow(inc_bapc_rates_all) == 0) {
      .bapc_verbose("[.rebuild_one_sex] Warning: Pure BAPC incidence not found. Falling back to freeze incidence for Orange line.")
      inc_bapc_rates_all <- base_full %>% dplyr::select(sex, age, period, rate_hat)
    }

    replace_cols <- c(
      "z_prev","p_cur","delta_p_cur","quit_flow","p_never","p_former_total","noncurrent_rescale",
      "offset_prev_rr","rr_inc","quit_horizon_years","prev_source",
      "within_prev_age_support","within_prev_period_support","within_prev_observed_support","within_prev_support",
      "prev_scenario_name","prev_scenario_applied",
      "coef_fc_offset_I_epi","coef_fc_offset_I_apc","coef_fc_offset_I",
      "coef_fc_offset_I_raw","coef_fc_offset_I_effective","coef_fc_posthoc_adj","coef_fc_posthoc_lock_mode",
      "eta_offset_manual","eta_total_manual","rate_manual","mu_manual","lp_gap_manual",
      "lp_mean","lp_lwr","lp_upr","fv_mean","fv_lwr","fv_upr",
      "mu_hat","mu_lwr","mu_upr","rate_from_fv_over_E","rate_from_fv_over_E_times_offset",
      "rate_from_lp_over_E","rate_from_lp_over_E_times_offset","rate_blend_geom","rate_blend_geom_times_offset",
      "rate_blend_arith","rate_blend_logmid","rate_hat","rate_lwr","rate_upr"
    )
    
    for (nm in replace_cols) {
      new_nm <- paste0(nm, "_new")
      if (new_nm %in% names(scen_full)) {
        if (nm %in% names(scen_full)) {
          scen_full[[nm]] <- dplyr::coalesce(scen_full[[new_nm]], scen_full[[nm]])
        } else {
          scen_full[[nm]] <- scen_full[[new_nm]]
        }
        scen_full[[new_nm]] <- NULL
      }
    }
    if ("q_eff_new" %in% names(scen_full)) scen_full$q_eff_new <- NULL

    scen_rates <- scen_full %>% dplyr::select(sex, age, period, rate_hat, rate_lwr, rate_upr)
    pop_sex <- pop_all_tbl %>% dplyr::mutate(sex = as.character(sex)) %>% dplyr::filter(sex == sex_sel)
    inc_annual_cond <- .annualise_inc(scen_rates, pop_sex, sex_sel = sex_sel, last_hist = last_hist)

    zf_rebuilt <- fut
    for (col in stock_diag_cols) {
      col_scen <- paste0(col, "_scen")
      if (col_scen %in% names(zf_rebuilt)) {
        # target_col already exists in zf_rebuilt (from the base object)
        zf_rebuilt[[col]] <- dplyr::coalesce(.safe_num(zf_rebuilt[[col_scen]]), .safe_num(zf_rebuilt[[col]]))
      }
    }
    
    zf_rebuilt <- zf_rebuilt %>%
      dplyr::select(dplyr::any_of(c("sex", "age", "period", "cohort", "q_eff", "z_prev",
                                     "p_cur", "delta_p_cur", "quit_flow", "p_never", "p_former_total",
                                     "noncurrent_rescale", "offset_prev_rr", "rr_inc", "quit_horizon_years",
                                     "prev_source", "within_prev_age_support", "within_prev_period_support",
                                     "within_prev_observed_support", "within_prev_support",
                                     "prev_scenario_name", "prev_scenario_applied",
                                     "coef_fc_signal_I", "coef_fc_recenter_I",
                                     "coef_fc_offset_I_epi", "coef_fc_offset_I_apc", "coef_fc_offset_I",
                                     "period_raw", "mapped_period", "period_is_clamped", "period_shift",
                                     "cohort_raw", "mapped_cohort", "cohort_is_edge", "cohort_clamp_low", "cohort_clamp_high", "cohort_shift",
                                     "support_n", "support_frac", "horizon", "horizon_block")))

    scen_sex$inc_fit$rates_all <- scen_rates
    scen_sex$inc_fit$rates_all_full <- scen_full
    scen_sex$inc_annual_cond <- inc_annual_cond
    scen_sex$diag$z_prev_future <- zf_rebuilt
    scen_sex$diag$scenario_build_mode <- "freeze_benchmark_rebuild"
    scen_sex$diag$scenario_build_benchmark <- "freeze"
    scen_sex$diag$scenario_name <- as.character(prev_cfg_scen$scenario %||% NA_character_)

    # 2. Get scenario offsets (Full and NoP)
    # Mortality Offset Calculation
    # mort_data_cond: Uses Scenario Incidence (I|P)
    # mort_data_noP: Uses Pure BAPC Incidence (I)
    
    # Setup join keys for mortality
    key_candidates_m <- c("sex", "age", "period", "cohort")
    join_keys_m <- intersect(key_candidates_m, names(scen_full))

    # Construct mortality template: Start with tech offsets from base, add prev from scen
    mort_template_attach <- .prep_mort_template_attach(base_sex$mort_anchor_data_cond) %>%
      dplyr::left_join(
        scen_full %>% dplyr::select(dplyr::all_of(join_keys_m), mort_offset_epi = offset_prev_rr),
        by = join_keys_m
      )

    mort_data_cond <- attach_external_mortality_offset(
      mort_all = mort_template_attach,
      inc_rates_all = scen_rates %>% dplyr::select(sex, age, period, rate_hat),
      pop_all_tbl = pop_sex,
      cause_id = cause_id_cur,
      sex_sel = sex_sel
    )
    
    mort_data_noP <- attach_external_mortality_offset(
      mort_all = mort_template_attach,
      inc_rates_all = inc_bapc_rates_all,
      pop_all_tbl = pop_sex,
      cause_id = cause_id_cur,
      sex_sel = sex_sel
    )

    # Recover keys in pred_base and mort_data_cond if missing
    pred_base <- recover_demographic_keys(pred_base)
    mort_data_cond <- recover_demographic_keys(mort_data_cond)

    key_candidates_m <- c("sex", "age", "period", "cohort")
    join_keys_m <- intersect(key_candidates_m, intersect(names(mort_data_cond), names(pred_base)))
