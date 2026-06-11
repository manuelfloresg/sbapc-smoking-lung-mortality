# 4) Modelos APC por etapa
# =============================================================

# --- Anclaje de nivel histórico para incidencia (ajuste de intercepto en LP) ---
apply_incidence_level_anchor <- function(grid_all, last_hist_year, period_col = "period", y_col = "cases", E_col = "E") {
  if (!is.data.frame(grid_all) || nrow(grid_all) == 0) return(grid_all)
  if (!(period_col %in% names(grid_all)) || !(y_col %in% names(grid_all)) || !(E_col %in% names(grid_all))) {
    return(grid_all)
  }
  req_lp <- c("lp_mean", "lp_lwr", "lp_upr")
  if (!all(req_lp %in% names(grid_all))) return(grid_all)

  hist_idx <- is.finite(grid_all[[period_col]]) & (grid_all[[period_col]] <= last_hist_year) &
    is.finite(grid_all[[y_col]]) & !is.na(grid_all[[y_col]])
  if (!any(hist_idx)) {
    grid_all$inc_level_anchor_adj <- 1
    grid_all$inc_level_anchor_log_adj <- 0
    grid_all$inc_level_anchor_basis <- "no_hist_rows"
    return(grid_all)
  }

  obs_hist <- sum(as.numeric(grid_all[[y_col]][hist_idx]), na.rm = TRUE)
  pred_hist <- sum(exp(as.numeric(grid_all$lp_mean[hist_idx])), na.rm = TRUE)

  if (!is.finite(obs_hist) || !is.finite(pred_hist) || pred_hist <= 0 || obs_hist <= 0) {
    grid_all$inc_level_anchor_adj <- 1
    grid_all$inc_level_anchor_log_adj <- 0
    grid_all$inc_level_anchor_basis <- "invalid_hist_totals"
    return(grid_all)
  }

  adj <- obs_hist / pred_hist
  log_adj <- log(adj)
  E_safe <- pmax(as.numeric(grid_all[[E_col]]), 1e-12)

  grid_all$inc_level_anchor_adj <- adj
  grid_all$inc_level_anchor_log_adj <- log_adj
  grid_all$inc_level_anchor_basis <- "hist_cases_over_hist_pred"

  grid_all$lp_mean <- as.numeric(grid_all$lp_mean) + log_adj
  grid_all$lp_lwr  <- as.numeric(grid_all$lp_lwr) + log_adj
  grid_all$lp_upr  <- as.numeric(grid_all$lp_upr) + log_adj

  if ("eta_apc_manual" %in% names(grid_all)) grid_all$eta_apc_manual <- as.numeric(grid_all$eta_apc_manual) + log_adj
  if ("eta_total_manual" %in% names(grid_all)) grid_all$eta_total_manual <- as.numeric(grid_all$eta_total_manual) + log_adj

  lp_rate_mean <- exp(as.numeric(grid_all$lp_mean)) / E_safe
  lp_rate_lwr  <- exp(as.numeric(grid_all$lp_lwr)) / E_safe
  lp_rate_upr  <- exp(as.numeric(grid_all$lp_upr)) / E_safe

  if ("rate_from_lp_over_E" %in% names(grid_all)) grid_all$rate_from_lp_over_E <- lp_rate_mean
  if ("rate_from_lp_over_E_times_offset" %in% names(grid_all)) {
    if ("eta_total_manual" %in% names(grid_all)) {
      grid_all$rate_from_lp_over_E_times_offset <- exp(as.numeric(grid_all$eta_total_manual))
    } else {
      grid_all$rate_from_lp_over_E_times_offset <- lp_rate_mean
    }
  }

  if (all(c("rate_from_fv_over_E","rate_from_lp_over_E") %in% names(grid_all))) {
    grid_all$rate_blend_geom <- sqrt(pmax(as.numeric(grid_all$rate_from_fv_over_E), 0) * pmax(as.numeric(grid_all$rate_from_lp_over_E), 0))
    grid_all$rate_blend_arith <- 0.5 * (as.numeric(grid_all$rate_from_fv_over_E) + as.numeric(grid_all$rate_from_lp_over_E))
    grid_all$rate_blend_logmid <- exp(0.5 * (log(pmax(as.numeric(grid_all$rate_from_fv_over_E), 1e-300)) + log(pmax(as.numeric(grid_all$rate_from_lp_over_E), 1e-300))))
  }
  if ("rate_blend_geom" %in% names(grid_all)) {
    if ("eta_total_manual" %in% names(grid_all)) {
      grid_all$rate_blend_geom_times_offset <- sqrt(pmax(as.numeric(grid_all$rate_from_fv_over_E_times_offset), 0) * pmax(exp(as.numeric(grid_all$eta_total_manual)), 0))
    } else {
      grid_all$rate_blend_geom_times_offset <- grid_all$rate_blend_geom
    }
  }

  if (all(c("lp_mean","eta_total_manual") %in% names(grid_all))) {
    manual_rate <- exp(as.numeric(grid_all$eta_total_manual))
    adj_center <- manual_rate / pmax(lp_rate_mean, 1e-12)

    if ("rate_manual" %in% names(grid_all)) grid_all$rate_manual <- manual_rate
    if ("mu_manual" %in% names(grid_all)) grid_all$mu_manual <- manual_rate * E_safe
    if ("mu_hat" %in% names(grid_all)) grid_all$mu_hat <- manual_rate * E_safe
    if ("mu_lwr" %in% names(grid_all)) grid_all$mu_lwr <- lp_rate_lwr * adj_center * E_safe
    if ("mu_upr" %in% names(grid_all)) grid_all$mu_upr <- lp_rate_upr * adj_center * E_safe
    grid_all$rate_hat <- manual_rate
    grid_all$rate_lwr <- lp_rate_lwr * adj_center
    grid_all$rate_upr <- lp_rate_upr * adj_center
    grid_all$lp_gap_manual <- as.numeric(grid_all$lp_mean) - (log(E_safe) + as.numeric(grid_all$eta_total_manual))
  } else {
    if ("mu_hat" %in% names(grid_all)) grid_all$mu_hat <- exp(as.numeric(grid_all$lp_mean))
    if ("mu_lwr" %in% names(grid_all)) grid_all$mu_lwr <- exp(as.numeric(grid_all$lp_lwr))
    if ("mu_upr" %in% names(grid_all)) grid_all$mu_upr <- exp(as.numeric(grid_all$lp_upr))
    if ("rate_manual" %in% names(grid_all)) grid_all$rate_manual <- lp_rate_mean
    if ("mu_manual" %in% names(grid_all)) grid_all$mu_manual <- as.numeric(grid_all$rate_manual) * E_safe
    if (all(c("lp_mean","eta_total_manual") %in% names(grid_all))) {
      grid_all$lp_gap_manual <- as.numeric(grid_all$lp_mean) - (log(E_safe) + as.numeric(grid_all$eta_total_manual))
    }
    grid_all$rate_hat <- lp_rate_mean
    grid_all$rate_lwr <- lp_rate_lwr
    grid_all$rate_upr <- lp_rate_upr
  }
  grid_all
}

# ---------------------------
# 4A) Incidencia APC simple (Poisson) con anti-esquinas
# ---------------------------
fit_apc_incidence <- function(inc_hist, pop_all,
                              age_min_i = AGE_I_MIN, age_max_i = AGE_I_MAX,
                              period_min_i = PERIOD_M_MIN, period_max_i = PERIOD_M_MAX,
                              use_weighted_cohort = USE_WEIGHTED_COHORT,
                              inc_scenario = INC_TREND_SCENARIO, delta_inc = DELTA_INC,
                              sd_beta_I   = SD_BETA_I,
                              inc_degree   = INC_TREND_DEGREE) {
  if (identical(INC_TREND_ON, "none")) inc_degree <- 0L
  pop_all <- ensure_exposure(pop_all)
  
  incH <- inc_hist %>%
    dplyr::filter(age >= age_min_i, age <= age_max_i,
                  period >= period_min_i, period <= period_max_i) %>%
    dplyr::arrange(sex, period, age) %>%
    dplyr::left_join(pop_all, by = c("period","age","sex")) %>%
    ensure_exposure() %>%
    dplyr::mutate(cohort = period - age,
                  y = as.integer(round(cases)),
                  E = pmax(as.numeric(exposure), 1e-12),
                  logE = log(E))
  
  if (!exists("incH", inherits = FALSE) || is.null(incH) || nrow(incH) == 0 || !any(is.finite(incH$period))) {
    stop("[fit_apc_incidence] empty historical incidence for this sex (n=0). Check filters/windows or whether the cause belongs only to the other sex.")
  }
  
  last_year_inc <- max(incH$period)
  lev_age <- sort(unique(incH$age)); lev_per <- sort(unique(incH$period)); lev_coh <- sort(unique(incH$cohort))
  mu_perI <- mean(incH$period)
  Xconstr_per_I <- make_slope_constr(lev_per, mu_perI)
  
  inc_inla <- incH %>%
    dplyr::mutate(age_id = match(age, lev_age),
                  period_id = match(period, lev_per),
                  cohort_id = match(cohort, lev_coh)) %>%
    dplyr::bind_cols(make_trend_vars(.$period, mu_perI, inc_degree, prefix = "inc_")) %>%
    dplyr::mutate(inc_tech_offset = 0, coef_fc_signal_I = 0, coef_fc_recenter_I = 0, coef_fc_offset_I = 0)
  
  inc_inla <- inc_inla %>%
    dplyr::mutate(period_iid = if (isTRUE(INC_PER_EXTRA_IID)) period_id else NA_integer_)
  
  inc_inla <- attach_edge_weights_hist(inc_inla, stage = "inc")
  
  hyper_age <- pc_hyper(INC_AGE_PC_U, INC_AGE_PC_A)
  hyper_per <- pc_hyper(INC_PER_PC_U, INC_PER_PC_A)
  hyper_coh <- pc_hyper(INC_COH_PC_U, INC_COH_PC_A)
  
  # --- Fórmula condicional según la perilla INC_TREND_ON ---
  if (INC_TREND_ON == "period") {
    .bapc_verbose("[fit_apc_incidence] Assigning linear trend to PERIOD.")
    base_formula <- y ~ 1 + inc_trend_t +
      f(age_id,    model = INC_AGE_MODEL, constr=TRUE, scale.model=TRUE, hyper=hyper_age) +
      f(period_id, model = INC_PER_MODEL, constr=TRUE, scale.model=TRUE, hyper=hyper_per)
    
    ctrl_fix <- list(
      mean = list(inc_trend_t = 0),
      prec = list(inc_trend_t = 1 / sd_beta_I^2)
    )
    
  } else if (INC_TREND_ON == "cohort") {
    .bapc_verbose("[fit_apc_incidence] Assigning trend to COHORT (rw2).")
    base_formula <- y ~ 1 +
      f(age_id,    model = INC_AGE_MODEL, constr=TRUE, scale.model=TRUE, hyper=hyper_age) +
      f(period_id, model = INC_PER_MODEL, constr=TRUE, scale.model=TRUE, hyper=hyper_per)
    
    ctrl_fix <- list()
    use_weighted_cohort <- FALSE  # rw2 de cohorte no es compatible con la matriz de pesos
    
  } else if (INC_TREND_ON == "none") {
    .bapc_verbose("[fit_apc_incidence] No explicit trend (absorbed by APC components).")
    base_formula <- y ~ 1 +
      f(age_id,    model = INC_AGE_MODEL, constr=TRUE, scale.model=TRUE, hyper=hyper_age) +
      f(period_id, model = INC_PER_MODEL, constr=TRUE, scale.model=TRUE, hyper=hyper_per)
    
    ctrl_fix <- list()
    
  } else {
    stop("INC_TREND_ON must be 'period', 'cohort' or 'none'.")
  }
  
  # --- Añadir el efecto de cohorte a la fórmula base ---
  if (INC_TREND_ON == "cohort") {
    form_inc <- update(base_formula, . ~ . + f(cohort_id, model=INC_COH_MODEL, constr=TRUE, scale.model=TRUE, hyper=hyper_coh))
  } else {
    if (use_weighted_cohort) {
      cov_coh <- incH %>% dplyr::count(cohort, name="n") %>%
        dplyr::right_join(tibble::tibble(cohort = lev_coh), by="cohort") %>%
        dplyr::mutate(w = dplyr::coalesce(n, 1L))
      Qw_inc <- make_Q_rw1_weighted(length(lev_coh), cov_coh$w)
      form_inc <- update(base_formula, . ~ . + f(cohort_id, model="generic0", Cmatrix=Qw_inc, constr=TRUE, hyper=hyper_coh))
    } else {
      form_inc <- update(base_formula, . ~ . + f(cohort_id, model=INC_COH_MODEL, constr=TRUE, scale.model=TRUE, hyper=hyper_coh))
    }
  }
  
  # --- Añadir el offset final ---
  form_inc <- update(form_inc, . ~ . + offset(logE + inc_tech_offset + coef_fc_offset_I))
  
  if (isTRUE(INC_PER_EXTRA_IID)) {
    form_inc <- update(form_inc, . ~ . +
                         f(period_iid, model = "iid",
                           hyper = list(prec = list(prior = "pc.prec",
                                                    param = c(INC_PER_IID_PC_U, INC_PER_IID_PC_A)))))
  }
  
  .check_no_mode_in_f(form_inc)
  fit_inc <- inla_tag("INC-HIST", formula = form_inc, family = "poisson", data = inc_inla,
                      weights = inc_inla$edge_weight,
                      control.fixed = ctrl_fix,
                      control.predictor = list(compute = TRUE),
                      control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE))
  
  # ---- Backtest temporal (Poisson) ----
  bt <- NULL
  if (isTRUE(BT_ENABLE) && isTRUE(BT_HOLDOUT_YEARS > 0)) {
    train_end <- period_max_i - as.integer(BT_HOLDOUT_YEARS)
    tag_bt <- sprintf("BT-INC-%s", as.character(unique(inc_inla$sex)[1]))
    bt <- backtest_inla_poisson(
      formula = form_inc,
      data    = inc_inla,
      train_end = train_end,
      tag = tag_bt,
      control.fixed = ctrl_fix
    )
  }
  
  # Futuro: escenarios + (opcional) pronóstico de APC→offset
  fut_grid <- pop_all %>%
    dplyr::filter(age >= age_min_i, age <= age_max_i, period > period_max_i) %>%
    ensure_exposure() %>%
    dplyr::mutate(cohort = period - age,
                  age_id = match(age, lev_age),
                  period_id = match(pmin(period, max(lev_per)), lev_per),
                  cohort_id = match(pmin(pmax(cohort, min(lev_coh)), max(lev_coh)), lev_coh),
                  E = pmax(as.numeric(exposure), 1e-12),
                  logE = log(E))
  fut_grid <- fut_grid %>%
    dplyr::mutate(period_iid = if (isTRUE(INC_PER_EXTRA_IID)) period_id else NA_integer_)

  border_diag_future <- make_incidence_border_diag(fut_grid, lev_per = lev_per, lev_coh = lev_coh,
                                                   hist_df = incH, last_hist_year = last_year_inc)
  if (nrow(border_diag_future)) {
    border_diag_join <- border_diag_future %>%
      dplyr::select(dplyr::any_of(c(
        "sex", "age", "period", "cohort",
        "period_raw", "mapped_period", "period_is_clamped", "period_shift",
        "cohort_raw", "mapped_cohort", "cohort_is_edge", "cohort_clamp_low", "cohort_clamp_high", "cohort_shift",
        "support_n", "support_frac", "horizon", "horizon_block"
      )))
    fut_grid <- fut_grid %>%
      dplyr::left_join(border_diag_join, by = c("sex", "age", "period", "cohort"))
  }
  
  # --- Aplica escenario de tendencia SÓLO si inc_degree > 0 ---
  if (inc_degree > 0) {
    fut_trend <- apply_trend_scenario_future(fut_grid$period, last_year_inc, mu_perI,
                                             degree = inc_degree, scenario = inc_scenario,
                                             delta = delta_inc, prefix = "inc_")
  } else {
    # Si no hay tendencia, crea las columnas con ceros para mantener la estructura del data.frame
    fut_trend <- tibble::tibble(inc_trend_t = 0, inc_trend_t2 = 0, inc_tech_offset = 0)
  }
  
  fut_grid <- dplyr::bind_cols(fut_grid, fut_trend) %>%
    dplyr::mutate(
      coef_fc_signal_I = 0,
      coef_fc_recenter_I = 0,
      coef_fc_offset_I = 0,
      y = NA_integer_
    )
  
  if (!identical(INC_COEF_FC_TARGET, "none")) {
    if (INC_COEF_FC_TARGET == "cohort") {
      fc_parts <- build_coef_fc_components(
        fit_inc$summary.random$cohort_id,
        lev_coh,
        fut_levels = fut_grid$cohort,
        ref_levels = pmin(pmax(fut_grid$cohort, min(lev_coh)), max(lev_coh)),
        method = INC_COEF_FC_METHOD
      )
    } else {
      fc_parts <- build_coef_fc_components(
        fit_inc$summary.random$period_id,
        lev_per,
        fut_levels = pmin(fut_grid$period, max(lev_per)),
        ref_levels = pmin(fut_grid$period, max(lev_per)),
        method = INC_COEF_FC_METHOD
      )
    }
    rec_locked <- apply_coef_fc_recenter_lock(fc_parts$recenter)
    off <- as.numeric(fc_parts$signal + rec_locked)
    off <- apply_coef_fc_lock(off)
    fut_grid$coef_fc_signal_I <- as.numeric(fc_parts$signal)
    fut_grid$coef_fc_recenter_I <- as.numeric(rec_locked)
    fut_grid$coef_fc_offset_I <- as.numeric(off)
  }
  
  data_all <- dplyr::bind_rows(
    inc_inla %>% dplyr::select(sex, age, period, cohort, E, logE, age_id, period_id, period_iid, cohort_id,
                               dplyr::any_of(c("period_raw", "mapped_period", "period_is_clamped", "period_shift",
                                               "cohort_raw", "mapped_cohort", "cohort_is_edge", "cohort_clamp_low", "cohort_clamp_high", "cohort_shift",
                                               "support_n", "support_frac", "horizon", "horizon_block")),
                               dplyr::any_of(c("edge_stage", "edge_geometry", "d_age_edge", "d_period_edge", "e_age", "e_period", "edge_score", "edge_weight")),
                               inc_trend_t, inc_trend_t2, inc_tech_offset, coef_fc_signal_I, coef_fc_recenter_I, coef_fc_offset_I, y),    
    fut_grid %>% dplyr::select(sex, age, period, cohort, E, logE, age_id, period_id, period_iid, cohort_id,
                               dplyr::any_of(c("period_raw", "mapped_period", "period_is_clamped", "period_shift",
                                               "cohort_raw", "mapped_cohort", "cohort_is_edge", "cohort_clamp_low", "cohort_clamp_high", "cohort_shift",
                                               "support_n", "support_frac", "horizon", "horizon_block")),
                               inc_trend_t, inc_trend_t2, inc_tech_offset, coef_fc_signal_I, coef_fc_recenter_I, coef_fc_offset_I, y)
  )
  
  .check_no_mode_in_f(form_inc)
  fit_pred <- inla_tag("INC-PRED", formula = form_inc, family = "poisson",
                       data = data_all, weights = data_all$edge_weight, control.fixed = ctrl_fix,
                       control.predictor = list(compute = TRUE),
                       control.compute = list(dic = FALSE, waic = FALSE, cpo = FALSE, config = TRUE))
  .chk("post INC-PRED #1")
  
  fv <- fit_pred$summary.fitted.values
  lp <- fit_pred$summary.linear.predictor
  grid_all <- data_all %>%
    dplyr::mutate(lp_mean = as.numeric(lp$mean),
                  lp_lwr  = as.numeric(lp$`0.025quant`),
                  lp_upr  = as.numeric(lp$`0.975quant`),
                  fv_mean = pmax(as.numeric(fv$mean), 0),
                  fv_lwr  = pmax(0, as.numeric(fv$`0.025quant`)),
                  fv_upr  = pmax(0, as.numeric(fv$`0.975quant`)),
                  mu_hat  = fv_mean,
                  mu_lwr  = fv_lwr,
                  mu_upr  = fv_upr,
                  .E_safe = pmax(exp(logE), 1e-12),
                  rate_total_fv = mu_hat / .E_safe,
                  rate_total_lp = exp(lp_mean) / .E_safe,
                  rate_hat = rate_total_lp,
                  rate_lwr = exp(lp_lwr) / .E_safe,
                  rate_upr = exp(lp_upr) / .E_safe) %>%
    dplyr::select(-.E_safe)

  grid_all <- apply_incidence_level_anchor(grid_all, last_hist_year = last_year_inc)
  grid_all <- apply_inc_coef_fc_posthoc_lock(grid_all, last_hist_year = last_year_inc)
  
  list(
    fit_inc   = fit_inc, 
    rates_all = grid_all %>% dplyr::select(sex, age, period, rate_hat, rate_lwr, rate_upr),
    rates_all_full = grid_all %>% dplyr::select(dplyr::any_of(c("sex", "age", "period", "cohort", "E", "logE", "p_cur", "delta_p_cur", "quit_flow", "p_never", "p_former_total", "noncurrent_rescale", "q_eff", "offset_prev_rr", "prev_source", "within_prev_age_support", "within_prev_period_support", "within_prev_observed_support", "within_prev_support", "prev_scenario_name", "prev_scenario_applied", "prev_inc_channel_requested", "prev_inc_channel_used", "used_A_I_in_main_channel", "prev_inc_channel_requested", "prev_inc_channel_used", "used_A_I_in_main_channel",
                                                                "period_raw", "mapped_period", "period_is_clamped", "period_shift",
                                                                "cohort_raw", "mapped_cohort", "cohort_is_edge", "cohort_clamp_low", "cohort_clamp_high", "cohort_shift",
                                                                "support_n", "support_frac", "horizon", "horizon_block",
                                                                "inc_tech_offset", "coef_fc_signal_I", "coef_fc_recenter_I", "coef_fc_offset_I", "coef_fc_offset_I_raw", "coef_fc_offset_I_effective", "coef_fc_posthoc_adj", "coef_fc_posthoc_lock_mode", "lp_mean", "lp_lwr", "lp_upr", "fv_mean", "fv_lwr", "fv_upr", "mu_hat", "mu_lwr", "mu_upr", "rate_total_fv", "rate_total_lp", "eta_apc_manual", "eta_offset_manual", "eta_total_manual", "rate_manual", "mu_manual", "lp_gap_manual", "rate_hat", "rate_lwr", "rate_upr"))),
    border_diag_future = border_diag_future,
    border_diag_summary = summarise_incidence_border_diag(border_diag_future, exposure_col = "E"),
    lev_inc   = list(age = lev_age, period = sort(unique(grid_all$period)), cohort = lev_coh),
    last_year_inc = last_year_inc, bt = bt,
    trend_meta = list(center = mu_perI, degree = inc_degree)
  )
}




# --------- APC de INCIDENCIA | PREVALENCIA (corregida) ----------

# =========================================================
# Hard-coded excess-risk reversal schedules after smoking cessation
# Target quantity:
#   w(s) = (RR(s) - 1) / (RR(0) - 1)
# Curated for model use:
# - non-negative
# - non-increasing in years since quit
# - linear interpolation between anchors, clamped at endpoints
# =========================================================


fit_apc_incidence_cond_prev <- function(inc_hist, pop_all, fit_prev,
                                        age_min_i = 35, age_max_i = 89,
                                        age_min_p = AGE_P_MIN, age_max_p = AGE_P_MAX,
                                        period_min_i = 1998, period_max_i = 2022,
                                        prev_sign = 1,
                                        A_I = 30, w_I = 1,
                                        use_weighted_cohort = TRUE,
                                        sd_theta_IP = SD_THETA_IP,
                                        gammaP_method = GAMMAP_METHOD, trend_type = TREND_TYPE,
                                        inc_trend_scenario = INC_TREND_SCENARIO,
                                        delta_inc = DELTA_INC,
                                        sd_beta_I = SD_BETA_I,
                                        inc_degree = INC_TREND_DEGREE,
                                        prev_scenario = PREV_SCENARIO,
                                        prev_scenario_axis = PREV_SCENARIO_AXIS,
                                        prev_annual_rate = PREV_ANNUAL_RATE,
                                        prev_annual_rate_up = PREV_ANNUAL_RATE_UP,
                                        prev_annual_rate_down = PREV_ANNUAL_RATE_DOWN,
                                        prev_annual_rate_down3 = PREV_ANNUAL_RATE_DOWN3,
                                        prev_base_year = PREV_BASE_YEAR,
                                        quit_mode = QUIT_MODE,
                                        quit_floor_sd = QUIT_FLOOR_SD,
                                        quit_floor_sd_M = QUIT_FLOOR_SD_M,
                                        quit_floor_sd_F = QUIT_FLOOR_SD_F,
                                        quit_half_life = QUIT_HALF_LIFE,
                                        quit_ramp_years = QUIT_RAMP_YEARS,
                                        prev_base_M = PREV_BASE_M,
                                        prev_base_F = PREV_BASE_F,
                                        prev_base_default = PREV_BASE_DEFAULT,
                                        rr_inc = NA_real_,
                                        prev_base_prob = NA_real_,
                                        cause_id = NA_character_,
                                        BETA_P_POSTFIT_RULE = "floor0",
                                        prev_inc_channel_mode = PREV_INC_CHANNEL_MODE) {
  
  if (identical(INC_TREND_ON, "none")) inc_degree <- 0L
  beta_mode <- "fixed_rr_offset"
  rr_use <- suppressWarnings(as.numeric(rr_inc))[1]
  if (!is.finite(rr_use) || rr_use <= 1) rr_use <- 2.0
  pop_all <- ensure_exposure(pop_all)
  prev_cfg <- make_prev_config(
    scenario = prev_scenario,
    axis = prev_scenario_axis,
    annual_rate = prev_annual_rate,
    annual_rate_up = prev_annual_rate_up,
    annual_rate_down = prev_annual_rate_down,
    annual_rate_down3 = prev_annual_rate_down3,
    base_year = prev_base_year,
    backbone = PREV_BACKBONE,
    quit_mode = quit_mode,
    quit_floor_sd = quit_floor_sd,
    quit_floor_sd_M = quit_floor_sd_M,
    quit_floor_sd_F = quit_floor_sd_F,
    quit_half_life = quit_half_life,
    quit_ramp_years = quit_ramp_years,
    prev_base_M = prev_base_M,
    prev_base_F = prev_base_F,
    prev_base_default = prev_base_default
  )
  prev_inc_channel_requested <- "stock_former"
  prev_inc_channel_used <- "stock_former"
  used_A_I_in_main_channel <- FALSE
  scenario_embedded_in_prev <- TRUE
  
  incH <- inc_hist %>%
    dplyr::filter(age >= age_min_i, age <= age_max_i,
                  period >= period_min_i, period <= period_max_i) %>%
    dplyr::arrange(sex, period, age)         %>%
    dplyr::left_join(pop_all, by = c("period","age","sex")) %>%
    ensure_exposure()                        %>%
    dplyr::mutate(cohort = period - age,
                  y = as.integer(round(cases)),
                  E = pmax(as.numeric(exposure), 1e-12),
                  logE = log(E))
  
  # γ^P (histórico + forecast)
  scP <- fit_prev$summary.random$cohort_id
  lev_coh_prev <- sort(unique(as.integer(fit_prev$.args$data$cohort)))
  gammaP_hist <- tibble::tibble(cohort = lev_coh_prev[scP$ID],
                                gammaP = scP$mean - mean(scP$mean)) %>%
    dplyr::arrange(cohort)
  coh_all_needed <- sort(unique((pop_all$period - pop_all$age)))
  gammaP_fut <- forecast_gammaP(gammaP_hist,
                                setdiff(coh_all_needed, gammaP_hist$cohort),
                                method = gammaP_method, trend_type = trend_type)
  gammaP_fut <- adjust_gammaP_future(gammaP_hist, gammaP_fut,
                                     scenario = "freeze",
                                     annual_rate = prev_annual_rate,
                                     annual_rate_up = prev_annual_rate_up,
                                     annual_rate_down = prev_annual_rate_down,
                                     annual_rate_down3 = prev_annual_rate_down3,
                                     base_year = prev_base_year)
  gammaP_all <- dplyr::bind_rows(gammaP_hist, gammaP_fut)
  
  lev_age <- sort(unique(incH$age))
  lev_per <- sort(unique(incH$period))
  lev_coh <- sort(unique(incH$cohort))
  mu_perI <- mean(incH$period)
  last_hist_I <- max(incH$period)
  Xconstr_per_IP <- make_slope_constr(lev_per, mu_perI)
  
  # ---------- HISTÓRICO y FUTURO (canal PREV -> INC) ----------
  inc_base_raw <- incH %>%
    dplyr::mutate(age_id = match(age, lev_age),
                  period_id = match(period, lev_per),
                  cohort_id = match(cohort, lev_coh))
  fut_base_raw <- pop_all %>%
    dplyr::filter(age >= age_min_i, age <= age_max_i, period > period_max_i) %>%
    ensure_exposure() %>%
    dplyr::mutate(cohort = period - age,
                  age_id = match(age, lev_age),
                  period_id = match(pmin(period, max(lev_per)), lev_per),
                  cohort_id = match(pmin(pmax(cohort, min(lev_coh)), max(lev_coh)), lev_coh),
                  E = pmax(as.numeric(exposure), 1e-12),
                  logE = log(E))

  combined_prev <- dplyr::bind_rows(
    inc_base_raw %>% dplyr::mutate(.is_future_previnc = FALSE),
    fut_base_raw %>% dplyr::mutate(.is_future_previnc = TRUE)
  ) %>%
    build_prev_rr_offset_stock_for_inc(
      fit_prev = fit_prev,
      cause_id = cause_id,
      rr_inc = rr_use,
      prev_inla = tryCatch(fit_prev$.args$data, error = function(e) NULL),
      sex_sel = tryCatch(as.character(unique(incH$sex)[1]), error = function(e) NA_character_),
      gammaP_method = gammaP_method,
      trend_type = trend_type,
      prev_cfg = prev_cfg,
      age_min_p = age_min_p,
      age_max_p = age_max_p,
      backcast_period_mode = PREV_BACKCAST_MODE,
      backcast_cohort_mode = PREV_BACKCAST_COHORT_MODE,
      post65_mode = PREV_POST65_MODE,
      quit_horizon_years = PREV_INC_MAX_QUIT_YEARS,
    )
  
  inc_base_z <- combined_prev %>% dplyr::filter(!.data$.is_future_previnc) %>% dplyr::select(-.data$.is_future_previnc)
  fut_base_z <- combined_prev %>% dplyr::filter(.data$.is_future_previnc) %>% dplyr::select(-.data$.is_future_previnc)
  s_histP_val <- NA_real_
  base_year_val <- prev_base_year
  
  # Centrado del offset para estabilizar INLA (centrado en el último año para evitar saltos)
  offset_raw  <- dplyr::coalesce(as.numeric(inc_base_z$offset_prev_rr), 0)
  offset_mean <- mean(inc_base_z$offset_prev_rr[inc_base_z$period == max(inc_base_z$period, na.rm=TRUE)], na.rm = TRUE)
  offset_epi_hist <- offset_raw - offset_mean
  
  # Limpiar columnas técnicas previas usando R base para evitar errores de duplicados
  tech_cols <- c("inc_tech_offset", "coef_fc_signal_I", "coef_fc_recenter_I", 
                 "coef_fc_offset_I_epi", "coef_fc_offset_I_apc", "coef_fc_offset_I")
  
  
  inc_inla <- inc_base_z[, !(colnames(inc_base_z) %in% tech_cols), drop = FALSE]
  
  # Asegurar que inc_inla es un data frame limpio
  inc_inla <- as.data.frame(inc_inla)
  
  inc_inla <- inc_inla %>%
    dplyr::bind_cols(make_trend_vars(.$period, mu_perI, inc_degree, prefix = "inc_")) %>%
    dplyr::mutate(
      y = as.integer(round(dplyr::coalesce(cases, 0))),
      inc_tech_offset = 0,
      coef_fc_signal_I = 0,
      coef_fc_recenter_I = offset_mean,
      coef_fc_offset_I_epi = offset_epi_hist,
      coef_fc_offset_I_apc = 0,
      coef_fc_offset_I = coef_fc_offset_I_epi
    )
  if (!identical(beta_mode, "fixed_rr_offset")) inc_inla$z_prev <- prev_sign * inc_inla$z_prev
  attr(inc_inla, "s_histP") <- s_histP_val
  attr(inc_inla, "base_year") <- base_year_val
  inc_inla <- inc_inla %>% dplyr::mutate(period_iid = period_id)

  fut_grid <- fut_base_z %>%
    dplyr::mutate(
      # No añadimos columnas técnicas aquí todavía para evitar duplicados con fut_trend
      period_iid = period_id
    )
  
  if (!identical(beta_mode, "fixed_rr_offset")) {
     # En modos antiguos z_prev es la señal cruda
     fut_grid$z_prev <- prev_sign * fut_grid$z_prev
  }
  
  attr(fut_grid, "s_histP") <- s_histP_val
  attr(fut_grid, "base_year") <- base_year_val

  border_diag_future <- make_incidence_border_diag(fut_grid, lev_per = lev_per, lev_coh = lev_coh,
                                                   hist_df = incH, last_hist_year = last_hist_I)
  if (nrow(border_diag_future)) {
    border_diag_join <- border_diag_future %>%
      dplyr::select(dplyr::any_of(c(
        "sex", "age", "period", "cohort",
        "period_raw", "mapped_period", "period_is_clamped", "period_shift",
        "cohort_raw", "mapped_cohort", "cohort_is_edge", "cohort_clamp_low", "cohort_clamp_high", "cohort_shift",
        "support_n", "support_frac", "horizon", "horizon_block"
      )))
    fut_grid <- fut_grid %>%
      dplyr::left_join(border_diag_join, by = c("sex", "age", "period", "cohort"))
  }
  
  # ---------- PRIMER AJUSTE ----------
  hyper_age <- pc_hyper(INC_AGE_PC_U, INC_AGE_PC_A)
  hyper_per <- pc_hyper(INC_PER_PC_U, INC_PER_PC_A)
  hyper_coh <- pc_hyper(INC_COH_PC_U, INC_COH_PC_A)
  
  # --- Fórmula condicional según la perilla INC_TREND_ON ---
  if (INC_TREND_ON == "period") {
    .bapc_verbose("[fit_apc_incidence_cond_prev] Linear trend on PERIOD.")
    if (identical(beta_mode, "fixed_rr_offset")) {
      base_formula <- y ~ 1 + inc_trend_t +
        f(age_id, model=INC_AGE_MODEL, constr=TRUE, scale.model=TRUE, hyper=hyper_age) +
        f(period_id, model = INC_PER_MODEL, constr = TRUE, extraconstr = Xconstr_per_IP, scale.model = TRUE, hyper = hyper_per)
      ctrl_fix <- list(
        mean = list(inc_trend_t = 0),
        prec = list(inc_trend_t = 1 / sd_beta_I^2)
      )
    } else {
      base_formula <- y ~ 1 + z_prev + inc_trend_t +
        f(age_id, model=INC_AGE_MODEL, constr=TRUE, scale.model=TRUE, hyper=hyper_age) +
        f(period_id, model = INC_PER_MODEL, constr = TRUE, extraconstr = Xconstr_per_IP, scale.model = TRUE, hyper = hyper_per)
      ctrl_fix <- list(
        mean = list(z_prev = 0, inc_trend_t = 0),
        prec = list(z_prev = 1 / sd_theta_IP^2,
                    inc_trend_t = 1 / sd_beta_I^2)
      )
    }

  } else if (INC_TREND_ON == "cohort") {
    .bapc_verbose("[fit_apc_incidence_cond_prev] Trend on COHORT (rw2).")
    if (identical(beta_mode, "fixed_rr_offset")) {
      base_formula <- y ~ 1 +
        f(age_id, model=INC_AGE_MODEL, constr=TRUE, scale.model=TRUE, hyper=hyper_age) +
        f(period_id, model = INC_PER_MODEL, constr = TRUE, extraconstr = Xconstr_per_IP, scale.model = TRUE, hyper = hyper_per)
      ctrl_fix <- list()
    } else {
      base_formula <- y ~ 1 + z_prev +
        f(age_id, model=INC_AGE_MODEL, constr=TRUE, scale.model=TRUE, hyper=hyper_age) +
        f(period_id, model = INC_PER_MODEL, constr = TRUE, extraconstr = Xconstr_per_IP, scale.model = TRUE, hyper = hyper_per)
      ctrl_fix <- list(
        mean = list(z_prev = 0),
        prec = list(z_prev = 1 / sd_theta_IP^2)
      )
    }
    use_weighted_cohort <- FALSE  # rw2 cohorte

  } else if (INC_TREND_ON == "none") {
    message("[fit_apc_incidence_cond_prev] Sin tendencia explícita (APC libre).")
    if (identical(beta_mode, "fixed_rr_offset")) {
      base_formula <- y ~ 1 +
        f(age_id, model=INC_AGE_MODEL, constr=TRUE, scale.model=TRUE, hyper=hyper_age) +
        f(period_id, model = INC_PER_MODEL, constr = TRUE, scale.model = TRUE, hyper = hyper_per) +
        f(cohort_id, model = INC_COH_MODEL, constr = TRUE, scale.model = TRUE, hyper = hyper_coh)
      ctrl_fix <- list()
    } else {
      base_formula <- y ~ 1 + z_prev +
        f(age_id, model=INC_AGE_MODEL, constr=TRUE, scale.model=TRUE, hyper=hyper_age) +
        f(period_id, model = INC_PER_MODEL, constr = TRUE, scale.model = TRUE, hyper = hyper_per) +
        f(cohort_id, model = INC_COH_MODEL, constr = TRUE, scale.model = TRUE, hyper = hyper_coh)
      ctrl_fix <- list(
        mean = list(z_prev = 0),
        prec = list(z_prev = 1 / sd_theta_IP^2)
      )
    }

  } else {
    stop("INC_TREND_ON must be 'period', 'cohort' or 'none'.")
  }
  
  # --- Añadir el efecto de cohorte a la fórmula base ---
  if (INC_TREND_ON == "cohort") {
    form_inc <- update(base_formula, . ~ . + f(cohort_id, model=INC_COH_MODEL, constr=TRUE, scale.model=TRUE, hyper=hyper_coh))
  } else {
    if (use_weighted_cohort) {
      cov_coh <- incH %>% dplyr::count(cohort, name="n") %>%
        dplyr::right_join(tibble::tibble(cohort = lev_coh), by="cohort") %>%
        dplyr::mutate(w = dplyr::coalesce(n, 1L))
      Qw_I <- make_Q_rw1_weighted(length(lev_coh), cov_coh$w)
      form_inc <- update(base_formula, . ~ . + f(cohort_id, model="generic0", Cmatrix=Qw_I, constr=TRUE, hyper=hyper_coh))
    } else {
      form_inc <- update(base_formula, . ~ . + f(cohort_id, model=INC_COH_MODEL, constr=TRUE, scale.model=TRUE, hyper=hyper_coh))
    }
  }
  
  # --- Añadir el offset final (z_prev es ahora el log-riesgo epidemiológico completo) ---
  if (identical(beta_mode, "fixed_rr_offset")) {
    form_inc <- update(form_inc, . ~ . + offset(logE + inc_tech_offset + coef_fc_offset_I))
  } else {
    form_inc <- update(form_inc, . ~ . + offset(logE + inc_tech_offset + coef_fc_offset_I + z_prev))
  }
  
  if (isTRUE(INC_PER_EXTRA_IID)) {
    form_inc <- update(form_inc, . ~ . +
                         f(period_iid, model = "iid",
                           hyper = list(prec = list(prior = "pc.prec",
                                                    param = c(INC_PER_IID_PC_U, INC_PER_IID_PC_A)))))
  }
  
  .chk("pre INC-2")
  .check_no_mode_in_f(form_inc)
  fit_inc <- inla_tag("INC-2",formula = form_inc, family = "poisson", data = inc_inla,
                        control.fixed = ctrl_fix,
                        control.predictor = list(compute = TRUE),
                        control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE))
  
  # ---- Backtest temporal (Poisson) ----
  bt <- NULL
  if (isTRUE(BT_ENABLE) && isTRUE(BT_HOLDOUT_YEARS > 0)) {
    train_end <- period_max_i - as.integer(BT_HOLDOUT_YEARS)
    tag_bt <- sprintf("BT-INCIP-%s", as.character(unique(inc_inla$sex)[1]))
    bt <- backtest_inla_poisson(
      formula = form_inc,
      data    = inc_inla,
      train_end = train_end,
      tag = tag_bt,
      control.fixed = ctrl_fix
    )
  }
  
  # ---------- Auto-arreglo de dirección del efecto de P (1 refit si hace falta) ----------
  betaP <- if (identical(beta_mode, "fixed_rr_offset")) 1 else tryCatch(fit_inc$summary.fixed["z_prev","mean"], error = function(e) NA_real_)
  if (is.finite(betaP) && betaP < 0) {
  }
  
  # === APLICAR SHOCKS DE ESCENARIO EN z_prev YA ORIENTADO ===
  if (!isTRUE(scenario_embedded_in_prev)) {
    # Recuperar s_histP (normaliza magnitud del paso) y base_year
  }
  # ---------- Completar futuro: términos de tendencia y offsets ----------
  # --- Aplica escenario de tendencia SÓLO si inc_degree > 0 ---
  if (inc_degree > 0) {
    fut_trend <- apply_trend_scenario_future(fut_grid$period, last_hist_I, mu_perI,
                                             degree = inc_degree, scenario = inc_trend_scenario,
                                             delta = delta_inc, prefix = "inc_")
  } else {
    # Si no hay tendencia, crea las columnas con ceros para mantener la estructura del data.frame
    fut_trend <- tibble::tibble(inc_trend_t = 0, inc_trend_t2 = 0, inc_tech_offset = 0)
  }
  
  offset_epi_fut <- .safe_num(fut_grid$offset_prev_rr) - offset_mean
  fut_grid <- dplyr::bind_cols(fut_grid, fut_trend) %>%
    dplyr::mutate(
      coef_fc_signal_I = 0,
      coef_fc_recenter_I = offset_mean,
      coef_fc_offset_I_epi = offset_epi_fut,
      coef_fc_offset_I_apc = 0,
      coef_fc_offset_I = coef_fc_offset_I_epi,
      y = NA_integer_,
      edge_stage = "inc",
      edge_geometry = EDGE_WEIGHT_GEOMETRY,
      d_age_edge = NA_integer_,
      d_period_edge = NA_integer_,
      e_age = 0,
      e_period = 0,
      edge_score = 0,
      edge_weight = 1
    )
  
  if (!identical(INC_COEF_FC_TARGET, "none")) {
    if (INC_COEF_FC_TARGET == "cohort") {
      fc_parts <- build_coef_fc_components(
        fit_inc$summary.random$cohort_id,
        lev_coh,
        fut_levels = fut_grid$cohort,
        ref_levels = pmin(pmax(fut_grid$cohort, min(lev_coh)), max(lev_coh)),
        method = gammaP_method,
        trend_type = trend_type
      )
    } else {
      fc_parts <- build_coef_fc_components(
        fit_inc$summary.random$period_id,
        lev_per,
        fut_levels = pmin(fut_grid$period, max(lev_per)),
        ref_levels = pmin(fut_grid$period, max(lev_per)),
        method = gammaP_method,
        trend_type = trend_type
      )
    }
    rec_locked <- apply_coef_fc_recenter_lock(fc_parts$recenter)
    off <- as.numeric(fc_parts$signal + rec_locked)
    off <- apply_coef_fc_lock(off)
    fut_grid$coef_fc_signal_I <- as.numeric(fc_parts$signal)
    fut_grid$coef_fc_recenter_I <- as.numeric(rec_locked)
    fut_grid$coef_fc_offset_I_apc <- as.numeric(off)
    fut_grid$coef_fc_offset_I <- fut_grid$coef_fc_offset_I_epi + fut_grid$coef_fc_offset_I_apc
  }
  
  # ---------- Predicción en historia + futuro ----------
  req_cols <- c("sex", "age", "period", "cohort", "E", "logE", "age_id", "period_id", "period_iid", "cohort_id",
                "inc_trend_t", "inc_trend_t2", "inc_tech_offset", "z_prev",
                "coef_fc_signal_I", "coef_fc_recenter_I", "coef_fc_offset_I_epi", "coef_fc_offset_I_apc", "coef_fc_offset_I", "y")
  
  data_all <- dplyr::bind_rows(
    inc_inla %>% dplyr::select(dplyr::any_of(c(req_cols, "period_raw", "mapped_period", "period_is_clamped", "period_shift",
                                               "cohort_raw", "mapped_cohort", "cohort_is_edge", "cohort_clamp_low", "cohort_clamp_high", "cohort_shift",
                                               "support_n", "support_frac", "horizon", "horizon_block",
                                               "z_prev", "p_cur", "delta_p_cur", "quit_flow", "p_never", "p_former_total",
                                               "noncurrent_rescale", "q_eff", "offset_prev_rr", "prev_source", "within_prev_support"))),
    fut_grid %>% dplyr::select(dplyr::any_of(c(req_cols, "period_raw", "mapped_period", "period_is_clamped", "period_shift",
                                               "cohort_raw", "mapped_cohort", "cohort_is_edge", "cohort_clamp_low", "cohort_clamp_high", "cohort_shift",
                                               "support_n", "support_frac", "horizon", "horizon_block",
                                               "z_prev", "p_cur", "delta_p_cur", "quit_flow", "p_never", "p_former_total",
                                               "noncurrent_rescale", "q_eff", "offset_prev_rr", "prev_source", "within_prev_support")))
  )
  
  # Validar integridad de data_all
  check_cols <- c("logE", "inc_trend_t", "inc_trend_t2", "inc_tech_offset", "z_prev", "coef_fc_offset_I")
  for (col in check_cols) {
    if (col %in% names(data_all)) {
      vec <- data_all[[col]]
      if (any(is.na(vec))) message("!!! ADVERTENCIA: ", col, " tiene NAs en data_all")
    }
  }

  # Offset continuity diagnostic
  if (isTRUE(BAPC_VERBOSE)) {
    hist_2022 <- data_all %>% dplyr::filter(period == 2022) %>% dplyr::pull(coef_fc_offset_I)
    fut_2023  <- data_all %>% dplyr::filter(period == 2023) %>% dplyr::pull(coef_fc_offset_I)
    message(sprintf("OFFSET: Mean=%f | 2022(avg)=%f | 2023(avg)=%f", 
                    offset_mean, mean(hist_2022), mean(fut_2023)))
  }

  .check_no_mode_in_f(form_inc)
  fit_inc_pred <- inla_tag("INC-PRED", formula = form_inc, family = "poisson",
                             data = data_all, control.fixed = ctrl_fix,
                             control.predictor = list(compute = TRUE),
                             control.compute = list(dic = FALSE, waic = FALSE, cpo = FALSE, config = TRUE))

  fv_all <- fit_inc_pred$summary.fitted.values
  lp_all <- fit_inc_pred$summary.linear.predictor

  .map_re_mean <- function(sr, ids) {
    if (is.null(sr) || length(ids) == 0L) return(rep(0, length(ids)))
    key <- as.character(sr[[1]])
    val <- as.numeric(sr[["mean"]])
    out <- val[match(as.character(ids), key)]
    out[is.na(out)] <- 0
    out
  }
  .fix_mean <- function(sf, nm, default = 0) {
    if (is.null(sf) || is.null(rownames(sf)) || !(nm %in% rownames(sf))) return(default)
    as.numeric(sf[nm, "mean"])
  }

  sf_pred <- fit_inc_pred$summary.fixed
  age_re_mean <- .map_re_mean(fit_inc_pred$summary.random$age_id, data_all$age_id)
  per_re_mean <- .map_re_mean(fit_inc_pred$summary.random$period_id, data_all$period_id)
  coh_re_mean <- .map_re_mean(fit_inc_pred$summary.random$cohort_id, data_all$cohort_id)
  iid_re_mean <- .map_re_mean(fit_inc_pred$summary.random$period_iid, data_all$period_iid)

  b0_mean  <- .fix_mean(sf_pred, "(Intercept)", 0)
  bt_mean  <- .fix_mean(sf_pred, "inc_trend_t", 0)
  bt2_mean <- .fix_mean(sf_pred, "inc_trend_t2", 0)
  bz_mean  <- .fix_mean(sf_pred, "z_prev", 0)

  grid_all <- data_all %>%
    dplyr::mutate(lp_mean = as.numeric(lp_all$mean),
                  lp_lwr  = as.numeric(lp_all$`0.025quant`),
                  lp_upr  = as.numeric(lp_all$`0.975quant`),
                  fv_mean = pmax(as.numeric(fv_all$mean), 0),
                  fv_lwr  = pmax(0, as.numeric(fv_all$`0.025quant`)),
                  fv_upr  = pmax(0, as.numeric(fv_all$`0.975quant`)),
                  mu_hat = fv_mean,
                  mu_lwr = fv_lwr,
                  mu_upr = fv_upr,
                  .E_safe = pmax(exp(logE), 1e-12),
                  .offset_mult = exp(dplyr::coalesce(inc_tech_offset, 0) + dplyr::coalesce(coef_fc_offset_I, 0)),
                  rate_from_fv_over_E = mu_hat / .E_safe,
                  rate_from_fv_over_E_times_offset = rate_from_fv_over_E * .offset_mult,
                  rate_from_lp_over_E = exp(lp_mean) / .E_safe,
                  rate_from_lp_over_E_times_offset = rate_from_lp_over_E * .offset_mult,
                  rate_blend_geom = sqrt(pmax(rate_from_fv_over_E_times_offset, 0) * pmax(rate_from_lp_over_E, 0)),
                  rate_blend_geom_times_offset = rate_blend_geom * .offset_mult,
                  rate_blend_arith = 0.5 * (rate_from_fv_over_E_times_offset + rate_from_lp_over_E),
                  rate_blend_logmid = exp(0.5 * (log(pmax(rate_from_fv_over_E_times_offset, 1e-300)) + log(pmax(rate_from_lp_over_E, 1e-300)))),
                  eta_apc_manual = b0_mean + age_re_mean + per_re_mean + coh_re_mean + iid_re_mean +
                                   bt_mean * dplyr::coalesce(inc_trend_t, 0) +
                                   bt2_mean * dplyr::coalesce(inc_trend_t2, 0) +
                                   bz_mean * dplyr::coalesce(z_prev, 0),
                  eta_offset_manual = dplyr::coalesce(inc_tech_offset, 0) + 
                                      dplyr::coalesce(coef_fc_offset_I, 0),
                  eta_total_manual = eta_apc_manual + eta_offset_manual,
                  rate_manual = exp(eta_total_manual),
                  mu_manual = rate_manual * .E_safe,
                  lp_gap_manual = lp_mean - (log(.E_safe) + eta_total_manual),
                  .adj_center_manual = 1, 
                  rate_hat = rate_from_lp_over_E,
                  rate_lwr = exp(lp_lwr) / .E_safe,
                  rate_upr = exp(lp_upr) / .E_safe,
                  mu_hat = rate_hat * .E_safe,
                  mu_lwr = rate_lwr * .E_safe,
                  mu_upr = rate_upr * .E_safe) %>%
    dplyr::select(-.E_safe, -.offset_mult, -.adj_center_manual)

  if (isTRUE(BAPC_VERBOSE)) {
    message(sprintf("LP GAP: Mean Gap=%f | Max Gap=%f", 
                    mean(grid_all$lp_gap_manual, na.rm=TRUE), 
                    max(abs(grid_all$lp_gap_manual), na.rm=TRUE)))
  }

  grid_all <- grid_all %>%
    dplyr::mutate(
      prev_inc_channel_requested = prev_inc_channel_requested,
      prev_inc_channel_used = prev_inc_channel_used,
      used_A_I_in_main_channel = used_A_I_in_main_channel
    )

  grid_all <- apply_incidence_level_anchor(grid_all, last_hist_year = last_hist_I, y_col = "cases")
  grid_all <- apply_inc_coef_fc_posthoc_lock(grid_all, last_hist_year = last_hist_I)
  
  # Forzar que los resultados finales usen el predictor de INLA anclado
  E_safe <- pmax(grid_all$E, 1e-12)
  grid_all$rate_hat <- exp(grid_all$lp_mean) / E_safe
  grid_all$rate_lwr <- exp(grid_all$lp_lwr) / E_safe
  grid_all$rate_upr <- exp(grid_all$lp_upr) / E_safe
  grid_all$mu_hat   <- grid_all$rate_hat * E_safe
  grid_all$mu_lwr   <- grid_all$rate_lwr * E_safe
  grid_all$mu_upr   <- grid_all$rate_upr * E_safe

  beta_P_hat <- if (identical(beta_mode, "fixed_rr_offset")) 1 else tryCatch(as.numeric(fit_inc$summary.fixed["z_prev","mean"]), error = function(e) NA_real_)
  
  # -- Asegurar β_P no-negativo según la regla post-fit (sin re-ajustar) --
  beta_post <- .beta_postfit_transform(beta_P_hat, rule = BETA_P_POSTFIT_RULE)
  beta_P_eff <- beta_post$eff
  beta_P_zeroed <- beta_post$zeroed
  if (!identical(beta_mode, 'fixed_rr_offset') && isTRUE(is.finite(beta_P_eff)) && isTRUE(is.finite(beta_P_hat)) && !isTRUE(all.equal(beta_P_eff, beta_P_hat))) {
    grid_all <- grid_all %>%
      dplyr::mutate(
        .adj_beta = exp((beta_P_eff - beta_P_hat) * dplyr::coalesce(z_prev, 0)),
        rate_hat  = pmax(rate_hat * .adj_beta, 1e-12),
        rate_lwr  = ifelse(is.finite(rate_lwr), pmax(rate_lwr * .adj_beta, 1e-12), NA_real_),
        rate_upr  = ifelse(is.finite(rate_upr), pmax(rate_upr * .adj_beta, 1e-12), NA_real_)
      ) %>%
      dplyr::select(-.adj_beta)
  }

  if (!is.data.frame(grid_all)) {
    message("!!! ERROR: grid_all no es data.frame. Clase: ", paste(class(grid_all), collapse=", "))
    warning("fit_apc_incidence_cond_prev: grid_all is not a dataframe at the end. Returning NULL rates.")
    rates_all_res <- NULL
    rates_all_full_res <- NULL
  } else {
    .bapc_verbose("grid_all OK. Filas: ", nrow(grid_all), " Cols: ", ncol(grid_all))
    rates_all_res <- grid_all %>% dplyr::select(dplyr::any_of(c("sex", "age", "period", "rate_hat", "rate_lwr", "rate_upr")))
    rates_all_full_res <- grid_all
  }

  list(
    fit_inc   = fit_inc,
    rates_all = rates_all_res,
    rates_all_full = rates_all_full_res,
    border_diag_future = border_diag_future,
    border_diag_summary = summarise_incidence_border_diag(border_diag_future, exposure_col = "E"),
    lev_inc   = list(age = lev_age, period = sort(unique(grid_all$period)), cohort = lev_coh),
    last_year_inc = last_hist_I,
    trend_meta = list(center = mu_perI, degree = inc_degree),
    prev_inc_channel_requested = prev_inc_channel_requested,
    prev_inc_channel_used = prev_inc_channel_used,
    used_A_I_in_main_channel = used_A_I_in_main_channel,
    beta_P = beta_P_hat, bt = bt,
    beta_P_pos = beta_P_eff,
    beta_P_eff = beta_P_eff,
    beta_P_zeroed = beta_P_zeroed,
    beta_P_rule = if (identical(beta_mode, "fixed_rr_offset")) "fixed_rr_offset" else BETA_P_POSTFIT_RULE,
    s_histP = s_histP_val,
    z_hist = as.data.frame(inc_inla %>% dplyr::select(dplyr::any_of(c("sex", "age", "period", "cohort", "q_eff", "z_prev", "p_cur", "delta_p_cur", "quit_flow", "p_never", "p_former_total", "noncurrent_rescale", "offset_prev_rr", "prev_source", "within_prev_age_support", "within_prev_period_support", "within_prev_observed_support", "within_prev_support", "prev_scenario_name", "prev_scenario_applied")))),
    z_future = as.data.frame(fut_grid %>% dplyr::select(dplyr::any_of(c("sex", "age", "period", "cohort", "q_eff", "z_prev", "p_cur", "delta_p_cur", "quit_flow", "p_never", "p_former_total", "noncurrent_rescale", "offset_prev_rr", "prev_source", "within_prev_age_support", "within_prev_period_support", "within_prev_observed_support", "within_prev_support", "prev_scenario_name", "prev_scenario_applied", "coef_fc_signal_I", "coef_fc_recenter_I", "coef_fc_offset_I_epi", "coef_fc_offset_I_apc", "coef_fc_offset_I",
                                                           "period_raw", "mapped_period", "period_is_clamped", "period_shift",
                                                           "cohort_raw", "mapped_cohort", "cohort_is_edge", "cohort_clamp_low", "cohort_clamp_high", "cohort_shift",
                                                           "support_n", "support_frac", "horizon", "horizon_block")))),
    gammaP_all = gammaP_all
  )
}

