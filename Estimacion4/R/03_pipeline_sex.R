# 5) Sex-specific Pipeline (Prev → Inc → Mort)
# =============================================================

run_pipeline_sex <- function(
    sex_sel = c("M","F"),
    period_min_m = PERIOD_M_MIN, period_max_m = PERIOD_M_MAX,
    period_min_p = PERIOD_M_MIN, period_max_p = PERIOD_M_MAX,
    age_min_m = AGE_M_MIN, age_max_m = AGE_M_MAX,
    age_min_p = AGE_P_MIN, age_max_p = AGE_P_MAX,
    age_min_i = AGE_I_MIN, age_max_i = AGE_I_MAX,
    L_I = L_I_DEFAULT, Da_I = DA_I, bridge_inc_years = BRIDGE_INC_YEARS,
    tech_scenario = MORT_TREND_SCENARIO, mort_bapc_trend_scenario = MORT_BAPC_TREND_SCENARIO, delta_tech = DELTA_TECH,
    inc_include_trend = (INC_TREND_DEGREE > 0),
    inc_trend_scenario = INC_TREND_SCENARIO,
    delta_inc = DELTA_INC,
    anchor_pseudo_w = ANCHOR_PSEUDO_W,
    sd_beta_I = SD_BETA_I,
    use_weighted_cohort = USE_WEIGHTED_COHORT,
    sd_cohort_resid = SD_COHORT_RESID, sd_beta_fixed = SD_BETA_FIXED,
    use_age_slope = FALSE,
    beta_force = NULL,
    gammaP_method = GAMMAP_METHOD, trend_type = TREND_TYPE,
    path_prev_dta = PATH_PREV_DTA,
    prev_micro_df = NULL,
    prev_cfg = NULL,
    L_I_max_years = L_I_MAX_YEARS,
    mort_period_shock_years = integer(0),
    mort_downweight_years = integer(0),
    mort_downweight_weight = MORT_DOWNWEIGHT_WEIGHT_F,
    mort_hist_tbl = mort_hist, pop_all_tbl = pop_all,
    inc_hist_tbl  = NULL,
    cause_id_override = NA_character_,
    rr_inc = NA_real_,
    rr_mort = NA_real_,
    sd_theta_IP = SD_THETA_IP
){
  if (is.null(prev_cfg)) prev_cfg <- make_prev_config()
  if (!"cause" %in% names(mort_hist_tbl)) {
    rlang::abort("run_pipeline_sex(): 'mort_hist_tbl' must contain a `cause` column. Load cause-specific histories.")
  }
  sex_sel             <- match.arg(sex_sel)
  beta_mode           <- "fixed_rr_offset"
  gammaP_method       <- match.arg(gammaP_method, c("freeze","arima","trend","damped_trend"))
  mort_trend_scenario <- match.arg(tech_scenario, c("freeze","continue","delta"))
  mort_bapc_trend_scenario <- match.arg(mort_bapc_trend_scenario, c("freeze","continue","delta"))
  params <- list(rr_inc = rr_inc, rr_mort = rr_mort) 
  
  # ---------- 1) Historia mortalidad & población (ventana)
  .bapc_verbose("run_pipeline_sex: mort_hist_tbl is NULL? ", is.null(mort_hist_tbl))
  .bapc_verbose("run_pipeline_sex: pop_all_tbl is NULL? ", is.null(pop_all_tbl))
  .bapc_verbose("run_pipeline_sex: inc_hist_tbl is NULL? ", is.null(inc_hist_tbl))
  mortH <- mort_hist_tbl %>% dplyr::filter(
    sex == sex_sel, age >= age_min_m, age <= age_max_m,
    period >= period_min_m, period <= period_max_m
  )
  last_hist_year <- max(mortH$period)

  pop_future <- pop_all_tbl %>% dplyr::filter(
    sex == sex_sel, period >= last_hist_year + 1,
    age >= age_min_m, age <= age_max_m
  )
  
  mortH <- mortH %>%
    dplyr::left_join(pop_all_tbl, by = c("period","age","sex"), suffix = c("",".pop")) %>%
    dplyr::mutate(exposure = dplyr::coalesce(exposure, exposure.pop)) %>%
    dplyr::select(-exposure.pop)
  
  # ---------- 2) PREVALENCIA APC 
  prev_agg <- if (!is.null(prev_micro_df)) {
    build_prev_from_micro_df(
      micro_df = prev_micro_df, sex_sel = sex_sel,
      period_min = period_min_p, period_max = period_max_p,
      age_min = age_min_p, age_max = age_max_p, min_neff = 3
    )
  } else {
    build_prev_from_micro(
      path_dta = path_prev_dta, sex_sel = sex_sel,
      period_min = period_min_p, period_max = period_max_p,
      age_min = age_min_p, age_max = age_max_p, min_neff = 3
    )
  }
  # light diagnostic of PREV inputs is omitted here to avoid cluttering the console in simulations
  tab_inst <- table(prev_agg$inst); keep_levels <- names(tab_inst)[tab_inst > 2]
  prev_agg <- prev_agg %>% dplyr::filter(inst %in% keep_levels) %>% droplevels()
  
  prev_inla <- prev_agg %>%
    dplyr::rename(y = y_eff, N = neff) %>%
    dplyr::mutate(y = as.integer(y), N = as.integer(N)) %>%
    dplyr::filter(is.finite(y), is.finite(N), N > 0L, y >= 0L, y <= N) %>%
    dplyr::mutate(
      age_id = match(age, sort(unique(age))),
      period_id = match(period, sort(unique(period))),
      cohort_id = match(cohort, sort(unique(cohort)))
    )
  if ("inst" %in% names(prev_inla)) {
    one_level_factors <- names(Filter(function(x) is.factor(x) && nlevels(x) < 2, prev_inla))
    if (length(one_level_factors)) prev_inla[one_level_factors] <- NULL
  }
  
  prev_inla <- attach_edge_weights_hist(prev_inla, stage = "prev")
  
  mu_perP <- mean(sort(unique(prev_inla$period)))
  prev_inla <- prev_inla %>%
    dplyr::bind_cols(make_trend_vars(.$period, mu_perP, PREV_TREND_DEGREE, prefix = "prev_"))
  
  ctrl_fix_prev <- if (PREV_TREND_DEGREE > 0) {
    list(
      mean = list(prev_trend_t = 0),
      prec = list(prev_trend_t = 1 / PREV_TREND_PRIOR_SD^2)
    )
  } else {
    list()
  }
  
  hyper_age_prev <- pc_hyper(PREV_AGE_PC_U, PREV_AGE_PC_A)
  hyper_per_prev <- pc_hyper(PREV_PER_PC_U, PREV_PER_PC_A)
  hyper_coh_prev <- pc_hyper(PREV_COH_PC_U, PREV_COH_PC_A)
  
  use_inst <- "inst" %in% names(prev_inla)
  if (use_inst) prev_inla$inst <- factor(prev_inla$inst) 
  
  if (use_inst) {
    form_prev <- if (PREV_TREND_DEGREE > 0) {
      y ~ -1 + inst + prev_trend_t +
        f(age_id,    model = PREV_AGE_MODEL, constr=TRUE, scale.model=TRUE, hyper=hyper_age_prev) +
        f(period_id, model = PREV_PER_MODEL, constr=TRUE, scale.model=TRUE, hyper=hyper_per_prev) +
        f(cohort_id, model = PREV_COH_MODEL, constr=TRUE, scale.model=TRUE, hyper=hyper_coh_prev)
    } else {
      y ~ -1 + inst +
        f(age_id,    model = PREV_AGE_MODEL, constr=TRUE, scale.model=TRUE, hyper=hyper_age_prev) +
        f(period_id, model = PREV_PER_MODEL, constr=TRUE, scale.model=TRUE, hyper=hyper_per_prev) +
        f(cohort_id, model = PREV_COH_MODEL, constr=TRUE, scale.model=TRUE, hyper=hyper_coh_prev)
    }
  } else {
    form_prev <- if (PREV_TREND_DEGREE > 0) {
      y ~ 1 + prev_trend_t +
        f(age_id,    model = PREV_AGE_MODEL, constr=TRUE, scale.model=TRUE, hyper=hyper_age_prev) +
        f(period_id, model = PREV_PER_MODEL, constr=TRUE, scale.model=TRUE, hyper=hyper_per_prev) +
        f(cohort_id, model = PREV_COH_MODEL, constr=TRUE, scale.model=TRUE, hyper=hyper_coh_prev)
    } else {
      y ~ 1 +
        f(age_id,    model = PREV_AGE_MODEL, constr=TRUE, scale.model=TRUE, hyper=hyper_age_prev) +
        f(period_id, model = PREV_PER_MODEL, constr=TRUE, scale.model=TRUE, hyper=hyper_per_prev) +
        f(cohort_id, model = PREV_COH_MODEL, constr=TRUE, scale.model=TRUE, hyper=hyper_coh_prev)
    }
  }
  
  .chk("pre PREV try")
  .check_no_mode_in_f(form_prev)
  fit_prev <- try(
    inla_tag("PREV", formula = form_prev, family="binomial", Ntrials=prev_inla$N, data=prev_inla,
    weights = prev_inla$edge_weight,
    control.fixed=ctrl_fix_prev,
    control.family=list(link="logit"),
    control.predictor=list(compute=TRUE),
    control.compute=list(dic=TRUE, waic=TRUE, cpo=TRUE, config = TRUE),
    control.inla=list(diagonal=1e-5)
  ), silent = TRUE)
  if (inherits(fit_prev, "try-error")) {
    form_prev <- update(form_prev, . ~ . - inst)
    .chk("pre PREV-2")
    .check_no_mode_in_f(form_prev)
    fit_prev <- inla_tag("PREV-2", formula = form_prev, family="binomial", Ntrials=prev_inla$N, data=prev_inla,
      weights = prev_inla$edge_weight,
      control.fixed=ctrl_fix_prev,
      control.family=list(link="logit"),
      control.predictor=list(compute=TRUE),
      control.compute=list(dic=TRUE, waic=TRUE, cpo=TRUE, config = TRUE),
      control.inla=list(diagonal=1e-5)
    )
  }
  
  scP <- fit_prev$summary.random$cohort_id
  lev_coh_prev <- sort(unique(prev_inla$cohort))
  gammaP_hist <- tibble::tibble(cohort = lev_coh_prev[scP$ID],
                                gammaP = scP$mean - mean(scP$mean)) %>% dplyr::arrange(cohort)
  coh_all_needed <- sort(unique((pop_all_tbl$period - pop_all_tbl$age)))
  gammaP_fut <- forecast_gammaP(
    gammaP_hist, setdiff(coh_all_needed, gammaP_hist$cohort),
    method = gammaP_method, trend_type = trend_type
  )
  gammaP_fut <- adjust_gammaP_future(gammaP_hist, gammaP_fut,
                                     scenario = if (!is.null(prev_cfg$backbone) && identical(prev_cfg$backbone, "forecast")) prev_cfg$scenario else "freeze",
                                     annual_rate = prev_cfg$annual_rate,
                                     annual_rate_up = prev_cfg$annual_rate_up %||% prev_cfg$annual_rate,
                                     annual_rate_down = prev_cfg$annual_rate_down %||% prev_cfg$annual_rate,
                                     annual_rate_down3 = prev_cfg$annual_rate_down3,
                                     base_year = prev_cfg$base_year)
  gammaP_all <- dplyr::bind_rows(gammaP_hist, gammaP_fut)
  
  # ---------- 3) Pure BAPC INCIDENCE (first), to extract historical period effect per(·)
  if (is.null(inc_hist_tbl)) inc_hist_tbl <- load_incidence_lung(PATH_INC_CSV, pop_all_tbl)
  inc_fit_bapc <- fit_apc_incidence(
    inc_hist = inc_hist_tbl %>% dplyr::filter(sex == sex_sel),
    pop_all  = pop_all_tbl %>% dplyr::filter(sex == sex_sel),
    age_min_i = age_min_i, age_max_i = age_max_i,       # <- FIX: usar edades de INCIDENCIA
    period_min_i = 1998, period_max_i = 2022,
    use_weighted_cohort = use_weighted_cohort,
    inc_scenario = inc_trend_scenario,
    delta_inc   = delta_inc
  )
  
  # ---------- 3.b) Unique operational policy PREV->INC
  engine_require_stock_former_policy(
    beta_mode = beta_mode,
    prev_cfg = prev_cfg,
    prev_inc_channel_mode = PREV_INC_CHANNEL_MODE
  )
  use_stock_channel_main <- TRUE
  prev_sign_ <- 1L
  cause_id_cur <- suppressWarnings(as.character(cause_id_override)[1])
  if (!length(cause_id_cur) || is.na(cause_id_cur) || !nzchar(cause_id_cur)) {
    cause_id_cur <- suppressWarnings(as.character(unique(mortH$cause))[1])
  }
  if (!length(cause_id_cur) || is.na(cause_id_cur) || !nzchar(cause_id_cur)) {
    cause_id_cur <- tryCatch(suppressWarnings(as.character(unique(inc_hist_tbl$causa)[1])), error = function(e) NA_character_)
  }
  if (!length(cause_id_cur) || is.na(cause_id_cur) || !nzchar(cause_id_cur)) {
    cause_id_cur <- tryCatch(suppressWarnings(as.character(unique(inc_hist_tbl$cause)[1])), error = function(e) NA_character_)
  }
  cause_id_cur <- tryCatch(normalize_cause_id(cause_id_cur), error = function(e) as.character(cause_id_cur)[1])

  if (!exists("make_engine_method_policy", inherits = TRUE)) {
    make_engine_method_policy <- function(...) {
      list(...)
    }
  }
  if (!exists("engine_method_policy_table", inherits = TRUE)) {
    engine_method_policy_table <- function(policy) {
      if (is.null(policy) || !length(policy)) {
        return(tibble::tibble(key = character(), value = character()))
      }
      tibble::tibble(
        key = names(policy),
        value = vapply(policy, function(v) {
          if (length(v) == 1 && (is.atomic(v) || is.null(v))) {
            as.character(v)
          } else {
            paste(utils::capture.output(str(v, max.level = 1, give.attr = FALSE)), collapse = " ")
          }
        }, character(1))
      )
    }
  }
  method_policy <- make_engine_method_policy(
    sex = sex_sel,
    age_min_p = age_min_p, age_max_p = age_max_p,
    age_min_i = age_min_i, age_max_i = age_max_i,
    age_min_m = age_min_m, age_max_m = age_max_m,
    prev_cfg = prev_cfg,
    gammaP_method = gammaP_method,
    trend_type = trend_type,
    use_weighted_cohort = use_weighted_cohort,
    inc_include_trend = inc_include_trend,
    inc_trend_scenario = inc_trend_scenario,
    mort_trend_scenario = mort_trend_scenario,
    mort_bapc_trend_scenario = mort_bapc_trend_scenario,
    beta_mode = beta_mode,
    cause_id = cause_id_cur
  )
  method_policy_tbl <- engine_method_policy_table(method_policy)
  rr_inc_cur <- if (!is.null(rr_inc) && is.finite(rr_inc)) as.numeric(rr_inc) else get_inc_rr_by_cause_sex(cause_id_cur, sex_sel)
  prev_base_prob <- prev_agg %>%
    dplyr::filter(period == (prev_cfg$base_year %||% PREV_BASE_YEAR)) %>%
    dplyr::summarise(p = sum(y_eff, na.rm = TRUE) / pmax(sum(neff, na.rm = TRUE), 1e-12), .groups = "drop") %>%
    dplyr::pull(p)
  if (!length(prev_base_prob) || !is.finite(prev_base_prob)) {
    prev_base_prob <- if (identical(sex_sel, "M")) prev_cfg$prev_base_M %||% prev_cfg$prev_base_default else prev_cfg$prev_base_F %||% prev_cfg$prev_base_default
  }
  .bapc_verbose(sprintf("[%s] PREV->INC channel=stock_former | beta_mode=%s | RR_I=%.3f | prev_base_prob=%.4f", sex_sel, beta_mode, rr_inc_cur, prev_base_prob))
  
  # ---------- 3.c) INCIDENCE conditional on PREV
  inc_fit <- tryCatch({
    fit_apc_incidence_cond_prev(
      inc_hist = inc_hist_tbl %>% dplyr::filter(sex == sex_sel),
      pop_all  = pop_all_tbl %>% dplyr::filter(sex == sex_sel),
      fit_prev = fit_prev,
      age_min_i = age_min_i, age_max_i = age_max_i,       
      period_min_i = 1998, period_max_i = 2022,
      A_I = NA_integer_, w_I = PREV_W_I,                 
      use_weighted_cohort = use_weighted_cohort,
      sd_theta_IP = sd_theta_IP,
      gammaP_method = gammaP_method, trend_type = trend_type,
      inc_trend_scenario = if (isTRUE(inc_include_trend)) inc_trend_scenario else "none",
      delta_inc = delta_inc,
      sd_beta_I = sd_beta_I,
      prev_sign = prev_sign_,
      prev_scenario = prev_cfg$scenario,
      prev_scenario_axis = prev_cfg$axis,
      prev_annual_rate = prev_cfg$annual_rate,
      prev_annual_rate_up = prev_cfg$annual_rate_up,
      prev_annual_rate_down = prev_cfg$annual_rate_down,
      prev_annual_rate_down3 = prev_cfg$annual_rate_down3,
      prev_base_year = prev_cfg$base_year,
      quit_mode = prev_cfg$quit_mode,
      quit_floor_sd = prev_cfg$quit_floor_sd,
      quit_floor_sd_M = prev_cfg$quit_floor_sd_M,
      quit_floor_sd_F = prev_cfg$quit_floor_sd_F,
      quit_half_life = prev_cfg$quit_half_life,
      quit_ramp_years = prev_cfg$quit_ramp_years,
      prev_base_M = prev_cfg$prev_base_M,
      prev_base_F = prev_cfg$prev_base_F,
      prev_base_default = prev_cfg$prev_base_default,
      rr_inc = rr_inc_cur,
      prev_base_prob = prev_base_prob,
      cause_id = cause_id_cur
    )
  }, error = function(e) {
    message(sprintf("[run_pipeline_sex] CRITICAL ERROR in incidence for %s: %s", sex_sel, conditionMessage(e)))
    NULL
  })
  
  if (is.null(inc_fit) || is.null(inc_fit$rates_all)) {
    warning(sprintf("[run_pipeline_sex] fit_apc_incidence_cond_prev returned NULL for %s. Aborting this sex safely.", sex_sel))
    return(NULL)
  }

  betaP_raw <- tryCatch(as.numeric(inc_fit$beta_P), error = function(e) NA_real_)
  beta_post <- .beta_postfit_transform(betaP_raw, rule = BETA_P_POSTFIT_RULE)
  betaP_eff <- dplyr::coalesce(tryCatch(inc_fit$beta_P_eff, error = function(e) NA_real_), beta_post$eff, 0)
  betaP_zeroed <- dplyr::coalesce(tryCatch(inc_fit$beta_P_zeroed, error = function(e) NA), beta_post$zeroed, FALSE)
  # VERY IMPORTANT: DO NOT modify prev_sign_ here. Retain the sign from the calibration stage (cal$sign).
  if (identical(beta_mode, "fixed_rr_offset")) {
    .bapc_verbose(sprintf("[run_pipeline_sex] beta_mode=fixed_rr_offset | RR_I=%.4f | prev_sign_=%+d", rr_inc_cur, prev_sign_))
  } else {
    .bapc_verbose(sprintf("[run_pipeline_sex] beta_P (raw)=%.4f -> beta_P_eff=%.4f | rule=%s | zeroed=%s | prev_sign_=%+d",
                    betaP_raw, betaP_eff, BETA_P_POSTFIT_RULE, as.character(betaP_zeroed), prev_sign_))
  }
  
  inc_rates_hist <- inc_fit$rates_all %>% dplyr::filter(period <= last_hist_year)
  
  # --- Save actually used A (P -> I) ---
  A_I_used <- NA_integer_                            # stock_former: A_I is not used
  params$A_I <- A_I_used
  .expose_selected(sex_sel, A_I_used = A_I_used)
  
  # ---- Join rates with exposure for expected counts (conditional on PREV)
  inc_withE <- inc_fit$rates_all %>%
    dplyr::left_join(
      pop_all_tbl %>% dplyr::filter(sex == sex_sel) %>%
        dplyr::select(age, period, sex, exposure),
      by = c("sex","age","period")
    ) %>%
    dplyr::mutate(
      exposure  = pmax(as.numeric(exposure), 1e-12),
      cases_hat = rate_hat * exposure,
      cases_lwr = dplyr::coalesce(rate_lwr, rate_hat) * exposure,
      cases_upr = dplyr::coalesce(rate_upr, rate_hat) * exposure
    )
  
  # ---- Annual totals (conditional on PREV)
  inc_annual <- inc_withE %>%
    dplyr::group_by(period) %>%
    dplyr::summarise(
      cases_hat = sum(cases_hat, na.rm = TRUE),
      cases_lwr = sum(cases_lwr, na.rm = TRUE),
      cases_upr = sum(cases_upr, na.rm = TRUE),
      .groups = "drop"
    )
  
  # ---- Observed historical totals
  inc_obs_annual <- inc_hist_tbl %>% dplyr::filter(sex == sex_sel) %>%
    dplyr::group_by(period) %>%
    dplyr::summarise(obs = sum(cases, na.rm = TRUE), .groups = "drop")
  
  # ---- Helper for BAPC/Conditional annual totals (used below for both)
  .as_annual_inc <- function(rates_all, pop_tbl) {
    rates_all %>%
      dplyr::left_join(pop_tbl %>% dplyr::select(age, period, sex, exposure),
                       by = c("age","period","sex")) %>%
      dplyr::mutate(
        cases_hat = rate_hat * exposure,
        cases_lwr = rate_lwr * exposure,
        cases_upr = rate_upr * exposure
      ) %>%
      dplyr::group_by(period) %>%
      dplyr::summarise(cases_hat = sum(cases_hat),
                       cases_lwr = sum(cases_lwr),
                       cases_upr = sum(cases_upr), .groups="drop")
  }
  
  # === Pure I (BAPC without prevalence channel) ===
  inc_annual_bapc <- .as_annual_inc(inc_fit_bapc$rates_all, pop_all_tbl %>% dplyr::filter(sex == sex_sel))

  # === I|P (main conditional) ===
  inc_annual_cond <- .as_annual_inc(inc_fit$rates_all, pop_all_tbl %>% dplyr::filter(sex == sex_sel))
  
  # === Counterfactual I without P ===
  # In fixed_rr_offset, remove the epidemiological offset exp(log1p(q_eff*(RR_I-1))).
  # In older modes, remove exp(beta_P * z_prev).
  z_grid <- tryCatch({
    dplyr::bind_rows(
      tibble::as_tibble(inc_fit$z_hist %||% tibble::tibble()),
      tibble::as_tibble(inc_fit$z_future %||% tibble::tibble())
    ) %>%
      dplyr::select(dplyr::any_of(c("sex","age","period","q_eff","z_prev"))) %>%
      dplyr::distinct()
  }, error = function(e) tibble::tibble())
  if (!nrow(z_grid)) {
    z_grid <- tidyr::expand_grid(
      age    = seq.int(age_min_i, age_max_i),
      period = inc_fit$lev_inc$period
    ) %>%
      dplyr::mutate(sex = factor(sex_sel, levels = c("M","F")))

    if (identical(beta_mode, "fixed_rr_offset")) {
      z_grid <- build_prev_rr_offset_for_inc(
        df_inc_grid = z_grid,
        gammaP_all = gammaP_all,
        A_I = NA_integer_,
        rr_inc = rr_inc_cur,
        prev_base_prob = prev_base_prob,
        base_year = prev_cfg$base_year
      ) %>%
        dplyr::select(sex, age, period, q_eff, z_prev)
    } else {
      # Fallback legacy (without build_prev_index_for_inc, which has been removed)
      z_grid <- z_grid %>% dplyr::mutate(z_prev = 0, q_eff = 0)
    }
  }
  
  betaP <- betaP_eff
  inc_rates_noP <- inc_fit$rates_all %>%
    dplyr::left_join(z_grid, by = c("sex","age","period")) %>%
    dplyr::mutate(
      adj = if (identical(beta_mode, "fixed_rr_offset")) {
        exp(dplyr::coalesce(z_prev, 0))
      } else {
        exp(betaP * dplyr::coalesce(z_prev, 0))
      },
      rate_hat = pmax(rate_hat / pmax(adj, 1e-12), 1e-12),
      rate_lwr = ifelse(is.finite(rate_lwr), pmax(rate_lwr / pmax(adj, 1e-12), 1e-12), NA_real_),
      rate_upr = ifelse(is.finite(rate_upr), pmax(rate_upr / pmax(adj, 1e-12), 1e-12), NA_real_)
    ) %>%
    dplyr::select(sex, age, period, rate_hat, rate_lwr, rate_upr)
  
  inc_annual_noP <- .as_annual_inc(inc_rates_noP, pop_all_tbl %>% dplyr::filter(sex == sex_sel))
  
  # ---------- 4) Mortality grid historical + future (including exposure)
  future_grid <- pop_future %>%
    dplyr::mutate(
      cohort = period - age,
      sex = factor(sex_sel, levels = c("M","F")),
      cause = unique(mortH$cause)[1],
      deaths = NA_real_
    ) %>%
    dplyr::select(age, period, sex, cause, deaths, exposure, cohort)
  
  mort_all <- dplyr::bind_rows(
    mortH %>% dplyr::select(age, period, sex, cause, deaths, exposure, cohort),
    future_grid
  )
  check_apc_grid(mort_all)
  
  if (!"exposure" %in% names(mort_all)) {
    exposure_cols <- intersect(c("exposure","exposure.x","exposure.y","poblacion","E"), names(mort_all))
    mort_all <- mort_all %>%
      dplyr::mutate(exposure = dplyr::coalesce(!!!rlang::syms(exposure_cols))) %>%
      dplyr::select(-dplyr::any_of(setdiff(exposure_cols, "exposure")))
  }
  stopifnot("exposure" %in% names(mort_all))
  
  mort_all <- mort_all %>%
    dplyr::mutate(
      y    = as.integer(round(deaths)),
      E    = pmax(as.numeric(exposure), 1e-12),
      logE = log(E)
    )
  
  # levels/indices and trend terms
  lev_age <- sort(unique(mortH$age))
  lev_per <- sort(unique(mortH$period))
  lev_coh <- sort(unique(mortH$cohort))
  mu_perM <- mean(mortH$period)
  
  Xconstr_per_M <- make_slope_constr(lev_per, mu_perM)
  
  mort_all <- mort_all %>%
    dplyr::mutate(
      age_id    = match(age, lev_age),
      period_id = match(pmin(period, max(lev_per)), lev_per),
      cohort_true = cohort,
      cohort_clamped = pmin(pmax(cohort_true, min(lev_coh)), max(lev_coh)),
      cohort_ref_anchor = cohort_clamped,
      cohort_id_bapc   = match(cohort_clamped, lev_coh),
      cohort_id_anchor = match(cohort_ref_anchor, lev_coh),
      is_future = period > last_hist_year
    ) %>%
    dplyr::bind_cols(make_trend_vars(.$period, mu_perM, MORT_TREND_DEGREE, prefix = "mort_")) %>%
    dplyr::rename(mort_anchor_trend_t = mort_trend_t,
                  mort_anchor_trend_t2 = mort_trend_t2) %>%
    dplyr::mutate(mort_anchor_tech_offset = 0) %>%
    dplyr::bind_cols(make_trend_vars(.$period, mu_perM, MORT_TREND_DEGREE, prefix = "mort_bapc_")) %>%
    dplyr::mutate(mort_bapc_tech_offset = 0)
  
  # apply trend scenarios in future
  idxF <- which(mort_all$is_future)
  if (length(idxF)) {
    fut_terms_anchor <- apply_trend_scenario_future(
      period   = mort_all$period[idxF],
      last_hist= last_hist_year, center = mu_perM,
      degree   = MORT_TREND_DEGREE,
      scenario = mort_trend_scenario,
      delta    = delta_tech, prefix = "mort_anchor_"
    )
    mort_all$mort_anchor_trend_t[idxF]     <- fut_terms_anchor$mort_anchor_trend_t
    mort_all$mort_anchor_trend_t2[idxF]    <- fut_terms_anchor$mort_anchor_trend_t2
    mort_all$mort_anchor_tech_offset[idxF] <- fut_terms_anchor$mort_anchor_tech_offset

    fut_terms_bapc <- apply_trend_scenario_future(
      period   = mort_all$period[idxF],
      last_hist= last_hist_year, center = mu_perM,
      degree   = MORT_TREND_DEGREE,
      scenario = mort_bapc_trend_scenario,
      delta    = delta_tech, prefix = "mort_bapc_"
    )
    mort_all$mort_bapc_trend_t[idxF]     <- fut_terms_bapc$mort_bapc_trend_t
    mort_all$mort_bapc_trend_t2[idxF]    <- fut_terms_bapc$mort_bapc_trend_t2
    mort_all$mort_bapc_tech_offset[idxF] <- fut_terms_bapc$mort_bapc_tech_offset
  }
  
  # shocks de período (solo historia)
  # 1) años shock definidos a nivel global/por causa
  shock_vals <- sort(unique(as.integer(mort_period_shock_years)))
  
  # 2) quedarnos sólo con los que realmente existen en los datos históricos
  if (length(shock_vals) > 0) {
    shock_vals_in_data <- intersect(shock_vals, unique(mort_all$period[!mort_all$is_future]))
  } else {
    shock_vals_in_data <- integer(0)
  }
  
  # 3) construir shock_id solo si hay shocks válidos en datos
  if (length(shock_vals_in_data) > 0) {
    mort_all <- mort_all %>%
      dplyr::mutate(
        shock_id = ifelse(!is_future & period %in% shock_vals_in_data,
                          match(period, shock_vals_in_data), NA_integer_)
      )
  } else {
    # sin shocks: aseguramos la columna (o la retiramos más abajo)
    mort_all <- mort_all %>%
      dplyr::mutate(shock_id = NA_integer_)
  }
  
  
  # --- Nuevo canal externo INC -> MORT ---
  # Ya no se calibra un único L_I: la alineación temporal la trae el kernel
  # distribuido por años post-diagnóstico.
  L_I_eff <- NA_integer_

  params$L_I <- NA_integer_
  params$bridge_years <- 0L
  params$mort_link_mode <- MORT_I_LINK_MODE
  params$mort_bapc_trend_scenario <- mort_bapc_trend_scenario
  params$mort_bapc_future_mode <- "autonomous_apc"
  .expose_selected(sex_sel, A_I_used = NA_integer_, bridge_years = params$bridge_years)

  # Grid base de mortalidad para el BAPC puro (sin enlace impuesto)
  df_all <- mort_all %>%
    dplyr::mutate(
      y    = ifelse(is.na(deaths), NA_integer_, as.integer(round(deaths))),
      E    = pmax(as.numeric(exposure), 1e-12),
      logE = log(E),
      cohort_true = cohort,
      cohort_clamped = pmin(pmax(cohort_true, min(lev_coh)), max(lev_coh)),
      cohort_ref_anchor = cohort_clamped,
      cohort_id_bapc   = match(cohort_clamped, lev_coh),
      cohort_id_anchor = match(cohort_ref_anchor, lev_coh)
    )
  
  # ---- Priors BAPC mortalidad
  hyper_age      <- pc_hyper(MORT_AGE_PC_U, MORT_AGE_PC_A)
  hyper_per      <- pc_hyper(MORT_PER_PC_U, MORT_PER_PC_A)
  hyper_coh_bapc <- pc_hyper(MORT_COH_PC_U, MORT_COH_PC_A)
  
  # ---------- BAPC mortalidad ----------
  if (use_weighted_cohort) {
    cov_coh_m <- cohort_coverage_weights(mortH) %>% dplyr::right_join(tibble::tibble(cohort = lev_coh), by="cohort") %>%
      dplyr::mutate(w = dplyr::coalesce(w, min(w, na.rm = TRUE)))
    Qw_m <- make_Q_rw1_weighted(length(lev_coh), cov_coh_m$w)
    form_bapc <- y ~ 1 + mort_bapc_trend_t + 
      f(age_id,    model = MORT_AGE_MODEL, constr = TRUE, scale.model = TRUE, hyper = hyper_age) +
      f(period_id, model = MORT_PER_MODEL, constr = TRUE, extraconstr = Xconstr_per_M, scale.model = TRUE, hyper = hyper_per) +
      f(cohort_id_bapc, model = "generic0", Cmatrix = Qw_m, constr = TRUE, hyper = hyper_coh_bapc) +
      offset(logE + mort_bapc_tech_offset)
    if (any(is.finite(mort_all$shock_id))) {
      form_bapc <- update(
        form_bapc,
        . ~ . + f(shock_id, model = "iid",
                  hyper = pc_hyper(MORT_SHOCK_PC_U, MORT_SHOCK_PC_A))
      )
    }
  } else {
    form_bapc <- y ~ 1 + mort_bapc_trend_t + 
      f(age_id,    model = MORT_AGE_MODEL, constr = TRUE, scale.model = TRUE, hyper = hyper_age) +
      f(period_id, model = MORT_PER_MODEL, constr = TRUE, extraconstr = Xconstr_per_M, scale.model = TRUE, hyper = hyper_per) +
      f(cohort_id_bapc, model = MORT_COH_MODEL, constr = TRUE, scale.model = TRUE, hyper = hyper_coh_bapc) +
      offset(logE + mort_bapc_tech_offset)
    if (any(is.finite(mort_all$shock_id))) {
      form_bapc <- update(
        form_bapc,
        . ~ . + f(shock_id, model = "iid",
                  hyper = pc_hyper(MORT_SHOCK_PC_U, MORT_SHOCK_PC_A))
      )
    }
  }
  trend_mean <- if (sex_sel == "F") MORT_TREND_PRIOR_MEAN_F else MORT_TREND_PRIOR_MEAN_M
  trend_sd   <- if (sex_sel == "F") MORT_TREND_PRIOR_SD_F   else MORT_TREND_PRIOR_SD_M
  ctrl_fixed_bapc <- list(
    mean = list(mort_bapc_trend_t = trend_mean),
    prec = list(
      mort_bapc_trend_t  = 1 / (trend_sd^2)
    )
  )
  df_bapc <- df_all %>% dplyr::mutate(cohort_id = cohort_id_bapc)
  df_bapc <- attach_edge_weights_hist(df_bapc, stage = "mort")
  
  
  # Pesos por observación (default = 1)
  w_obs_bapc <- rep(1, nrow(df_bapc))
  if (length(mort_downweight_years) > 0) {
    idx_dw <- df_bapc$period %in% mort_downweight_years
    w_obs_bapc[idx_dw] <- mort_downweight_weight
  }
  w_fit_bapc <- as.numeric(w_obs_bapc) * dplyr::coalesce(df_bapc$edge_weight, 1)
  
  .check_no_mode_in_f(form_bapc)
  fit_bapc <- inla_tag("MORT-BAPC", formula = form_bapc, family = "poisson",
                       data = df_bapc,
                       weights = w_fit_bapc,
                       control.fixed     = ctrl_fixed_bapc,
                       control.predictor = list(compute = TRUE),
                       control.compute   = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
                       control.inla      = list(diagonal = 1e-5)
  )
  
  if (sex_sel == "F") {
    yrs_dw <- sort(unique(df_bapc$period[w_obs_bapc != 1]))
    if (length(yrs_dw)) {
      .bapc_verbose(paste0("Female downweight applied in BAPC for years: ", paste(yrs_dw, collapse = ", ")))
    } else {
      .bapc_verbose("Female downweight applied in BAPC: none")
    }
  }
  
  # ===== Incidencia para el enlace a mortalidad =====
  # El canal principal ahora transforma incidencia en muertes esperadas
  # mediante un kernel externo distribuido por años post-diagnóstico.
  mort_Icond <- attach_external_mortality_offset(
    mort_all = df_all,
    inc_rates_all = inc_fit$rates_all,
    pop_all_tbl = pop_all_tbl %>% dplyr::filter(.data$sex == sex_sel),
    cause_id = cause_id_cur,
    sex_sel = sex_sel
  )

  # Contrafactual I sin P
  inc_noP_rates_all <- inc_rates_noP %>%
    dplyr::select(sex, age, period, rate_hat)

  mort_InoP <- attach_external_mortality_offset(
    mort_all = df_all,
    inc_rates_all = inc_noP_rates_all,
    pop_all_tbl = pop_all_tbl %>% dplyr::filter(.data$sex == sex_sel),
    cause_id = cause_id_cur,
    sex_sel = sex_sel
  )

  kernel_tbl <- attr(mort_Icond, "kernel") %||% tibble::tibble()
  kernel_summary_cond <- attr(mort_Icond, "kernel_summary") %||% tibble::tibble()
  kernel_summary_noP <- attr(mort_InoP, "kernel_summary") %||% tibble::tibble()

  params$mort_kernel_max_lag <- if (nrow(kernel_tbl)) max(kernel_tbl$lag, na.rm = TRUE) else NA_integer_
  params$mort_kernel_total_mass <- if (nrow(kernel_tbl)) sum(kernel_tbl$weight, na.rm = TRUE) else NA_real_

  # ===== Formulación común del modelo anclado =====
  hyper_coh_anchor <- pc_hyper(MORT_COH_PC_U, MORT_COH_PC_A)
  
  trend_mean <- if (sex_sel == "F") MORT_TREND_PRIOR_MEAN_F else MORT_TREND_PRIOR_MEAN_M
  trend_sd   <- if (sex_sel == "F") MORT_TREND_PRIOR_SD_F   else MORT_TREND_PRIOR_SD_M

  form_anchor_core <- y ~ 1 + mort_anchor_trend_t + 
    f(age_id,    model = MORT_AGE_MODEL, constr = TRUE, scale.model = TRUE, hyper = hyper_age) +
    f(period_id, model = MORT_PER_MODEL, constr = TRUE, extraconstr = Xconstr_per_M, scale.model = TRUE, hyper = hyper_per) +
    offset(log_mort_ext + mort_anchor_tech_offset)

  ctrl_fixed_anchor <- list(
    mean = list(mort_anchor_trend_t = trend_mean),
    prec = list(mort_anchor_trend_t = 1 / (trend_sd^2))
  )
  
  if (use_weighted_cohort) {
    form_anchor <- update(form_anchor_core, . ~ . +
                            f(cohort_id_anchor, model = "generic0", Cmatrix = Qw_m, constr = TRUE, hyper = hyper_coh_anchor))
  } else {
    form_anchor <- update(form_anchor_core, . ~ . +
                            f(cohort_id_anchor, model = MORT_COH_MODEL, constr = TRUE, scale.model = TRUE, hyper = hyper_coh_anchor))
  }
  
  if (any(is.finite(mort_all$shock_id))) {
    form_anchor <- update(form_anchor, . ~ . + f(shock_id, model = "iid",
                                                 hyper = pc_hyper(MORT_SHOCK_PC_U, MORT_SHOCK_PC_A)))
  }
  
  
  # ===== Pseudo-observaciones de continuidad (último hist → primer futuro) =====
#  if (bridge_inc_years > 0) { 
#    w_pseudo <- anchor_pseudo_w
#    
#    pseudo_rows_cond <- mort_Icond %>%
#      dplyr::filter(period == last_hist_year) %>%
#      dplyr::mutate(period = last_hist_year + 1, y = y)
#    mort_Icond_w_pseudo <- mort_Icond %>%
#      dplyr::filter(period != last_hist_year + 1) %>%
#      dplyr::bind_rows(pseudo_rows_cond)
#    pesos_cond <- c(rep(1, nrow(mort_Icond_w_pseudo) - nrow(pseudo_rows_cond)),
#                    rep(w_pseudo, nrow(pseudo_rows_cond)))
#    
#    pseudo_rows_noP <- mort_InoP %>%
#      dplyr::filter(period == last_hist_year) %>%
#      dplyr::mutate(period = last_hist_year + 1, y = y)
#    mort_InoP_w_pseudo <- mort_InoP %>%
#      dplyr::filter(period != last_hist_year + 1) %>%
#      dplyr::bind_rows(pseudo_rows_noP)
#    pesos_noP <- c(rep(1, nrow(mort_InoP_w_pseudo) - nrow(pseudo_rows_noP)),
#                   rep(w_pseudo, nrow(pseudo_rows_noP)))
#  }
  
  # ========= Elegir data robusta =========
  data_cond <- mort_Icond
  data_noP  <- mort_InoP
  
  scale_cond <- if (exists("pesos_cond", inherits = FALSE)) pesos_cond else rep(1, nrow(data_cond))
  scale_noP  <- if (exists("pesos_noP",  inherits = FALSE)) pesos_noP  else rep(1, nrow(data_noP))
  
  # ========= Sanity: futuras como NA; offsets finitos =========
  # (1) Las filas futuras deben tener y = NA (jamás 0)
  is_future_cond <- is.na(data_cond$y)
  is_future_noP  <- is.na(data_noP$y)
  # Por las dudas, si en algún paso quedaron en 0, forzamos NA:
  data_cond$y[is_future_cond] <- NA_real_
  data_noP$y[is_future_noP]   <- NA_real_
  
  # (2) logE debe ser finito (si E==0 o faltante, lo arreglamos)
  fix_logE <- function(df) {
    if (!("logE" %in% names(df)) && "E" %in% names(df)) df$logE <- log(df$E)
    if ("logE" %in% names(df) && any(!is.finite(df$logE))) {
      if (!("E" %in% names(df))) stop("No tengo 'E' para recomputar logE")
      df$logE <- log(pmax(df$E, 1e-9))
    }
    df
  }
  data_cond <- fix_logE(data_cond)
  data_noP  <- fix_logE(data_noP)
  
  # (3) el offset externo debe ser finito
  fix_log_ext <- function(df) {
    if (!("log_mort_ext" %in% names(df))) stop("Falta 'log_mort_ext' en dataset de anclado")
    ok_hist <- is.finite(df$log_mort_ext) & !is.na(df$y)
    if (any(!is.finite(df$log_mort_ext))) {
      repl <- if (any(ok_hist)) mean(df$log_mort_ext[ok_hist]) else log(1e-12)
      df$log_mort_ext[!is.finite(df$log_mort_ext)] <- repl
    }
    df
  }
  data_cond <- fix_log_ext(data_cond)
  data_noP  <- fix_log_ext(data_noP)
  
  data_cond <- attach_edge_weights_hist(data_cond, stage = "mort_anchor")
  data_noP  <- attach_edge_weights_hist(data_noP,  stage = "mort_anchor")
  
  # --------- Downweight por años (solo F) ----------
  w_obs_cond <- rep(1, nrow(data_cond))
  w_obs_noP  <- rep(1, nrow(data_noP))
  
  if (length(mort_downweight_years) > 0) {
    w_obs_cond[data_cond$period %in% mort_downweight_years] <- mort_downweight_weight
    w_obs_noP [data_noP$period  %in% mort_downweight_years] <- mort_downweight_weight
    scale_cond <- as.numeric(scale_cond) * w_obs_cond
    scale_noP  <- as.numeric(scale_noP)  * w_obs_noP
  }
  
  w_fit_cond <- as.numeric(scale_cond) * dplyr::coalesce(data_cond$edge_weight, 1)
  w_fit_noP  <- as.numeric(scale_noP)  * dplyr::coalesce(data_noP$edge_weight, 1)
  
  # --------- FITS anclados ----------
  .chk("pre MORT-ANCHOR-cond")
  .check_no_mode_in_f(form_anchor)
  
  fit_anchor_cond <- inla_tag("MORT-ANCHOR-cond", formula = form_anchor, family = "poisson",
                              data = data_cond,
                              weights = w_fit_cond,
                              # OJO: sin 'E = ...' aquí, el offset va en la fórmula
                              control.fixed     = ctrl_fixed_anchor,
                              control.predictor = list(compute = TRUE),
                              control.compute   = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
                              control.inla      = list(diagonal = 1e-5)
  )
  
  .chk("pre MORT-ANCHOR-noP")
  .check_no_mode_in_f(form_anchor)
  fit_anchor_noP <- inla_tag("MORT-ANCHOR-noP", formula = form_anchor, family = "poisson",
                              data = data_noP,
                              weights = w_fit_noP,
                              control.fixed     = ctrl_fixed_anchor,
                              control.predictor = list(compute = TRUE),
                              control.compute   = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
                              control.inla      = list(diagonal = 1e-5)
  )
  
  # ---- Guardar beta del enlace Incidencia -> Mortalidad ----
  beta_cond <- extract_beta_I(fit_anchor_cond, MORT_I_LINK_MODE)
  beta_noP  <- if (exists("fit_anchor_noP")) extract_beta_I(fit_anchor_noP, MORT_I_LINK_MODE) else c(NA,NA,NA)
  
  # En params (para que viaje con el resultado)
  params$beta_I_mean      <- beta_cond[1]
  params$beta_I_lwr       <- beta_cond[2]
  params$beta_I_upr       <- beta_cond[3]
  params$beta_I_noP_mean  <- beta_noP[1]
  params$beta_I_noP_lwr   <- beta_noP[2]
  params$beta_I_noP_upr   <- beta_noP[3]
  
  # (Opcional) exponer al GlobalEnv para mirar rápido tras cada corrida
  if (isTRUE(EXPOSE_SELECTED_LAGS_TO_ENV)) {
    sx <- ifelse(sex_sel == "F", "F", "M")
    assign(paste0("BETA_I_SELECTED_", sx), params$beta_I_mean, envir = .GlobalEnv)
  }
  
  # ---------- 8) Predicciones anuales
  # Importante: en los modelos anclados de mortalidad el offset externo
  # (log_mort_ext + mort_anchor_tech_offset) debe reincorporarse explícitamente
  # al pasar del predictor lineal a medias en escala de conteos. En esta
  # arquitectura, usar exp(summary.linear.predictor) sin sumar el offset aplana
  # las diferencias entre escenarios aunque log_mort_ext sí cambie.
  .safe_off <- function(x) dplyr::coalesce(as.numeric(x), 0)
  
  .fixed_mean_multi_local <- function(fit, names_try, default = 0) {
    sf <- fit$summary.fixed
    if (is.null(sf) || is.null(rownames(sf))) return(default)
    hit <- intersect(as.character(names_try), rownames(sf))
    if (!length(hit)) return(default)
    as.numeric(sf[hit[1], "mean"])
  }
  
  .random_mean_map_local <- function(fit, name, id_vec) {
    sr <- fit$summary.random
    if (is.null(sr) || !name %in% names(sr)) return(rep(0, length(id_vec)))
    tbl <- sr[[name]]
    idx <- suppressWarnings(as.integer(id_vec))
    out <- rep(0, length(idx))
    ok <- !is.na(idx) & idx >= 1L & idx <= nrow(tbl)
    out[ok] <- as.numeric(tbl$mean[idx[ok]])
    out
  }
  
  .rebuild_anchor_future <- function(fit, data_df) {
    intercept <- .fixed_mean_multi_local(fit, c("(Intercept)"), 0)
    beta_t  <- .fixed_mean_multi_local(fit, c("mort_anchor_trend_t", "mort_trend_t"), 0)
    beta_t2 <- .fixed_mean_multi_local(fit, c("mort_anchor_trend_t2", "mort_trend_t2"), 0)
    
    age_eff <- .random_mean_map_local(fit, "age_id", data_df$age_id)
    per_eff <- .random_mean_map_local(fit, "period_id", data_df$period_id)
    coh_eff <- .random_mean_map_local(fit, "cohort_id_anchor", data_df$cohort_id_anchor)
    shock_eff <- if ("shock_id" %in% names(data_df)) {
      .random_mean_map_local(fit, "shock_id", data_df$shock_id)
    } else {
      rep(0, nrow(data_df))
    }
    
    mort_trend_t <- if ("mort_anchor_trend_t" %in% names(data_df)) {
      as.numeric(data_df$mort_anchor_trend_t)
    } else if ("mort_trend_t" %in% names(data_df)) {
      as.numeric(data_df$mort_trend_t)
    } else {
      rep(0, nrow(data_df))
    }
    
    mort_trend_t2 <- if ("mort_anchor_trend_t2" %in% names(data_df)) {
      as.numeric(data_df$mort_anchor_trend_t2)
    } else if ("mort_trend_t2" %in% names(data_df)) {
      as.numeric(data_df$mort_trend_t2)
    } else {
      rep(0, nrow(data_df))
    }
    
    other_term <- intercept +
      beta_t * mort_trend_t +
      beta_t2 * mort_trend_t2 +
      age_eff + per_eff + coh_eff + shock_eff
    
    input_term <- as.numeric(data_df$log_mort_ext) +
      .safe_off(data_df$mort_anchor_tech_offset)
    
    list(
      other_term = other_term,
      input_term = input_term,
      mu_hat = pmax(exp(other_term + input_term), 1e-12)
    )
  }
  
  lp_b <- fit_bapc$summary.linear.predictor
  pred_bapc <- df_all %>% dplyr::transmute(
    sex = sex_sel, period, hist_flag = !is.na(y),
    cohort_true = cohort_true,
    cohort_ref = cohort_clamped,
    mu_hat = exp(lp_b$mean),
    mu_lwr = pmax(0, exp(lp_b$`0.025quant`)),
    mu_upr = pmax(0, exp(lp_b$`0.975quant`))
  )
  pred_bapc <- apply_mort_cohort_fc_posthoc(
    pred_df = pred_bapc,
    summary_df = tryCatch(fit_bapc$summary.random$cohort_id_bapc, error = function(e) NULL),
    levels_vec = lev_coh,
    target_levels = df_all$cohort_true,
    ref_levels = df_all$cohort_clamped
  )
  
  lp_ac <- fit_anchor_cond$summary.linear.predictor
  fv_ac <- tryCatch(fit_anchor_cond$summary.fitted.values, error = function(e) NULL)
  use_fv_ac <- !is.null(fv_ac) && is.data.frame(fv_ac) && nrow(fv_ac) == nrow(data_cond)
  fv_mean_ac <- if (use_fv_ac) as.numeric(fv_ac$mean) else rep(NA_real_, nrow(data_cond))
  fv_lwr_ac  <- if (use_fv_ac) as.numeric(fv_ac$`0.025quant`) else rep(NA_real_, nrow(data_cond))
  fv_upr_ac  <- if (use_fv_ac) as.numeric(fv_ac$`0.975quant`) else rep(NA_real_, nrow(data_cond))
  rebuild_ac <- .rebuild_anchor_future(fit_anchor_cond, data_cond)
  
  pred_anchor_cond <- data_cond %>% dplyr::transmute(
    sex = sex_sel, period, hist_flag = !is.na(y),
    cohort_true = cohort_true,
    cohort_ref = cohort_ref_anchor,
    log_mort_ext = as.numeric(log_mort_ext),
    mort_anchor_tech_offset = .safe_off(mort_anchor_tech_offset),
    offset_total = log_mort_ext + mort_anchor_tech_offset,
    eta_lp_hat = as.numeric(lp_ac$mean),
    eta_lp_lwr = as.numeric(lp_ac$`0.025quant`),
    eta_lp_upr = as.numeric(lp_ac$`0.975quant`),
    mu_hat_legacy = exp(eta_lp_hat),
    mu_hat = dplyr::if_else(
      hist_flag & !is.null(fv_ac) & is.data.frame(fv_ac) & nrow(fv_ac) == n(),
      as.numeric(fv_ac$mean),
      rebuild_ac$mu_hat
    ),
    mu_lwr = dplyr::if_else(
      hist_flag & !is.null(fv_ac) & is.data.frame(fv_ac) & nrow(fv_ac) == n(),
      pmax(0, as.numeric(fv_ac$`0.025quant`)),
      pmax(rebuild_ac$mu_hat * exp(eta_lp_lwr - eta_lp_hat), 1e-12)
    ),
    mu_upr = dplyr::if_else(
      hist_flag & !is.null(fv_ac) & is.data.frame(fv_ac) & nrow(fv_ac) == n(),
      pmax(0, as.numeric(fv_ac$`0.975quant`)),
      pmax(rebuild_ac$mu_hat * exp(eta_lp_upr - eta_lp_hat), 1e-12)
    ),
    eta_resid_hat = log(pmax(mu_hat, 1e-12)) - offset_total,
    eta_resid_lwr = log(pmax(mu_lwr, 1e-12)) - offset_total,
    eta_resid_upr = log(pmax(mu_upr, 1e-12)) - offset_total,
    mort_ext_deaths = dplyr::coalesce(as.numeric(mort_ext_deaths), 0)  )
  
  pred_anchor_cond <- apply_mort_cohort_fc_posthoc(
    pred_df = pred_anchor_cond,
    summary_df = tryCatch(fit_anchor_cond$summary.random$cohort_id_anchor, error = function(e) NULL),
    levels_vec = lev_coh,
    target_levels = data_cond$cohort_true,
    ref_levels = data_cond$cohort_ref_anchor
  )
  
  lp_an <- fit_anchor_noP$summary.linear.predictor
  fv_an <- tryCatch(fit_anchor_noP$summary.fitted.values, error = function(e) NULL)
  use_fv_an <- !is.null(fv_an) && is.data.frame(fv_an) && nrow(fv_an) == nrow(data_noP)
  fv_mean_an <- if (use_fv_an) as.numeric(fv_an$mean) else rep(NA_real_, nrow(data_noP))
  fv_lwr_an  <- if (use_fv_an) as.numeric(fv_an$`0.025quant`) else rep(NA_real_, nrow(data_noP))
  fv_upr_an  <- if (use_fv_an) as.numeric(fv_an$`0.975quant`) else rep(NA_real_, nrow(data_noP))
  rebuild_an <- .rebuild_anchor_future(fit_anchor_noP, data_noP)
  
  pred_anchor_noP <- data_noP %>% dplyr::transmute(
    period, hist_flag = !is.na(y),
    cohort_true = cohort_true,
    cohort_ref = cohort_ref_anchor,
    log_mort_ext = as.numeric(log_mort_ext),
    mort_anchor_tech_offset = .safe_off(mort_anchor_tech_offset),
    offset_total = log_mort_ext + mort_anchor_tech_offset,
    eta_lp_hat = as.numeric(lp_an$mean),
    eta_lp_lwr = as.numeric(lp_an$`0.025quant`),
    eta_lp_upr = as.numeric(lp_an$`0.975quant`),
    sex = sex_sel,
    mu_hat_legacy = exp(eta_lp_hat),
    mu_hat = dplyr::if_else(
      hist_flag & !is.null(fv_an) & is.data.frame(fv_an) & nrow(fv_an) == n(),
      as.numeric(fv_an$mean),
      rebuild_an$mu_hat
    ),
    mu_lwr = dplyr::if_else(
      hist_flag & !is.null(fv_an) & is.data.frame(fv_an) & nrow(fv_an) == n(),
      pmax(0, as.numeric(fv_an$`0.025quant`)),
      pmax(rebuild_an$mu_hat * exp(eta_lp_lwr - eta_lp_hat), 1e-12)
    ),
    mu_upr = dplyr::if_else(
      hist_flag & !is.null(fv_an) & is.data.frame(fv_an) & nrow(fv_an) == n(),
      pmax(0, as.numeric(fv_an$`0.975quant`)),
      pmax(rebuild_an$mu_hat * exp(eta_lp_upr - eta_lp_hat), 1e-12)
    ),
    eta_resid_hat = log(pmax(mu_hat, 1e-12)) - offset_total,
    eta_resid_lwr = log(pmax(mu_lwr, 1e-12)) - offset_total,
    eta_resid_upr = log(pmax(mu_upr, 1e-12)) - offset_total,
    mort_ext_deaths = dplyr::coalesce(as.numeric(mort_ext_deaths), 0)
  )
  pred_anchor_noP <- apply_mort_cohort_fc_posthoc(
    pred_df = pred_anchor_noP,
    summary_df = tryCatch(fit_anchor_noP$summary.random$cohort_id_anchor, error = function(e) NULL),
    levels_vec = lev_coh,
    target_levels = data_noP$cohort_true,
    ref_levels = data_noP$cohort_ref_anchor
  )
  
  annual_external_cond <- summarise_annual_external(pred_anchor_cond)
  
  annual_external_noP <- summarise_annual_external(pred_anchor_noP)
  
  annual_anchor_components_cond <- summarise_annual_anchor_components(pred_anchor_cond)
  
  annual_anchor_components_noP <- summarise_annual_anchor_components(pred_anchor_noP)
  
  annual_bapc <- summarise_annual_prediction(pred_bapc)
  
  annual_anchor <- summarise_annual_prediction(pred_anchor_cond)
  
  annual_anchor_noP <- summarise_annual_prediction(pred_anchor_noP)
  
  annual_bapc_raw       <- annual_bapc
  annual_anchor_raw     <- annual_anchor
  annual_anchor_noP_raw <- annual_anchor_noP
  
  obs_annual <- mortH %>% dplyr::group_by(period) %>%
    dplyr::summarise(obs = sum(deaths, na.rm=TRUE), .groups="drop") %>%
    dplyr::mutate(sex = sex_sel, .before = 1)
  
  obs_2022 <- obs_annual$obs[obs_annual$period == last_hist_year]
  hat_b_2022 <- annual_bapc$deaths_hat[annual_bapc$period == last_hist_year]
  k_b <- ifelse(is.finite(obs_2022) & obs_2022 > 0 & is.finite(hat_b_2022) & hat_b_2022 > 0, obs_2022 / hat_b_2022, 1)

  # Empalme anual: recurso técnico, no epidemiológico. Debe poder apagarse.
  if (isTRUE(MORT_ANNUAL_BRIDGE)) {
    annual_bapc <- annual_bapc %>%
      dplyr::mutate(
        deaths_hat = ifelse(period > last_hist_year, deaths_hat * k_b, deaths_hat),
        deaths_lwr = ifelse(period > last_hist_year, deaths_lwr * k_b, deaths_lwr),
        deaths_upr = ifelse(period > last_hist_year, deaths_upr * k_b, deaths_upr)
      )

    hat_an_2022 <- annual_anchor$deaths_hat[annual_anchor$period == last_hist_year]
    k_an <- ifelse(is.finite(obs_2022) & obs_2022 > 0 & is.finite(hat_an_2022) & hat_an_2022 > 0, obs_2022 / hat_an_2022, 1)
    annual_anchor <- annual_anchor %>%
      dplyr::mutate(
        deaths_hat = ifelse(period > last_hist_year, deaths_hat * k_an, deaths_hat),
        deaths_lwr = ifelse(period > last_hist_year, deaths_lwr * k_an, deaths_lwr),
        deaths_upr = ifelse(period > last_hist_year, deaths_upr * k_an, deaths_upr)
      )

    hat_noP_2022 <- annual_anchor_noP$deaths_hat[annual_anchor_noP$period == last_hist_year]
    k_noP <- ifelse(is.finite(obs_2022) & obs_2022 > 0 & is.finite(hat_noP_2022) & hat_noP_2022 > 0, obs_2022 / hat_noP_2022, 1)
    annual_anchor_noP <- annual_anchor_noP %>%
      dplyr::mutate(
        deaths_hat = ifelse(period > last_hist_year, deaths_hat * k_noP, deaths_hat),
        deaths_lwr = ifelse(period > last_hist_year, deaths_lwr * k_noP, deaths_lwr),
        deaths_upr = ifelse(period > last_hist_year, deaths_upr * k_noP, deaths_upr)
      )
  } else {
    k_b <- 1
    k_an <- 1
    k_noP <- 1
  }

  # ==== DIAGNÓSTICO PREV→INC (por sexo/causa) ====
  diag_prev <- tryCatch({
    .scalar_num <- function(x) {
      x <- suppressWarnings(as.numeric(x))
      x <- x[is.finite(x)]
      if (!length(x)) return(NA_real_)
      unname(x[[1]])
    }
    .median_num <- function(x) {
      x <- suppressWarnings(as.numeric(x))
      x <- x[is.finite(x)]
      if (!length(x)) return(NA_real_)
      stats::median(x, na.rm = TRUE)
    }
    .scalar_chr <- function(x, default = NA_character_) {
      x <- as.character(x)
      x <- x[!is.na(x) & nzchar(x)]
      if (!length(x)) return(default)
      x[[1]]
    }

    beta_rule_cur <- .scalar_chr(tryCatch(inc_fit$beta_P_rule, error = function(e) NA_character_),
                                 default = .scalar_chr(beta_mode, default = NA_character_))
    fixed_rr_mode <- identical(beta_rule_cur, "fixed_rr_offset") || identical(as.character(beta_mode)[1], "fixed_rr_offset")

    z_hist_src <- tryCatch(inc_fit$z_hist, error = function(e) NULL)
    if (is.null(z_hist_src) || !nrow(z_hist_src)) {
      z_hist_src <- inc_inla %>% dplyr::select(dplyr::any_of(c("age", "period", "q_eff", "z_prev")))
    }
    z_fut_src <- tryCatch(inc_fit$z_future, error = function(e) NULL)
    if (is.null(z_fut_src) || !nrow(z_fut_src)) {
      z_fut_src <- fut_grid %>% dplyr::select(dplyr::any_of(c("age", "period", "q_eff", "z_prev")))
    }

    z_hist <- z_hist_src |> dplyr::filter(period == last_hist_year) |> dplyr::pull(z_prev)
    z_proj_end <- z_fut_src |> dplyr::filter(period == PROJ_TO) |> dplyr::pull(z_prev)
    q_hist <- if ("q_eff" %in% names(z_hist_src)) z_hist_src |> dplyr::filter(period == last_hist_year) |> dplyr::pull(q_eff) else numeric()
    q_proj_end <- if ("q_eff" %in% names(z_fut_src)) z_fut_src |> dplyr::filter(period == PROJ_TO) |> dplyr::pull(q_eff) else numeric()

    z_base_med <- .median_num(z_hist)
    z_proj_end_med <- .median_num(z_proj_end)
    q_base_med <- .median_num(q_hist)
    q_proj_end_med <- .median_num(q_proj_end)
    rr_proj_end <- if (fixed_rr_mode && is.finite(z_base_med) && is.finite(z_proj_end_med)) {
      exp(z_proj_end_med - z_base_med)
    } else if (is.finite(betaP_eff) && is.finite(z_base_med) && is.finite(z_proj_end_med)) {
      exp(betaP_eff * (z_proj_end_med - z_base_med))
    } else {
      NA_real_
    }

    tibble::tibble(
      cause           = .scalar_chr(unique(mort_hist_tbl$cause)),
      sex             = .scalar_chr(sex_sel),
      last_hist_year  = as.integer(last_hist_year),
      beta_mode       = .scalar_chr(beta_mode),
      rr_inc          = .scalar_num(rr_inc_cur),
      inc_trend_on    = .scalar_chr(INC_TREND_ON),
      inc_trend_degree= as.integer(if (isTRUE(inc_include_trend)) INC_TREND_DEGREE else 0L),
      inc_trend_scenario = .scalar_chr(inc_trend_scenario),
      inc_trend_coef_t  = { inc_fix <- tryCatch(inc_fit$fit_inc$summary.fixed, error = function(e) NULL); if (!is.null(inc_fix) && "inc_trend_t" %in% rownames(inc_fix)) .scalar_num(inc_fix["inc_trend_t", "mean"]) else NA_real_ },
      inc_trend_coef_t2 = { inc_fix <- tryCatch(inc_fit$fit_inc$summary.fixed, error = function(e) NULL); if (!is.null(inc_fix) && "inc_trend_t2" %in% rownames(inc_fix)) .scalar_num(inc_fix["inc_trend_t2", "mean"]) else NA_real_ },
      beta_P_rule     = beta_rule_cur,
      prev_sign       = suppressWarnings(as.integer(prev_sign_)),
      s_histP         = .scalar_num(tryCatch(inc_fit$s_histP, error = function(e) NA_real_)),
      q_eff_base      = q_base_med,
      q_eff_proj_end      = q_proj_end_med,
      offset_inc_base = z_base_med,
      offset_inc_proj_end = z_proj_end_med,
      z_base          = z_base_med,
      z_proj_end          = z_proj_end_med,
      RR_proj_end_vs_base = rr_proj_end,
      prev_error      = NA_character_
    )
  }, error = function(e) {
    beta_rule_cur <- as.character(tryCatch(inc_fit$beta_P_rule, error = function(err) NA_character_))[1] %||% as.character(beta_mode)[1]
    fixed_rr_mode <- identical(beta_rule_cur, "fixed_rr_offset") || identical(as.character(beta_mode)[1], "fixed_rr_offset")
    tibble::tibble(
      cause           = as.character(unique(mort_hist_tbl$cause))[1] %||% NA_character_,
      sex             = as.character(sex_sel)[1] %||% NA_character_,
      last_hist_year  = as.integer(last_hist_year),
      beta_mode       = as.character(beta_mode)[1] %||% NA_character_,
      rr_inc          = suppressWarnings(as.numeric(rr_inc_cur))[1] %||% NA_real_,
      inc_trend_on    = as.character(INC_TREND_ON)[1] %||% NA_character_,
      inc_trend_degree= as.integer(if (isTRUE(inc_include_trend)) INC_TREND_DEGREE else 0L),
      inc_trend_scenario = as.character(inc_trend_scenario)[1] %||% NA_character_,
      inc_trend_coef_t  = NA_real_,
      inc_trend_coef_t2 = NA_real_,
      beta_P_rule     = beta_rule_cur %||% NA_character_,
      prev_sign       = suppressWarnings(as.integer(prev_sign_))[1] %||% NA_integer_,
      s_histP         = NA_real_,
      q_eff_base      = NA_real_,
      q_eff_proj_end      = NA_real_,
      offset_inc_base = NA_real_,
      offset_inc_proj_end = NA_real_,
      z_base          = NA_real_,
      z_proj_end          = NA_real_,
      RR_proj_end_vs_base = NA_real_,
      prev_error      = conditionMessage(e)
    )
  })

  # ==========================================================
  # ΔLCPO (alt - benchmark) con SE + backtest temporal (mortalidad)
  # ==========================================================
  d_inc_ip  <- delta_lcpo_se(inc_fit$fit_inc,      inc_fit_bapc$fit_inc)
  d_m_ip    <- delta_lcpo_se(fit_anchor_cond,      fit_bapc)
  d_m_i     <- delta_lcpo_se(fit_anchor_noP,       fit_bapc)
  
  bt_m_bench <- bt_m_ip <- bt_m_i <- NULL
  if (isTRUE(BT_ENABLE) && isTRUE(BT_HOLDOUT_YEARS > 0)) {
    train_end_m <- period_max_m - as.integer(BT_HOLDOUT_YEARS)
    
    bt_m_bench <- backtest_inla_poisson(
      formula = form_bapc, data = df_bapc,
      train_end = train_end_m,
      tag = sprintf("BT-MORT-BENCH-%s", sex_sel),
            control.fixed  = ctrl_fixed_bapc,
      control.inla   = list(diagonal = 1e-5)
    )
    
    bt_m_ip <- backtest_inla_poisson(
      formula = form_anchor, data = data_cond,
      train_end = train_end_m,
      tag = sprintf("BT-MORT-IP-%s", sex_sel),
            control.fixed  = ctrl_fixed_anchor,
      control.inla   = list(diagonal = 1e-5)
    )
    
    bt_m_i <- backtest_inla_poisson(
      formula = form_anchor, data = data_noP,
      train_end = train_end_m,
      tag = sprintf("BT-MORT-I-%s", sex_sel),
            control.fixed  = ctrl_fixed_anchor,
      control.inla   = list(diagonal = 1e-5)
    )
  }
  
  # ===== 8.5) Scores de ajuste (para tabla LaTeX) + extras =====
  fit_scores_base <- dplyr::bind_rows(
    score_row("Prevalence BAPC", fit_prev),
    score_row("Incidence benchmark (APC)", inc_fit_bapc$fit_inc),
    score_row("Incidence prevalence-informed", inc_fit$fit_inc),
    score_row("Mortality benchmark (APC)", fit_bapc),
    score_row("Mortality anchored on I|P", fit_anchor_cond),
    score_row("Mortality anchored on I only", fit_anchor_noP)
  ) %>%
    dplyr::mutate(sex = sex_sel) %>%
    dplyr::select(sex, model, WAIC, DIC, LCPO)
  
  extras <- dplyr::bind_rows(
    data.frame(model = "Incidence benchmark (APC)",
               dLCPO = NA_real_, se_dLCPO = NA_real_,
               BT_LPD = inc_fit_bapc$bt$BT_LPD %||% NA_real_,
               BT_RMSE = inc_fit_bapc$bt$BT_RMSE %||% NA_real_,
               stringsAsFactors = FALSE),
    data.frame(model = "Incidence prevalence-informed",
               dLCPO = d_inc_ip$delta, se_dLCPO = d_inc_ip$se,
               BT_LPD = inc_fit$bt$BT_LPD %||% NA_real_,
               BT_RMSE = inc_fit$bt$BT_RMSE %||% NA_real_,
               stringsAsFactors = FALSE),
    data.frame(model = "Mortality benchmark (APC)",
               dLCPO = NA_real_, se_dLCPO = NA_real_,
               BT_LPD = bt_m_bench$BT_LPD %||% NA_real_,
               BT_RMSE = bt_m_bench$BT_RMSE %||% NA_real_,
               stringsAsFactors = FALSE),
    data.frame(model = "Mortality anchored on I|P",
               dLCPO = d_m_ip$delta, se_dLCPO = d_m_ip$se,
               BT_LPD = bt_m_ip$BT_LPD %||% NA_real_,
               BT_RMSE = bt_m_ip$BT_RMSE %||% NA_real_,
               stringsAsFactors = FALSE),
    data.frame(model = "Mortality anchored on I only",
               dLCPO = d_m_i$delta, se_dLCPO = d_m_i$se,
               BT_LPD = bt_m_i$BT_LPD %||% NA_real_,
               BT_RMSE = bt_m_i$BT_RMSE %||% NA_real_,
               stringsAsFactors = FALSE)
  )
  
  fit_scores <- fit_scores_base %>%
    dplyr::left_join(extras, by = "model") %>%
    dplyr::select(sex, model, WAIC, DIC, LCPO, dLCPO, se_dLCPO, BT_LPD, BT_RMSE)
  
  
  
  
  proj_horizon_info <- tryCatch(
    projection_horizon_from_border_diag(inc_fit$border_diag_future, exposure_col = "E"),
    error = function(e) list(year_diag = tibble::tibble(), frontier = tibble::tibble())
  )

  ## ===== 9) Salida =====
  out <- list(
    sex = sex_sel,
    params = params, 
    annual_bapc   = annual_bapc,
    annual_bapc_raw = annual_bapc_raw,
    annual_external_cond = annual_external_cond,
    annual_external_noP = annual_external_noP,
    annual_anchor_components_cond = annual_anchor_components_cond,
    annual_anchor_components_noP = annual_anchor_components_noP,
    annual_anchor = annual_anchor,           # I|P (principal)
    annual_anchor_raw = annual_anchor_raw,
    annual_anchor_noP = annual_anchor_noP,   # I sin P (contrafactual)
    annual_anchor_noP_raw = annual_anchor_noP_raw,
    mort_anchor_pred_detail = pred_anchor_cond,
    mort_anchor_noP_pred_detail = pred_anchor_noP,
    mort_anchor_data_cond = data_cond,
    mort_anchor_data_noP = data_noP,
    obs_annual    = obs_annual,
    
    inc_fit = list(
      fit           = inc_fit$fit_inc,
      rates_all     = inc_fit$rates_all,
      rates_all_full= inc_fit$rates_all_full,
      lev_inc       = inc_fit$lev_inc,
      last_year_inc = inc_fit$last_year_inc,
      border_diag_future = tryCatch(inc_fit$border_diag_future, error = function(e) tibble::tibble()),
      border_diag_summary = tryCatch(inc_fit$border_diag_summary, error = function(e) tibble::tibble()),
      projection_horizon_year = proj_horizon_info$year_diag,
      projection_horizon_frontier = proj_horizon_info$frontier,
      beta_P        = if (identical(beta_mode, "fixed_rr_offset")) NA_real_ else betaP_raw,
      beta_P_pos    = if (identical(beta_mode, "fixed_rr_offset")) NA_real_ else betaP_eff,
      beta_P_eff    = if (identical(beta_mode, "fixed_rr_offset")) NA_real_ else betaP_eff,
      beta_P_zeroed = if (identical(beta_mode, "fixed_rr_offset")) NA else betaP_zeroed,
      beta_P_rule   = tryCatch(inc_fit$beta_P_rule, error = function(e) if (identical(beta_mode, "fixed_rr_offset")) "fixed_rr_offset" else BETA_P_POSTFIT_RULE)
    ),
    inc_fit_bapc = list(
      fit           = inc_fit_bapc$fit_inc,
      rates_all     = inc_fit_bapc$rates_all,
      rates_all_full = inc_fit_bapc$rates_all_full,
      lev_inc       = inc_fit_bapc$lev_inc,
      last_year_inc = inc_fit_bapc$last_year_inc
    ),
    inc_annual_bapc = inc_annual_bapc,
    inc_annual_cond = inc_annual_cond,
    inc_annual_noP = inc_annual_noP,
    inc_obs_annual  = inc_obs_annual,
    
    fit_prev   = fit_prev,
    fit_bapc   = fit_bapc,
    fit_anchor = fit_anchor_cond,
    fit_anchor_noP = fit_anchor_noP,
    fit_scores = fit_scores,
    
    lev_prev = list(age = sort(unique(prev_inla$age)),
                    period = sort(unique(prev_inla$period)),
                    cohort = sort(unique(prev_inla$cohort))),
    lev_mort = list(age = lev_age, period = lev_per, cohort = lev_coh),
    
    diag = list(
      last_hist_year = last_hist_year,
      L_I = as.integer(L_I_eff), Da_I = Da_I,
      L_I_input = L_I,
      mort_link_mode = MORT_I_LINK_MODE,
      mort_kernel = kernel_tbl,
      mort_kernel_summary_cond = kernel_summary_cond,
      mort_kernel_summary_noP = kernel_summary_noP,
      prev_sign = prev_sign_,
      mort_trend_scenario = mort_trend_scenario,
      mort_bapc_trend_scenario = mort_bapc_trend_scenario,
      mort_bapc_future_mode = "autonomous_apc",
      delta_tech = delta_tech,
      df_all = df_all,
      z_prev_hist = as.data.frame(inc_fit$z_hist %||% tibble::tibble()),
      z_prev_future = as.data.frame(inc_fit$z_future %||% tibble::tibble()),
      s_histP = tryCatch(inc_fit$s_histP, error = function(e) NA_real_),
      prev = diag_prev,
      bridge_factors = tibble::tibble(
        sex = sex_sel,
        k_b = as.numeric(k_b),
        k_an = as.numeric(k_an),
        k_inc_noP_cf = as.numeric(k_noP)
      ),
      projection_horizon_year = proj_horizon_info$year_diag,
      projection_horizon_frontier = proj_horizon_info$frontier,
      max_projection_year_endogenous = projection_max_year_from_frontier(proj_horizon_info$frontier, policy = "endogenous_max"),
      method_policy_table = method_policy_tbl
    )
  )
  # out <- sanitize_pipeline_output(out)
  return(out)
}


# =============================================================
