# 2) Helpers (utilidades generales)
# =============================================================

# A partir del último año observado calcula horizonte y años futuros
proj_horizon <- function(last_hist_year, proj_to = PROJ_TO) {
  if (proj_to < last_hist_year + 1)
    stop("proj_to debe ser >= último observado + 1")
  yrs <- seq.int(last_hist_year + 1, proj_to)
  list(n_ahead = length(yrs), years = yrs)
}

# Filtro práctico por año máximo (para tablas o plots)
clip_to_proj <- function(df, year_var = "period", proj_to = PROJ_TO) {
  df[df[[year_var]] <= proj_to, , drop = FALSE]
}


.bapc_verbose <- function(..., verbose = BAPC_VERBOSE) {
  if (isTRUE(verbose)) message(...)
  invisible(NULL)
}

.bapc_info <- function(..., verbose = BAPC_VERBOSE) {
  if (!isTRUE(verbose)) message(...)
  invisible(NULL)
}



# Exposure-weighted annual Lexis/border diagnostics used to define the
# endogenous reporting horizon for projected outputs.
make_projection_horizon_year_diag <- function(border_df, exposure_col = NULL) {
  if (!is.data.frame(border_df) || !nrow(border_df)) {
    return(tibble::tibble())
  }

  exposure_col <- exposure_col %||% intersect(c("exposure", "E"), names(border_df))[1]
  if (is.na(exposure_col) || !nzchar(exposure_col) || !(exposure_col %in% names(border_df))) {
    border_df$.w_h <- 1
  } else {
    border_df$.w_h <- suppressWarnings(as.numeric(border_df[[exposure_col]]))
    border_df$.w_h[!is.finite(border_df$.w_h) | border_df$.w_h <= 0] <- 1
  }

  .wmean_h <- function(x, w) {
    ok <- is.finite(x) & is.finite(w) & w > 0
    if (!any(ok)) return(NA_real_)
    stats::weighted.mean(x[ok], w[ok])
  }

  border_df %>%
    dplyr::mutate(
      sex = as.character(sex),
      period = suppressWarnings(as.integer(period)),
      horizon = suppressWarnings(as.numeric(horizon))
    ) %>%
    dplyr::filter(is.finite(period), is.finite(horizon)) %>%
    dplyr::group_by(sex, period, horizon) %>%
    dplyr::summarise(
      n_cells = dplyr::n(),
      share_cohort_edge = .wmean_h(as.numeric(cohort_is_edge %in% TRUE), .w_h),
      share_cohort_low = .wmean_h(as.numeric(cohort_clamp_low %in% TRUE), .w_h),
      share_cohort_high = .wmean_h(as.numeric(cohort_clamp_high %in% TRUE), .w_h),
      mean_abs_cohort_shift = .wmean_h(abs(suppressWarnings(as.numeric(cohort_shift))), .w_h),
      mean_support_frac = .wmean_h(suppressWarnings(as.numeric(support_frac)), .w_h),
      .groups = "drop"
    ) %>%
    dplyr::arrange(sex, horizon, period)
}

build_projection_horizon_frontier <- function(year_diag,
                                              support_floor_credible = HORIZON_SUPPORT_FLOOR_CREDIBLE,
                                              support_floor_caution  = HORIZON_SUPPORT_FLOOR_CAUTION,
                                              support_floor_max      = HORIZON_SUPPORT_FLOOR_MAX,
                                              edge_share_credible    = HORIZON_EDGE_SHARE_CREDIBLE,
                                              edge_share_caution     = HORIZON_EDGE_SHARE_CAUTION,
                                              edge_share_max         = HORIZON_EDGE_SHARE_MAX) {
  if (!is.data.frame(year_diag) || !nrow(year_diag)) {
    return(list(year_diag = tibble::tibble(), frontier = tibble::tibble()))
  }

  d <- year_diag %>%
    dplyr::mutate(
      mild_flag = (is.finite(mean_support_frac) & mean_support_frac < support_floor_credible) |
        (is.finite(share_cohort_edge) & share_cohort_edge >= edge_share_credible),
      severe_flag = (is.finite(mean_support_frac) & mean_support_frac < support_floor_caution) |
        (is.finite(share_cohort_edge) & share_cohort_edge >= edge_share_caution),
      terminal_flag = (is.finite(mean_support_frac) & mean_support_frac < support_floor_max) |
        (is.finite(share_cohort_edge) & share_cohort_edge >= edge_share_max)
    )

  first_hit <- function(x_horizon, cond) {
    cond <- cond %in% TRUE
    idx <- which(cond)
    if (!length(idx)) return(NA_real_)
    suppressWarnings(as.numeric(x_horizon[min(idx)]))
  }

  frontier <- d %>%
    dplyr::group_by(sex) %>%
    dplyr::summarise(
      last_hist_year = suppressWarnings(as.integer(stats::median(period - horizon, na.rm = TRUE))),
      max_horizon_available = suppressWarnings(as.numeric(max(horizon, na.rm = TRUE))),
      first_horizon_caution = first_hit(horizon, mild_flag),
      first_horizon_risky = first_hit(horizon, severe_flag),
      first_horizon_beyond_max = first_hit(horizon, terminal_flag),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      end_horizon_credible = dplyr::if_else(is.finite(first_horizon_caution), pmax(0, first_horizon_caution - 1), max_horizon_available),
      end_horizon_caution  = dplyr::if_else(is.finite(first_horizon_risky), pmax(end_horizon_credible, first_horizon_risky - 1), max_horizon_available),
      end_horizon_risky    = dplyr::if_else(is.finite(first_horizon_beyond_max), pmax(end_horizon_caution, first_horizon_beyond_max - 1), max_horizon_available),
      end_year_credible = last_hist_year + end_horizon_credible,
      end_year_caution  = last_hist_year + end_horizon_caution,
      end_year_risky    = last_hist_year + end_horizon_risky,
      max_projection_year_endogenous = end_year_risky,
      support_floor_credible = support_floor_credible,
      support_floor_caution  = support_floor_caution,
      support_floor_max      = support_floor_max,
      edge_share_credible    = edge_share_credible,
      edge_share_caution     = edge_share_caution,
      edge_share_max         = edge_share_max
    )

  d <- d %>%
    dplyr::left_join(frontier %>% dplyr::select(sex, end_horizon_credible, end_horizon_caution, end_horizon_risky), by = "sex") %>%
    dplyr::mutate(
      projection_zone = dplyr::case_when(
        horizon <= end_horizon_credible ~ "credible",
        horizon <= end_horizon_caution  ~ "caution",
        horizon <= end_horizon_risky    ~ "risky",
        TRUE ~ "beyond_max"
      ),
      projection_zone = factor(projection_zone, levels = c("credible", "caution", "risky", "beyond_max"))
    )

  list(year_diag = d, frontier = frontier)
}

projection_horizon_from_border_diag <- function(border_df, exposure_col = NULL) {
  yd <- make_projection_horizon_year_diag(border_df, exposure_col = exposure_col)
  build_projection_horizon_frontier(yd)
}

projection_max_year_from_frontier <- function(frontier_tbl, policy = c("endogenous_max", "risky", "caution", "credible")) {
  policy <- match.arg(policy)
  if (!is.data.frame(frontier_tbl) || !nrow(frontier_tbl)) return(NA_integer_)
  col <- switch(policy,
    endogenous_max = "max_projection_year_endogenous",
    risky = "end_year_risky",
    caution = "end_year_caution",
    credible = "end_year_credible"
  )
  out <- suppressWarnings(as.integer(frontier_tbl[[col]][1]))
  if (!is.finite(out)) NA_integer_ else out
}

projection_max_year_from_res_sex <- function(res_sex, policy = c("endogenous_max", "risky", "caution", "credible"), default = NA_integer_) {
  policy <- match.arg(policy)
  frontier_tbl <- tryCatch(res_sex$diag$projection_horizon_frontier, error = function(e) NULL)
  out <- projection_max_year_from_frontier(frontier_tbl, policy = policy)
  if (!is.finite(out)) suppressWarnings(as.integer(default)[1]) else out
}

projection_common_max_year_from_res_both <- function(res_both, policy = c("endogenous_max", "risky", "caution", "credible"), default = NA_integer_) {
  policy <- match.arg(policy)
  vals <- c(
    tryCatch(projection_max_year_from_res_sex(res_both$resM, policy = policy), error = function(e) NA_integer_),
    tryCatch(projection_max_year_from_res_sex(res_both$resF, policy = policy), error = function(e) NA_integer_)
  )
  vals <- vals[is.finite(vals)]
  if (!length(vals)) return(suppressWarnings(as.integer(default)[1]))
  suppressWarnings(as.integer(min(vals, na.rm = TRUE)))
}

clip_to_year <- function(df, max_year, year_var = "period") {
  if (!is.data.frame(df) || !nrow(df) || !nzchar(year_var) || !(year_var %in% names(df))) return(df)
  max_year <- suppressWarnings(as.integer(max_year))[1]
  if (!is.finite(max_year)) return(df)
  df[df[[year_var]] <= max_year, , drop = FALSE]
}

cs01 <- function(x) { m <- mean(x, na.rm = TRUE); s <- stats::sd(x, na.rm = TRUE); if (is.na(s)||s==0) return(rep(0,length(x))); (x-m)/s }
assert_no_na <- function(x, msg) { if (anyNA(x)) stop(msg) }

check_apc_grid <- function(df) {
  dfa <- df %>% mutate(check = period - age)
  if (any(dfa$check != dfa$cohort, na.rm = TRUE)) {
    bad <- dfa %>% filter(check != cohort) %>% head(10)
    stop("Cohort != period - age en algunas filas.\n", paste(capture.output(print(bad)), collapse = "\n"))
  }
  if (any(df$exposure < 0, na.rm = TRUE)) stop("exposure negativo detectado")
  dup <- df %>% count(age, period, sex, cause) %>% filter(n > 1)
  if (nrow(dup) > 0) stop("Celdas duplicadas (age,period,sex,cause).")
  invisible(TRUE)
}

make_indices <- function(df) {
  df %>% mutate(
    age_id    = as.integer(factor(age,    levels = sort(unique(age)))),
    period_id = as.integer(factor(period, levels = sort(unique(period)))),
    cohort_id = as.integer(factor(cohort, levels = sort(unique(cohort))))
  )
}

kish_neff <- function(w) { s1 <- sum(w); s2 <- sum(w^2); if (s2 <= 0) return(0); (s1^2) / s2 }

pc_hyper <- function(u, a) list(prec = list(prior = "pc.prec", param = c(u, a)))


.beta_postfit_transform <- function(beta_raw, rule = BETA_P_POSTFIT_RULE) {
  rule <- match.arg(rule, choices = c("floor0", "softplus", "identity"))
  if (length(beta_raw) == 0 || !is.finite(beta_raw)) {
    return(list(raw = if(length(beta_raw)==0) NA_real_ else beta_raw, eff = NA_real_, zeroed = NA, rule = rule))
  }
  eff <- switch(rule,
    floor0   = max(beta_raw, 0),
    softplus = log1p(exp(beta_raw)),
    identity = beta_raw
  )
  zeroed <- isTRUE(rule == "floor0") && is.finite(beta_raw) && beta_raw < 0 && isTRUE(all.equal(eff, 0))
  list(raw = beta_raw, eff = as.numeric(eff), zeroed = zeroed, rule = rule)
}

.beta_eff <- function(beta_raw, rule = BETA_P_POSTFIT_RULE) {
  .beta_postfit_transform(beta_raw = beta_raw, rule = rule)$eff
}


# ---------- Forecast γ^P y otros coeficientes APC
.forecast_apc_tail <- function(levels_hist, values_hist, future_levels,
                               method = c("freeze", "arima", "trend", "damped_trend"),
                               arima_args = list(stepwise = FALSE, approximation = FALSE, seasonal = FALSE),
                               trend_type = c("trend", "level"),
                               window = COHORT_FC_WINDOW,
                               damping = COHORT_FC_DAMPING) {
  method <- match.arg(method)
  trend_type <- match.arg(trend_type)
  
  levels_hist <- suppressWarnings(as.numeric(levels_hist))
  values_hist <- suppressWarnings(as.numeric(values_hist))
  future_levels <- suppressWarnings(as.numeric(future_levels))
  
  keep_hist <- is.finite(levels_hist) & is.finite(values_hist)
  levels_hist <- levels_hist[keep_hist]
  values_hist <- values_hist[keep_hist]
  if (!length(future_levels)) return(numeric(0))
  if (!length(levels_hist)) return(rep(0, length(future_levels)))
  
  ord <- order(levels_hist)
  levels_hist <- levels_hist[ord]
  values_hist <- values_hist[ord]
  
  last_val <- tail(values_hist, 1)
  h <- length(future_levels)
  
  if (method == "freeze") {
    return(rep(as.numeric(last_val), h))
  }
  
  yts <- ts(values_hist, start = min(levels_hist), frequency = 1)
  
  if (method == "arima") {
    fit <- do.call(forecast::auto.arima, c(list(yts), arima_args))
    fc <- forecast::forecast(fit, h = h)$mean
    return(as.numeric(fc))
  }
  
  if (method == "trend") {
    fit <- StructTS(yts, type = trend_type)
    pr <- predict(fit, n.ahead = h)$pred
    return(as.numeric(pr))
  }
  
  k <- suppressWarnings(as.integer(window))[1]
  if (!is.finite(k) || k < 2L) k <- 2L
  k <- min(k, length(values_hist))
  xk <- tail(levels_hist, k)
  yk <- tail(values_hist, k)
  
  if (length(unique(xk)) < 2L) {
    slope <- 0
  } else {
    slope <- tryCatch(unname(stats::coef(stats::lm(yk ~ xk))[2]), error = function(e) NA_real_)
    if (!is.finite(slope)) slope <- 0
  }
  
  delta <- suppressWarnings(as.numeric(damping))[1]
  if (!is.finite(delta)) delta <- 0.8
  delta <- max(min(delta, 1), 0)
  
  increments <- slope * cumsum(delta ^ (seq_len(h) - 1L))
  as.numeric(last_val + increments)
}

.lookup_apc_hist_effect <- function(levels_hist, values_hist, query_levels) {
  levels_hist <- suppressWarnings(as.numeric(levels_hist))
  values_hist <- suppressWarnings(as.numeric(values_hist))
  query_levels <- suppressWarnings(as.numeric(query_levels))
  out <- rep(NA_real_, length(query_levels))
  if (!length(levels_hist) || !length(values_hist) || !length(query_levels)) return(out)
  idx <- match(query_levels, levels_hist)
  keep <- is.finite(idx)
  out[keep] <- values_hist[idx[keep]]
  out
}

forecast_gammaP <- function(gammaP_hist_df, cohorts_future,
                            method = c("freeze", "arima", "trend", "damped_trend"),
                            arima_args = list(stepwise = FALSE, approximation = FALSE, seasonal = FALSE),
                            trend_type = c("trend", "level"),
                            window = COHORT_FC_WINDOW,
                            damping = COHORT_FC_DAMPING) {
  method <- match.arg(method)
  trend_type <- match.arg(trend_type)
  h <- length(cohorts_future)
  if (h <= 0) return(tibble(cohort = integer(0), gammaP = numeric(0)))
  
  g <- gammaP_hist_df %>% dplyr::arrange(cohort)
  hist_levels <- suppressWarnings(as.numeric(g$cohort))
  hist_values <- suppressWarnings(as.numeric(g$gammaP))
  future_levels <- suppressWarnings(as.numeric(cohorts_future))
  
  fc <- .forecast_apc_tail(
    levels_hist = hist_levels,
    values_hist = hist_values,
    future_levels = future_levels,
    method = method,
    arima_args = arima_args,
    trend_type = trend_type,
    window = window,
    damping = damping
  )
  
  tibble::tibble(cohort = future_levels, gammaP = fc)
}

# Ajusta γ^P futuro según el escenario (cohortes no observadas)
adjust_gammaP_future <- function(gammaP_hist_df, gammaP_fut_df,
                                 scenario = PREV_SCENARIO,
                                 annual_rate = PREV_ANNUAL_RATE,
                                 annual_rate_down3 = PREV_ANNUAL_RATE_DOWN3,
                                 base_year = PREV_BASE_YEAR) {
  if (!nrow(gammaP_fut_df)) return(gammaP_fut_df)
  gammaP_fut_df <- gammaP_fut_df %>% arrange(cohort)
  
  if (identical(scenario, "freeze")) {
    last_val <- tail(gammaP_hist_df$gammaP[order(gammaP_hist_df$cohort)], 1)
    gammaP_fut_df$gammaP <- last_val
    
  } else if (scenario %in% c("up1pc","down1pc","down3pc")) {
    step <- if (identical(scenario, "up1pc")) abs(annual_rate) else if (identical(scenario, "down1pc")) -abs(annual_rate) else -abs(annual_rate_down3)
    k <- pmax(0, gammaP_fut_df$cohort - base_year)  # años desde 2022
    # odds_t = odds_ref * (1+step)^k  =>  logit shift = k*log(1+step)
    gammaP_fut_df$gammaP <- gammaP_fut_df$gammaP + k * log1p(step)
  }
  # "quit": no tocamos γ^P aquí; se maneja poniendo z_prev = 0 en futuro
  gammaP_fut_df
}


# ---------- Infraestructura v1 del nuevo canal PREV -> INC (stock/flow por cohorte)
prev_inc_channel_assumptions <- function() {
  tibble::tribble(
    ~assumption_id, ~label, ~description,
    "P1", "Proportional re-entry from non-current pool", "Si la prevalencia de fumadores actuales aumenta entre dos edades consecutivas de una misma cohorte, el aumento se toma proporcionalmente del pool no fumador (never + former).",
    "P2", "Finite former-smoker memory", "Solo importan los exfumadores con antigüedad de cesación hasta L_q años; por encima de ese horizonte el exceso de riesgo remanente se toma como cero.",
    "P3", "Window-start initialization", "Al comienzo de la ventana retrospectiva usada para construir estados por cohorte, el stock de exfumadores se inicializa en cero; la memoria previa a esa ventana se aproxima por truncación finita.",
    "P4", "Frozen first-period backcast in PREV", "Para periodos anteriores al primer año observado de PREV, el efecto de período se congela en su primer valor observado.",
    "P5", "Frozen last-period forward PREV level", "Para periodos posteriores al último año observado de PREV, el efecto de período se congela en su último valor observado; los cambios futuros quedan en la cohorte/backbone que se defina.",
    "P6", "Carry states above PREV support", "Para edades mayores al máximo observado en PREV, los estados de tabaquismo se transportan por cohorte sin introducir nuevos cambios no observados por encuesta."
  )
}

get_prev_inc_quit_horizon <- function(cause_id = NA_character_, sex = NA_character_, fallback = QUIT_A_I_MAX) {
  fb <- suppressWarnings(as.integer(fallback))[1]
  if (is.finite(fb) && fb > 0L) return(fb)
  
  if (exists("risk_reversal_schedule_tbl", inherits = TRUE) && exists("normalize_rr_schedule_site", inherits = TRUE)) {
    site <- tryCatch(normalize_rr_schedule_site(cause_id), error = function(e) NA_character_)
    sx <- .normalize_rr_lookup_sex(sex)
    tab <- tryCatch(risk_reversal_schedule_tbl(), error = function(e) NULL)
    if (is.data.frame(tab) && nrow(tab)) {
      sub <- tab[tab$site == site & tab$sex == sx, , drop = FALSE]
      if (!nrow(sub)) sub <- tab[tab$site == site, , drop = FALSE]
      yrs <- suppressWarnings(as.integer(sub$years_since_quit))
      yrs <- yrs[is.finite(yrs)]
      if (length(yrs)) return(max(yrs, na.rm = TRUE))
    }
  }
  
  # Ultimate fallback to QUIT_A_I_MAX (L_q from paper logic)
  as.integer(fallback %||% QUIT_A_I_MAX %||% 50L)
}

# =========================================================
# Hard-coded excess-risk reversal schedules after smoking cessation
# Target quantity:
#   w(s) = (RR(s) - 1) / (RR(0) - 1)
# Curated for model use:
# - non-negative
# - non-increasing in years since quit
# - linear interpolation between anchors, clamped at endpoints
# =========================================================

risk_reversal_schedule_tbl <- function() {
  tibble::tribble(
    ~site, ~sex, ~years_since_quit, ~w_remaining,
    "bladder", "F", 0, 1.000000,
    "bladder", "F", 4, 0.813699,
    "bladder", "F", 9, 0.682192,
    "bladder", "F", 10, 0.295890,
    "bladder", "F", 50, 0.000000,
    "bladder", "M", 0, 1.000000,
    "bladder", "M", 4, 0.803419,
    "bladder", "M", 9, 0.639667,
    "bladder", "M", 10, 0.321799,
    "bladder", "M", 50, 0.000000,
    "cervix", "F", 0, 1.000000,
    "cervix", "F", 4, 0.108696,
    "cervix", "F", 9, 0.108696,
    "cervix", "F", 10, 0.000000,
    "esophagus", "F", 0, 1.000000,
    "esophagus", "F", 5, 0.460650,
    "esophagus", "F", 10, 0.237899,
    "esophagus", "F", 20, 0.132354,
    "esophagus", "F", 50, 0.000000,
    "esophagus", "M", 0, 1.000000,
    "esophagus", "M", 5, 0.460650,
    "esophagus", "M", 10, 0.237899,
    "esophagus", "M", 20, 0.132354,
    "esophagus", "M", 50, 0.000000,
    "kidney", "F", 0, 1.000000,
    "kidney", "F", 9, 0.622642,
    "kidney", "F", 19, 0.471698,
    "kidney", "F", 20, 0.283019,
    "kidney", "F", 50, 0.000000,
    "kidney", "M", 0, 1.000000,
    "kidney", "M", 9, 0.622642,
    "kidney", "M", 19, 0.471698,
    "kidney", "M", 20, 0.283019,
    "kidney", "M", 50, 0.000000,
    "larynx", "F", 0, 1.000000,
    "larynx", "F", 2, 1.000000,
    "larynx", "F", 9, 0.578947,
    "larynx", "F", 19, 0.210526,
    "larynx", "F", 20, 0.126316,
    "larynx", "F", 50, 0.000000,
    "larynx", "M", 0, 1.000000,
    "larynx", "M", 2, 1.000000,
    "larynx", "M", 9, 0.578947,
    "larynx", "M", 19, 0.210526,
    "larynx", "M", 20, 0.126316,
    "larynx", "M", 50, 0.000000,
    "lung", "F", 0, 1.000000,
    "lung", "F", 1, 0.814000,
    "lung", "F", 5, 0.572000,
    "lung", "F", 10, 0.369000,
    "lung", "F", 15, 0.267000,
    "lung", "F", 20, 0.197000,
    "lung", "F", 50, 0.000000,
    "lung", "M", 0, 1.000000,
    "lung", "M", 1, 0.814000,
    "lung", "M", 5, 0.572000,
    "lung", "M", 10, 0.369000,
    "lung", "M", 15, 0.267000,
    "lung", "M", 20, 0.197000,
    "lung", "M", 50, 0.000000,
    "oral_pharynx", "F", 0, 1.000000,
    "oral_pharynx", "F", 2, 0.709016,
    "oral_pharynx", "F", 5, 0.479564,
    "oral_pharynx", "F", 9, 0.332879,
    "oral_pharynx", "F", 14, 0.084469,
    "oral_pharynx", "F", 15, 0.057221,
    "oral_pharynx", "F", 50, 0.000000,
    "oral_pharynx", "M", 0, 1.000000,
    "oral_pharynx", "M", 2, 0.709016,
    "oral_pharynx", "M", 5, 0.479564,
    "oral_pharynx", "M", 9, 0.332879,
    "oral_pharynx", "M", 14, 0.084469,
    "oral_pharynx", "M", 15, 0.057221,
    "oral_pharynx", "M", 50, 0.000000,
    "pancreas", "F", 0, 1.000000,
    "pancreas", "F", 9, 0.648649,
    "pancreas", "F", 10, 0.202703,
    "pancreas", "F", 20, 0.000000,
    "pancreas", "M", 0, 1.000000,
    "pancreas", "M", 9, 0.648649,
    "pancreas", "M", 10, 0.202703,
    "pancreas", "M", 20, 0.000000,
    "stomach", "F", 0, 1.000000,
    "stomach", "F", 10, 0.632653,
    "stomach", "F", 19, 0.632653,
    "stomach", "F", 20, 0.632653,
    "stomach", "F", 50, 0.000000,
    "stomach", "M", 0, 1.000000,
    "stomach", "M", 10, 0.793103,
    "stomach", "M", 19, 0.551724,
    "stomach", "M", 20, 0.198276,
    "stomach", "M", 50, 0.000000
  )
}

normalize_rr_schedule_site <- function(cause_id) {
  x <- suppressWarnings(as.character(cause_id)[1])
  if (!length(x) || is.na(x) || !nzchar(x)) return(NA_character_)
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  if (x %in% c("oralphar", "oralpharynx", "oral_pharynx", "oropharynx", "cause_oralphar")) return("oral_pharynx")
  if (x %in% c("cause_lung", "c34", "lung_res_both")) return("lung")
  if (x %in% c("c67", "cause_bladder")) return("bladder")
  if (x %in% c("c32", "cause_larynx")) return("larynx")
  if (x %in% c("c15", "cause_esophagus", "oesophagus")) return("esophagus")
  if (x %in% c("c16", "cause_stomach")) return("stomach")
  if (x %in% c("c25", "cause_pancreas")) return("pancreas")
  if (x %in% c("c64", "cause_kidney")) return("kidney")
  if (x %in% c("c53", "cause_cervix")) return("cervix")
  x
}

get_risk_reversal_w <- function(cause_id = NULL, sex = NULL, years_since_quit = 0, half_life = QUIT_HALF_LIFE) {
  t <- suppressWarnings(as.numeric(years_since_quit))
  
  # Try table lookup first
  site <- normalize_rr_schedule_site(cause_id)
  sx <- toupper(substr(as.character(sex)[1], 1, 1))
  tab <- risk_reversal_schedule_tbl()
  sub <- tab[tab$site == site & tab$sex == sx, , drop = FALSE]
  
  if (nrow(sub) > 0) {
    # Linear interpolation from table
    sub <- sub[order(sub$years_since_quit), , drop = FALSE]
    x_coords <- as.numeric(sub$years_since_quit)
    y_coords <- pmax(0, pmin(1, as.numeric(sub$w_remaining)))
    y_coords <- cummin(y_coords)
    out <- stats::approx(x = x_coords, y = y_coords, xout = t, method = "linear", rule = 2, ties = "ordered")$y
    return(pmin(pmax(as.numeric(out), 0), 1))
  }
  
  # Fallback to Exponential decay: weight = 2^(-t/hl)
  hl <- suppressWarnings(as.numeric(half_life))[1]
  if (!is.finite(hl) || hl <= 0) hl <- 2.0
  out <- exp(-log(2) * t / hl)
  pmin(pmax(as.numeric(out), 0), 1)
}

.prev_rr_schedule_vec <- function(cause_id, sex, L_q) {
  # Robust handling of horizon
  L_q_val <- suppressWarnings(as.integer(L_q))[1]
  if (is.na(L_q_val)) L_q_val <- suppressWarnings(as.integer(get_prev_inc_quit_horizon(cause_id, sex)))
  if (is.na(L_q_val)) L_q_val <- 50L
  
  L_q <- max(0L, L_q_val)
  if (L_q <= 0L) return(numeric(0))
  # Use seq(0, L_q - 1) instead of seq_len(L_q) to include year 0 (w=1)
  # This avoids the immediate jump at the 2022/2023 boundary when p_curr changes.
  out <- get_risk_reversal_w(cause_id = cause_id, sex = sex, years_since_quit = seq(0, L_q - 1))
  out <- pmax(0, pmin(1, as.numeric(out)))
  out[!is.finite(out)] <- 0
  if (length(out) > 1L) out <- cummin(out)
  out
}

.prev_extract_random_mean <- function(sr, values) {
  if (is.null(sr) || !nrow(sr)) return(tibble::tibble(value = values, eff = rep(0, length(values))))
  out <- tibble::tibble(value = values, eff = rep(0, length(values)))
  if ("ID" %in% names(sr)) {
    id <- suppressWarnings(as.integer(sr$ID))
    keep <- is.finite(id) & id >= 1L & id <= length(values)
    out$eff[id[keep]] <- suppressWarnings(as.numeric(sr$mean[keep]))
  } else {
    n <- min(nrow(sr), length(values))
    out$eff[seq_len(n)] <- suppressWarnings(as.numeric(sr$mean[seq_len(n)]))
  }
  out
}


.prev_fixed_anchor <- function(fit_prev, prev_inla = NULL) {
  if (is.null(prev_inla)) prev_inla <- tryCatch(fit_prev$.args$data, error = function(e) NULL)
  sf <- tryCatch(fit_prev$summary.fixed, error = function(e) NULL)
  if (is.null(sf) || !nrow(sf)) return(0)
  nm <- rownames(sf)
  vv <- suppressWarnings(as.numeric(sf$mean))
  if (!length(vv)) return(0)
  if (!is.null(prev_inla) && is.data.frame(prev_inla) && "inst" %in% names(prev_inla)) {
    inst_rows <- grep('^inst', nm)
    if (length(inst_rows)) {
      inst_levels_fit <- sub('^inst', '', nm[inst_rows])
      inst_tbl <- prev_inla %>%
        dplyr::mutate(.inst_chr = as.character(inst), .w_inst = if ("N" %in% names(prev_inla)) pmax(as.numeric(N), 1) else 1) %>%
        dplyr::group_by(.inst_chr) %>%
        dplyr::summarise(w = sum(.w_inst, na.rm = TRUE), .groups = 'drop')
      map_tbl <- tibble::tibble(.inst_chr = inst_levels_fit, eff = vv[inst_rows])
      joined <- inst_tbl %>% dplyr::left_join(map_tbl, by = '.inst_chr')
      ok <- is.finite(joined$w) & joined$w > 0 & is.finite(joined$eff)
      if (any(ok)) return(stats::weighted.mean(joined$eff[ok], joined$w[ok]))
    }
  }
  idx0 <- which(nm %in% c('(Intercept)', 'Intercept'))
  if (length(idx0)) return(vv[idx0[1]])
  keep <- !(nm %in% c('prev_trend_t', 'prev_trend_t2'))
  if (any(keep & is.finite(vv))) return(mean(vv[keep & is.finite(vv)], na.rm = TRUE))
  0
}

.prev_trend_component <- function(period_vec, fit_prev, prev_inla = NULL) {
  if (is.null(prev_inla)) prev_inla <- tryCatch(fit_prev$.args$data, error = function(e) NULL)
  sf <- tryCatch(fit_prev$summary.fixed, error = function(e) NULL)
  if (is.null(sf) || !nrow(sf)) return(rep(0, length(period_vec)))
  nm <- rownames(sf)
  vv <- suppressWarnings(as.numeric(sf$mean))
  center_period <- if (!is.null(prev_inla) && is.data.frame(prev_inla) && "period" %in% names(prev_inla)) {
    mean(sort(unique(suppressWarnings(as.integer(prev_inla$period)))), na.rm = TRUE)
  } else {
    0
  }
  t <- suppressWarnings(as.numeric(period_vec)) - center_period
  beta1 <- if ("prev_trend_t" %in% nm) vv[match("prev_trend_t", nm)] else 0
  beta2 <- if ("prev_trend_t2" %in% nm) vv[match("prev_trend_t2", nm)] else 0
  beta1 <- dplyr::coalesce(as.numeric(beta1), 0)
  beta2 <- dplyr::coalesce(as.numeric(beta2), 0)
  out <- beta1 * t + beta2 * (t^2)
  dplyr::coalesce(as.numeric(out), 0) # Blindaje final contra NAs en tendencia
}

.prev_hist_surface_from_fit <- function.prev_hist_surface_from_fit <- function(fit_prev, prev_inla = NULL) {
  if (is.null(prev_inla)) prev_inla <- tryCatch(fit_prev$.args$data, error = function(e) NULL)
  if (is.null(prev_inla) || !is.data.frame(prev_inla) || !nrow(prev_inla)) return(tibble::tibble())
  lp <- tryCatch(fit_prev$summary.linear.predictor, error = function(e) NULL)
  if (is.null(lp) || !nrow(lp)) return(tibble::tibble())
  n <- min(nrow(prev_inla), nrow(lp))
  df <- prev_inla[seq_len(n), , drop = FALSE]
  if (!"sex" %in% names(df)) df$sex <- NA_character_
  if (!"cohort" %in% names(df) && all(c("period", "age") %in% names(df))) df$cohort <- df$period - df$age
  df$.p_row <- plogis(suppressWarnings(as.numeric(lp$mean[seq_len(n)])))
  df$.w_row <- if ("N" %in% names(df)) pmax(suppressWarnings(as.numeric(df$N)), 1) else 1
  df %>%
    dplyr::group_by(sex, age, period, cohort) %>%
    dplyr::summarise(
      p_cur = stats::weighted.mean(.p_row, .w_row, na.rm = TRUE),
      prev_source = 'hist_fit',
      within_prev_age_support = TRUE,
      within_prev_period_support = TRUE,
      within_prev_observed_support = TRUE,
      within_prev_support = TRUE,
      prev_scenario_name = NA_character_,
      prev_scenario_applied = FALSE,
      .groups = 'drop'
    )
}

backcast_gammaP <- function(gammaP_hist_df, cohorts_past, method = c('freeze_oldest', 'trend'), trend_type = c('level', 'trend')) {
  method <- match.arg(method)
  trend_type <- match.arg(trend_type)
  cohorts_past <- sort(unique(suppressWarnings(as.integer(cohorts_past))))
  cohorts_past <- cohorts_past[is.finite(cohorts_past)]
  if (!length(cohorts_past)) return(tibble::tibble(cohort = integer(0), gammaP = numeric(0)))
  g <- gammaP_hist_df %>% dplyr::arrange(cohort)
  if (!nrow(g)) return(tibble::tibble(cohort = cohorts_past, gammaP = rep(0, length(cohorts_past))))
  if (identical(method, 'freeze_oldest') || nrow(g) < 3) {
    return(tibble::tibble(cohort = cohorts_past, gammaP = rep(g$gammaP[1], length(cohorts_past))))
  }
  rev_hist <- g %>% dplyr::arrange(dplyr::desc(cohort)) %>% dplyr::transmute(cohort = -cohort, gammaP = gammaP)
  rev_future <- sort(unique(-cohorts_past))
  rev_fc <- forecast_gammaP(rev_hist, rev_future, method = 'trend', trend_type = trend_type) %>%
    dplyr::transmute(cohort = -cohort, gammaP = gammaP) %>%
    dplyr::arrange(cohort)
  tibble::tibble(cohort = cohorts_past) %>% dplyr::left_join(rev_fc, by = 'cohort')
}

build_prev_current_surface_for_inc <- function(target_grid,
                                               fit_prev,
                                               prev_inla = NULL,
                                               sex_sel = NULL,
                                               gammaP_method = GAMMAP_METHOD,
                                               trend_type = TREND_TYPE,
                                               prev_cfg = NULL,
                                               backcast_period_mode = PREV_BACKCAST_MODE,
                                               backcast_cohort_mode = PREV_BACKCAST_COHORT_MODE,
                                               post65_mode = PREV_POST65_MODE,
                                               age_min_p = AGE_P_MIN,
                                               age_max_p = AGE_P_MAX) {
  if (is.null(prev_cfg)) prev_cfg <- make_prev_config()
  if (is.null(prev_inla)) prev_inla <- tryCatch(fit_prev$.args$data, error = function(e) NULL)
  if (is.null(prev_inla) || !is.data.frame(prev_inla) || !nrow(prev_inla)) {
    stop('build_prev_current_surface_for_inc: no pude recuperar prev_inla desde fit_prev.')
  }
  req <- c('age', 'period')
  miss <- setdiff(req, names(target_grid))
  if (length(miss)) stop('build_prev_current_surface_for_inc: faltan columnas: ', paste(miss, collapse = ', '))
  target_grid <- tibble::as_tibble(target_grid)
  if (!'sex' %in% names(target_grid)) target_grid$sex <- sex_sel %||% NA_character_
  if (!'cohort' %in% names(target_grid)) target_grid$cohort <- target_grid$period - target_grid$age
  target_grid$.row_id_prev_current <- seq_len(nrow(target_grid))
  sex_levels <- sort(unique(as.character(target_grid$sex)))
  target_grid$sex <- as.character(target_grid$sex)

  lev_age <- sort(unique(suppressWarnings(as.integer(prev_inla$age))))
  lev_per <- sort(unique(suppressWarnings(as.integer(prev_inla$period))))
  lev_coh <- sort(unique(suppressWarnings(as.integer(prev_inla$cohort))))
  age_min_obs <- min(lev_age, na.rm = TRUE)
  age_max_obs <- max(lev_age, na.rm = TRUE)
  per_min_obs <- min(lev_per, na.rm = TRUE)
  per_max_obs <- max(lev_per, na.rm = TRUE)

  age_re <- .prev_extract_random_mean(tryCatch(fit_prev$summary.random$age_id, error = function(e) NULL), lev_age) %>%
    dplyr::rename(age = value, age_eff = eff)
  per_re_hist <- .prev_extract_random_mean(tryCatch(fit_prev$summary.random$period_id, error = function(e) NULL), lev_per) %>%
    dplyr::rename(period = value, period_eff = eff)
  coh_hist <- .prev_extract_random_mean(tryCatch(fit_prev$summary.random$cohort_id, error = function(e) NULL), lev_coh) %>%
    dplyr::rename(cohort = value, cohort_eff = eff) %>%
    dplyr::arrange(cohort)

  # cohorts needed on PREV support ages only
  support_target <- target_grid %>% dplyr::filter(age >= age_min_p, age <= age_max_p)
  cohorts_needed <- sort(unique(suppressWarnings(as.integer(support_target$cohort))))
  if (!length(cohorts_needed)) cohorts_needed <- sort(unique(suppressWarnings(as.integer(target_grid$cohort))))
  coh_future_needed <- cohorts_needed[cohorts_needed > max(coh_hist$cohort, na.rm = TRUE)]
  coh_past_needed <- cohorts_needed[cohorts_needed < min(coh_hist$cohort, na.rm = TRUE)]
  coh_future <- forecast_gammaP(coh_hist %>% dplyr::transmute(cohort, gammaP = cohort_eff), coh_future_needed,
                                method = gammaP_method, trend_type = trend_type)
  coh_future <- adjust_gammaP_future(coh_hist %>% dplyr::transmute(cohort, gammaP = cohort_eff), coh_future,
                                     scenario = if (!is.null(prev_cfg$backbone) && identical(prev_cfg$backbone, 'forecast')) prev_cfg$scenario else 'freeze',
                                     annual_rate = prev_cfg$annual_rate,
                                     annual_rate_down3 = prev_cfg$annual_rate_down3,
                                     base_year = prev_cfg$base_year) %>%
    dplyr::rename(cohort_eff = gammaP)
  coh_past <- backcast_gammaP(coh_hist %>% dplyr::transmute(cohort, gammaP = cohort_eff), coh_past_needed,
                              method = backcast_cohort_mode, trend_type = trend_type) %>%
    dplyr::rename(cohort_eff = gammaP)
  coh_all <- dplyr::bind_rows(coh_hist, coh_future, coh_past) %>% dplyr::distinct(cohort, .keep_all = TRUE)

  hist_surface <- .prev_hist_surface_from_fit(fit_prev, prev_inla)
  fix_anchor <- .prev_fixed_anchor(fit_prev, prev_inla)
  first_per_eff <- per_re_hist$period_eff[match(per_min_obs, per_re_hist$period)]
  last_per_eff <- per_re_hist$period_eff[match(per_max_obs, per_re_hist$period)]
  first_per_eff <- ifelse(is.finite(first_per_eff), first_per_eff, 0)
  last_per_eff <- ifelse(is.finite(last_per_eff), last_per_eff, 0)

  support_grid <- target_grid %>%
    dplyr::filter(age >= age_min_p, age <= age_max_p) %>%
    dplyr::left_join(hist_surface, by = c('sex', 'age', 'period', 'cohort'))

  miss_support <- is.na(support_grid$p_cur)
  if (any(miss_support)) {
    fill_df <- support_grid[miss_support, , drop = FALSE] %>%
      dplyr::left_join(age_re, by = 'age') %>%
      dplyr::left_join(coh_all, by = 'cohort') %>%
      dplyr::left_join(per_re_hist, by = 'period') %>%
      dplyr::mutate(
        period_eff = dplyr::if_else(is.na(period_eff) & period < per_min_obs, as.numeric(first_per_eff), as.numeric(period_eff)),
        period_eff = dplyr::if_else(is.na(period_eff) & period > per_max_obs, as.numeric(last_per_eff), as.numeric(period_eff)),
        age_eff = dplyr::coalesce(age_eff, 0),
        cohort_eff = dplyr::coalesce(cohort_eff, 0),
        period_eff = dplyr::coalesce(period_eff, 0),
        prev_trend_component = .prev_trend_component(period, fit_prev, prev_inla),
        eta_prev = fix_anchor + prev_trend_component + age_eff + period_eff + cohort_eff,
        p_cur = plogis(eta_prev),
        prev_source = dplyr::case_when(
          period < per_min_obs ~ 'backcast_period',
          period > per_max_obs ~ 'forecast_period',
          TRUE ~ 'support_model_fill'
        ),
        within_prev_age_support = TRUE,
        within_prev_period_support = dplyr::between(period, per_min_obs, per_max_obs),
        within_prev_observed_support = within_prev_age_support & within_prev_period_support,
        within_prev_support = within_prev_observed_support,
        prev_scenario_name = NA_character_,
        prev_scenario_applied = FALSE
      ) %>%
      dplyr::select(.row_id_prev_current, p_cur, prev_source, within_prev_age_support, within_prev_period_support,
                    within_prev_observed_support, within_prev_support, prev_scenario_name, prev_scenario_applied)
    support_grid <- support_grid %>% dplyr::left_join(fill_df, by = '.row_id_prev_current', suffix = c('', '.fill')) %>%
      dplyr::mutate(
        p_cur = dplyr::coalesce(p_cur, p_cur.fill),
        prev_source = dplyr::coalesce(prev_source, prev_source.fill),
        within_prev_age_support = dplyr::coalesce(within_prev_age_support, within_prev_age_support.fill),
        within_prev_period_support = dplyr::coalesce(within_prev_period_support, within_prev_period_support.fill),
        within_prev_observed_support = dplyr::coalesce(within_prev_observed_support, within_prev_observed_support.fill),
        within_prev_support = dplyr::coalesce(within_prev_support, within_prev_support.fill),
        prev_scenario_name = dplyr::coalesce(prev_scenario_name, prev_scenario_name.fill),
        prev_scenario_applied = dplyr::coalesce(prev_scenario_applied, prev_scenario_applied.fill)
      ) %>%
      dplyr::select(-dplyr::any_of(c('p_cur.fill', 'prev_source.fill', 'within_prev_age_support.fill', 'within_prev_period_support.fill',
                                     'within_prev_observed_support.fill', 'within_prev_support.fill', 'prev_scenario_name.fill', 'prev_scenario_applied.fill')))
  }

  if (!is.null(prev_cfg) && identical(prev_cfg$axis, 'period')) {
    scen_prev <- tryCatch(normalize_prev_scenario_name(prev_cfg$scenario), error = function(e) 'freeze')
    k_prev <- pmax(0L, suppressWarnings(as.integer(support_grid$period)) - suppressWarnings(as.integer(prev_cfg$base_year))[1])
    idx_prev <- is.finite(k_prev) & k_prev > 0L
    if (any(idx_prev) && !identical(scen_prev, 'freeze')) {
      p_base <- pmin(pmax(as.numeric(support_grid$p_cur), 0), 1)
      if (identical(scen_prev, 'up1pc')) {
        p_base[idx_prev] <- pmin(pmax(p_base[idx_prev] * (1 + abs(prev_cfg$annual_rate))^k_prev[idx_prev], 0), 1)
      } else if (identical(scen_prev, 'down1pc')) {
        p_base[idx_prev] <- pmin(pmax(p_base[idx_prev] * (1 - abs(prev_cfg$annual_rate))^k_prev[idx_prev], 0), 1)
      } else if (identical(scen_prev, 'down3pc')) {
        p_base[idx_prev] <- pmin(pmax(p_base[idx_prev] * (1 - abs(prev_cfg$annual_rate_down3))^k_prev[idx_prev], 0), 1)
      } else if (identical(scen_prev, 'quit')) {
        p_base[idx_prev] <- 0
      }
      support_grid$p_cur <- p_base
      support_grid$prev_scenario_name <- dplyr::if_else(idx_prev, as.character(scen_prev), support_grid$prev_scenario_name)
      support_grid$prev_scenario_applied <- dplyr::if_else(idx_prev, TRUE, support_grid$prev_scenario_applied %||% FALSE)
    }
  }


  out <- target_grid %>% dplyr::left_join(
    support_grid %>% dplyr::select(.row_id_prev_current, p_cur, prev_source,
                                   within_prev_age_support, within_prev_period_support,
                                   within_prev_observed_support, within_prev_support,
                                   prev_scenario_name, prev_scenario_applied),
    by = '.row_id_prev_current'
  )

  if (identical(post65_mode, 'carry_states')) {
    over_idx <- which(out$age > age_max_p)
    if (length(over_idx)) {
      anchors <- out %>%
        dplyr::filter(age <= age_max_p, is.finite(p_cur)) %>%
        dplyr::group_by(sex, cohort) %>%
        dplyr::arrange(age, .by_group = TRUE) %>%
        dplyr::slice_tail(n = 1) %>%
        dplyr::ungroup() %>%
        dplyr::select(sex, cohort, p_cur_anchor = p_cur)
      out <- out %>% dplyr::left_join(anchors, by = c('sex', 'cohort')) %>%
        dplyr::mutate(
          p_cur = dplyr::if_else(age > age_max_p & is.na(p_cur), p_cur_anchor, p_cur),
          prev_source = dplyr::if_else(age > age_max_p & is.na(prev_source) & !is.na(p_cur_anchor), 'carried_post65', prev_source),
          within_prev_age_support = dplyr::if_else(age > age_max_p & !is.na(p_cur_anchor), FALSE, within_prev_age_support),
          within_prev_observed_support = dplyr::if_else(age > age_max_p & !is.na(p_cur_anchor), FALSE, within_prev_observed_support),
          within_prev_support = dplyr::if_else(age > age_max_p & !is.na(p_cur_anchor), FALSE, within_prev_support)
        ) %>%
        dplyr::select(-p_cur_anchor)
    }
  }

  out <- out %>%
    dplyr::mutate(
      p_cur = pmin(pmax(as.numeric(p_cur), 0), 1),
      prev_source = dplyr::coalesce(prev_source, 'unfilled'),
      within_prev_age_support = dplyr::coalesce(within_prev_age_support, age >= age_min_p & age <= age_max_p),
      within_prev_period_support = dplyr::coalesce(within_prev_period_support, dplyr::between(period, per_min_obs, per_max_obs)),
      within_prev_observed_support = dplyr::coalesce(within_prev_observed_support, within_prev_age_support & within_prev_period_support),
      within_prev_support = dplyr::coalesce(within_prev_support, within_prev_observed_support),
      prev_scenario_name = dplyr::coalesce(prev_scenario_name, NA_character_),
      prev_scenario_applied = dplyr::coalesce(prev_scenario_applied, FALSE)
    )

  rem_na <- which(!is.finite(out$p_cur))
  if (length(rem_na)) {
    fallback_df <- out[rem_na, , drop = FALSE] %>%
      dplyr::left_join(age_re, by = 'age') %>%
      dplyr::left_join(coh_all, by = 'cohort') %>%
      dplyr::left_join(per_re_hist, by = 'period') %>%
      dplyr::mutate(
        period_eff = dplyr::if_else(is.na(period_eff) & period < per_min_obs, as.numeric(first_per_eff), as.numeric(period_eff)),
        period_eff = dplyr::if_else(is.na(period_eff) & period > per_max_obs, as.numeric(last_per_eff), as.numeric(period_eff)),
        age_eff = dplyr::coalesce(age_eff, 0),
        cohort_eff = dplyr::coalesce(cohort_eff, 0),
        period_eff = dplyr::coalesce(period_eff, 0),
        prev_trend_component = .prev_trend_component(period, fit_prev, prev_inla),
        eta_prev = fix_anchor + prev_trend_component + age_eff + period_eff + cohort_eff,
        p_cur_fill = plogis(eta_prev),
        prev_source_fill = dplyr::case_when(
          age > age_max_p ~ 'carried_post65',
          period < per_min_obs ~ 'backcast_period',
          period > per_max_obs ~ 'forecast_period',
          TRUE ~ 'support_model_fill'
        ),
        within_prev_age_support_fill = age >= age_min_p & age <= age_max_p,
        within_prev_period_support_fill = dplyr::between(period, per_min_obs, per_max_obs),
        within_prev_observed_support_fill = within_prev_age_support_fill & within_prev_period_support_fill,
        within_prev_support_fill = within_prev_observed_support_fill
      ) %>%
      dplyr::select(.row_id_prev_current, p_cur_fill, prev_source_fill,
                    within_prev_age_support_fill, within_prev_period_support_fill,
                    within_prev_observed_support_fill, within_prev_support_fill)

    out <- out %>%
      dplyr::left_join(fallback_df, by = '.row_id_prev_current') %>%
      dplyr::mutate(
        p_cur = dplyr::coalesce(p_cur, p_cur_fill),
        prev_source = dplyr::coalesce(prev_source, prev_source_fill),
        within_prev_age_support = dplyr::coalesce(within_prev_age_support, within_prev_age_support_fill),
        within_prev_period_support = dplyr::coalesce(within_prev_period_support, within_prev_period_support_fill),
        within_prev_observed_support = dplyr::coalesce(within_prev_observed_support, within_prev_observed_support_fill),
        within_prev_support = dplyr::coalesce(within_prev_support, within_prev_support_fill)
      ) %>%
      dplyr::select(-dplyr::any_of(c('p_cur_fill','prev_source_fill','within_prev_age_support_fill',
                                     'within_prev_period_support_fill','within_prev_observed_support_fill',
                                     'within_prev_support_fill')))
  }

  if (isTRUE(any(!is.finite(out$p_cur)))) {
    bad <- out %>% dplyr::filter(!is.finite(p_cur)) %>% dplyr::slice_head(n = 10)
    stop(
      'build_prev_current_surface_for_inc: quedaron p_cur no finitos tras el fill. Ejemplos: ',
      paste(sprintf('[sex=%s age=%s period=%s cohort=%s src=%s]', bad$sex, bad$age, bad$period, bad$cohort, bad$prev_source), collapse = ' ')
    )
  }

  # (Scenario applied at support_grid level for consistent stock-model history)


  out %>%
    dplyr::arrange(.row_id_prev_current) %>%
    dplyr::select(-.row_id_prev_current)
}

.calculate_stock_former_q_eff <- function(p_cur_vec, w_vec, L_q, return_all = FALSE) {
  L_q <- max(0L, suppressWarnings(as.integer(L_q))[1])
  if (!is.finite(L_q)) L_q <- 0L
  
  n <- length(p_cur_vec)
  if (n == 0) return(numeric(0))
  
  former <- if (L_q > 0) matrix(0, nrow = n, ncol = L_q) else NULL
  q_eff <- rep(NA_real_, n)
  
  p0 <- pmin(pmax(as.numeric(p_cur_vec[1]), 0), 1)
  q_eff[1] <- p0
  
  if (n >= 2) {
    for (i in 2:n) {
      p_prev <- pmin(pmax(as.numeric(p_cur_vec[i - 1]), 0), 1)
      p_now  <- pmin(pmax(as.numeric(p_cur_vec[i]), 0), 1)
      dlt <- p_now - p_prev
      
      aged_former <- if (L_q > 0) {
        c(max(-dlt, 0), former[i - 1, seq_len(max(L_q - 1, 1))])[seq_len(L_q)]
      } else numeric(0)
      
      if (isTRUE(dlt > 0)) {
        prev_sum_former <- if (L_q > 0) sum(former[i - 1, ], na.rm = TRUE) else 0
        never_prev <- max(0, 1 - p_prev - prev_sum_former)
        noncurrent_prev <- never_prev + sum(aged_former, na.rm = TRUE)
        noncurrent_now <- max(0, 1 - p_now)
        eta <- if (is.finite(noncurrent_prev) && noncurrent_prev > 1e-12) min(max(noncurrent_now / noncurrent_prev, 0), 1) else 0
        if (L_q > 0) former[i, ] <- aged_former * eta
      } else {
        if (L_q > 0) former[i, ] <- aged_former
      }
      sum_w_former <- if (L_q > 0) sum(w_vec * former[i, ], na.rm = TRUE) else 0
      q_eff[i] <- p_now + sum_w_former
    }
  }
  
  if (isTRUE(return_all)) {
    p_former_total <- if (L_q > 0) rowSums(former, na.rm = TRUE) else rep(0, n)
    p_never <- pmax(0, 1 - p_cur_vec - p_former_total)
    p_prev_vec <- c(as.numeric(p_cur_vec[1]), as.numeric(p_cur_vec[1:(n-1)]))
    delta_p <- as.numeric(p_cur_vec) - p_prev_vec
    return(tibble::tibble(
      q_eff = pmin(pmax(q_eff, 0), 1),
      delta_p_cur = delta_p,
      p_former_total = p_former_total,
      p_never = p_never
    ))
  }
  
  q_eff <- pmin(pmax(q_eff, 0), 1)
  q_eff[!is.finite(q_eff)] <- p_cur_vec[!is.finite(q_eff)] # Fallback to current if error
  q_eff
}

build_prev_rr_offset_stock_for_inc <- function(df_inc_grid,
                                               fit_prev,
                                               cause_id,
                                               rr_inc,
                                               prev_inla = NULL,
                                               sex_sel = NULL,
                                               gammaP_method = GAMMAP_METHOD,
                                               trend_type = TREND_TYPE,
                                               prev_cfg = NULL,
                                               age_min_p = AGE_P_MIN,
                                               age_max_p = AGE_P_MAX,
                                               backcast_period_mode = PREV_BACKCAST_MODE,
                                               backcast_cohort_mode = PREV_BACKCAST_COHORT_MODE,
                                               post65_mode = PREV_POST65_MODE,
                                               quit_horizon_years = PREV_INC_MAX_QUIT_YEARS,
                                               return_internal = FALSE) {
  req <- c('age', 'period')
  miss <- setdiff(req, names(df_inc_grid))
  if (length(miss)) stop('build_prev_rr_offset_stock_for_inc: faltan columnas: ', paste(miss, collapse = ', '))
  if (is.null(prev_cfg)) prev_cfg <- make_prev_config()
  grid <- tibble::as_tibble(df_inc_grid)
  if (!'sex' %in% names(grid)) grid$sex <- sex_sel %||% NA_character_
  if (!'cohort' %in% names(grid)) grid$cohort <- grid$period - grid$age
  grid$.row_id_prev_stock <- seq_len(nrow(grid))

  res_list <- lapply(split(grid, grid$sex %||% 'NA', drop = TRUE), function(gsx) {
    sx <- as.character(gsx$sex[1])
    L_q <- get_prev_inc_quit_horizon(cause_id = cause_id, sex = sx, fallback = quit_horizon_years)
    rr_use <- suppressWarnings(as.numeric(rr_inc[1]))
    if (length(rr_inc) > 1 && !is.null(names(rr_inc))) {
      rr_try <- suppressWarnings(as.numeric(rr_inc[as.character(sx)]))
      if (is.finite(rr_try)) rr_use <- rr_try
    }
    if (!is.finite(rr_use) || rr_use <= 1) rr_use <- 2.0

    cohort_window <- gsx %>%
      dplyr::group_by(sex, cohort) %>%
      dplyr::summarise(
        # age_min_p is usually 20. We can't go below that.
        age_min_need = max(age_min_p, min(age, na.rm = TRUE) - L_q, na.rm = TRUE),
        age_max_need = max(age, na.rm = TRUE), 
        .groups = 'drop'
      ) %>%
      dplyr::filter(is.finite(age_min_need), is.finite(age_max_need))
    internal_grid <- cohort_window %>%
      dplyr::rowwise() %>%
      dplyr::do({
        tibble::tibble(sex = .$sex, cohort = .$cohort, age = seq.int(.$age_min_need, .$age_max_need))
      }) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(period = cohort + age)

    prev_surface <- build_prev_current_surface_for_inc(
      target_grid = internal_grid,
      fit_prev = fit_prev,
      prev_inla = prev_inla,
      sex_sel = sx,
      gammaP_method = gammaP_method,
      trend_type = trend_type,
      prev_cfg = prev_cfg,
      backcast_period_mode = backcast_period_mode,
      backcast_cohort_mode = backcast_cohort_mode,
      post65_mode = post65_mode,
      age_min_p = age_min_p,
      age_max_p = age_max_p
    )

    w_vec <- .prev_rr_schedule_vec(cause_id = cause_id, sex = sx, L_q = L_q)
    
    state_all <- prev_surface %>%
      dplyr::arrange(sex, cohort, age) %>%
      dplyr::group_by(sex, cohort) %>%
      dplyr::mutate(
        state_info = .calculate_stock_former_q_eff(p_cur, w_vec, L_q, return_all = TRUE)
      ) %>%
      tidyr::unpack(state_info) %>%
      dplyr::mutate(
        z_prev = log1p(q_eff * (rr_use - 1)),
        offset_prev_rr = log1p(q_eff * (rr_use - 1)),
        rr_inc = rr_use,
        quit_horizon_years = L_q
      ) %>%
      dplyr::ungroup()

    surface_out <- gsx %>%
      dplyr::left_join(state_all %>% dplyr::select(dplyr::any_of(c('sex', 'age', 'period', 'cohort', 'q_eff', 'z_prev', 'p_cur', 'delta_p_cur', 'quit_flow', 'p_never', 'p_former_total', 'noncurrent_rescale', 'offset_prev_rr', 'rr_inc', 'quit_horizon_years', 'prev_source', 'within_prev_support'))),
                       by = c('sex', 'age', 'period', 'cohort'))
    surface_out
  })

  dplyr::bind_rows(res_list) %>% dplyr::arrange(.row_id_prev_stock) %>% dplyr::select(-.row_id_prev_stock)
}

# ---------- Anti-"esquinas" (pesos de cobertura para RW1 de cohorte)
cohort_coverage_weights <- function(df, coh_var = "cohort", age_var = "age") {
  df %>% distinct(.data[[coh_var]], .data[[age_var]]) %>% count(.data[[coh_var]], name = "n_age") %>%
    mutate(w = n_age / max(n_age))
}
make_Q_rw1_weighted <- function(nC, w) {
  D <- diff(diag(nC), differences = 1); W <- diag(pmax(w[-1], 1e-6), nrow = nC - 1)
  t(D) %*% W %*% D
}

# ---------- Garantizar que existe variable exposición
ensure_exposure <- function(df) {
  # Busca cualquier columna plausible de exposición y la normaliza a 'exposure'
  exp_cols <- intersect(c("exposure","exposure.x","exposure.y","poblacion","E"), names(df))
  if (length(exp_cols) == 0)
    stop("No encuentro columna de exposición en el data.frame pasado.")
  df %>%
    dplyr::mutate(exposure = dplyr::coalesce(!!!rlang::syms(exp_cols))) %>%
    dplyr::select(-dplyr::any_of(setdiff(exp_cols, "exposure")))
}

# ---- Calibración del lag Incidencia -> Mortalidad
choose_inc_lag <- function(mort_hist, inc_rates_hist,
                           ages_use = 40:89, Lmax = 15, Da = 0) {
  mh <- mort_hist %>% dplyr::filter(age %in% ages_use)
  
  best <- list(dev = Inf, L = 0)
  for (L in 0:Lmax) {
    df <- mh %>%
      dplyr::transmute(ageL = age - Da, perL = period - L, sex, deaths, exposure) %>%
      dplyr::left_join(inc_rates_hist %>% dplyr::rename(rate = rate_hat),
                       by = c("ageL" = "age", "perL" = "period", "sex")) %>%
      dplyr::filter(is.finite(rate), rate > 0)
    
    if (nrow(df) < 200) next
    
    fit <- stats::glm(deaths ~ 1 + offset(log(exposure) + log(rate)),
                      family = stats::poisson(), data = df)
    
    if (fit$deviance < best$dev) best <- list(dev = fit$deviance, L = L)
  }
  best$L
}

# ---- Crear tendencias polinómicas
make_trend_vars <- function(period, center, degree = 1, prefix = "") {
  t  <- if (degree >= 1) (period - center) else 0
  t2 <- if (degree >= 2) (period - center)^2 else 0
  tibble(
    !!paste0(prefix,"trend_t")  := t,
    !!paste0(prefix,"trend_t2") := t2
  )
}

apply_trend_scenario_future <- function(period, last_hist, center, degree = 1,
                                        scenario = c("freeze","continue","delta"),
                                        delta = 0, prefix = "") {
  scenario <- match.arg(scenario)
  base <- dplyr::case_when(
    scenario == "freeze"   ~ make_trend_vars(last_hist, center, degree, prefix),
    TRUE                   ~ make_trend_vars(period,    center, degree, prefix) # "continue"/"delta"
  )
  off <- if (scenario == "delta") log1p(delta) * pmax(0, period - last_hist) else 0
  base[[paste0(prefix,"tech_offset")]] <- off
  base
}

# Pronóstico del efecto APC (cohorte/período) y descomposición del offset futuro
# en dos piezas:
#   signal   = pronóstico futuro del componente centrado
#   recenter = ajuste constante respecto del último histórico (=-último histórico centrado)
# de modo que: offset = signal + recenter
build_coef_fc_components <- function(summary_df, levels_vec, fut_levels,
                                     method = c("freeze", "arima", "trend", "damped_trend"),
                                     ref_levels = pmin(pmax(fut_levels, min(levels_vec, na.rm = TRUE)), max(levels_vec, na.rm = TRUE)),
                                     trend_type = TREND_TYPE,
                                     window = COHORT_FC_WINDOW,
                                     damping = COHORT_FC_DAMPING) {
  method <- match.arg(method)
  h <- length(fut_levels)
  if (h <= 0) {
    return(list(signal = rep(0, 0), recenter = rep(0, 0), offset = rep(0, 0), ref = numeric(0), fc = numeric(0)))
  }
  
  ord <- order(levels_vec[summary_df$ID])
  hist_levels <- suppressWarnings(as.numeric(levels_vec[summary_df$ID][ord]))
  hist_values <- suppressWarnings(as.numeric(summary_df$mean[ord]))
  hist_values <- hist_values - mean(hist_values)
  
  fut_levels <- suppressWarnings(as.numeric(fut_levels))
  ref_levels <- suppressWarnings(as.numeric(ref_levels))
  
  signal <- .lookup_apc_hist_effect(hist_levels, hist_values, fut_levels)
  need_hi <- is.na(signal) & is.finite(fut_levels) & fut_levels > max(hist_levels, na.rm = TRUE)
  need_lo <- is.na(signal) & is.finite(fut_levels) & fut_levels < min(hist_levels, na.rm = TRUE)
  
  if (any(need_hi)) {
    hi_levels <- sort(unique(fut_levels[need_hi]))
    hi_fc <- .forecast_apc_tail(
      levels_hist = hist_levels,
      values_hist = hist_values,
      future_levels = hi_levels,
      method = method,
      trend_type = trend_type,
      window = window,
      damping = damping
    )
    signal[need_hi] <- hi_fc[match(fut_levels[need_hi], hi_levels)]
  }
  if (any(need_lo)) {
    signal[need_lo] <- hist_values[1]
  }
  
  ref <- .lookup_apc_hist_effect(hist_levels, hist_values, ref_levels)
  ref[is.na(ref) & is.finite(ref_levels) & ref_levels > max(hist_levels, na.rm = TRUE)] <- tail(hist_values, 1)
  ref[is.na(ref) & is.finite(ref_levels) & ref_levels < min(hist_levels, na.rm = TRUE)] <- hist_values[1]
  ref[is.na(ref)] <- tail(hist_values, 1)
  
  recenter <- -ref
  offset <- signal + recenter
  list(signal = as.numeric(signal), recenter = as.numeric(recenter), offset = as.numeric(offset), ref = as.numeric(ref), fc = as.numeric(signal))
}

# Compatibilidad hacia atrás: devuelve el offset total legacy.
build_coef_fc_offset <- function(summary_df, levels_vec, fut_levels,
                                 method = c("freeze", "arima", "trend", "damped_trend"),
                                 ref_levels = pmin(pmax(fut_levels, min(levels_vec, na.rm = TRUE)), max(levels_vec, na.rm = TRUE)),
                                 trend_type = TREND_TYPE,
                                 window = COHORT_FC_WINDOW,
                                 damping = COHORT_FC_DAMPING) {
  parts <- build_coef_fc_components(
    summary_df = summary_df,
    levels_vec = levels_vec,
    fut_levels = fut_levels,
    method = method,
    ref_levels = ref_levels,
    trend_type = trend_type,
    window = window,
    damping = damping
  )
  as.numeric(parts$offset)
}

apply_mort_cohort_fc_posthoc <- function(pred_df,
                                         summary_df,
                                         levels_vec,
                                         target_levels,
                                         ref_levels,
                                         enabled = MORT_COHORT_FC_ON,
                                         method = COHORT_FC_METHOD,
                                         trend_type = TREND_TYPE,
                                         window = COHORT_FC_WINDOW,
                                         damping = COHORT_FC_DAMPING) {
  if (!isTRUE(enabled) || !is.data.frame(pred_df) || !nrow(pred_df)) return(pred_df)
  if (is.null(summary_df) || !is.data.frame(summary_df) || !nrow(summary_df)) return(pred_df)
  
  parts <- build_coef_fc_components(
    summary_df = summary_df,
    levels_vec = levels_vec,
    fut_levels = target_levels,
    method = method,
    ref_levels = ref_levels,
    trend_type = trend_type,
    window = window,
    damping = damping
  )
  
  out <- pred_df
  hist_flag <- if ("hist_flag" %in% names(out)) isTRUE(FALSE) | out$hist_flag else rep(FALSE, nrow(out))
  adj <- exp(ifelse(hist_flag, 0, dplyr::coalesce(parts$offset, 0)))
  
  for (nm in c("mu_hat", "mu_lwr", "mu_upr", "mu_hat_legacy")) {
    if (nm %in% names(out)) out[[nm]] <- as.numeric(out[[nm]]) * adj
  }
  for (nm in c("eta_lp_hat", "eta_lp_lwr", "eta_lp_upr")) {
    if (nm %in% names(out)) out[[nm]] <- as.numeric(out[[nm]]) + ifelse(hist_flag, 0, dplyr::coalesce(parts$offset, 0))
  }
  if (all(c("mu_hat", "offset_total") %in% names(out))) out$eta_resid_hat <- log(pmax(as.numeric(out$mu_hat), 1e-12)) - as.numeric(out$offset_total)
  if (all(c("mu_lwr", "offset_total") %in% names(out))) out$eta_resid_lwr <- log(pmax(as.numeric(out$mu_lwr), 1e-12)) - as.numeric(out$offset_total)
  if (all(c("mu_upr", "offset_total") %in% names(out))) out$eta_resid_upr <- log(pmax(as.numeric(out$mu_upr), 1e-12)) - as.numeric(out$offset_total)
  
  out$cohort_fc_target <- suppressWarnings(as.numeric(target_levels))
  out$cohort_fc_ref_level <- suppressWarnings(as.numeric(ref_levels))
  out$cohort_fc_signal <- parts$signal
  out$cohort_fc_ref <- parts$ref
  out$cohort_fc_offset <- ifelse(hist_flag, 0, dplyr::coalesce(parts$offset, 0))
  out$cohort_fc_adj <- adj
  out$cohort_fc_method <- method
  out
}

apply_coef_fc_recenter_lock <- function(recenter,
                                        mode = INC_COEF_FC_RECENTER_LOCK_MODE,
                                        fixed_value = INC_COEF_FC_RECENTER_LOCK_VALUE) {
  recenter <- as.numeric(recenter)
  if (!length(recenter)) return(recenter)
  mode <- as.character(mode)[1]
  if (is.na(mode) || !nzchar(mode)) mode <- "none"
  mode <- tolower(mode)
  if (identical(mode, "none")) return(recenter)
  if (identical(mode, "zero")) return(rep(0, length(recenter)))
  if (identical(mode, "fixed")) {
    val <- suppressWarnings(as.numeric(fixed_value))[1]
    if (!is.finite(val)) val <- 0
    return(rep(val, length(recenter)))
  }
  warning(sprintf("apply_coef_fc_recenter_lock: modo desconocido '%s'; se deja sin cambios.", mode))
  recenter
}


apply_coef_fc_lock <- function(off,
                               mode = INC_COEF_FC_LOCK_MODE,
                               fixed_value = INC_COEF_FC_LOCK_VALUE) {
  off <- as.numeric(off)
  if (!length(off)) return(off)
  mode <- as.character(mode)[1]
  if (is.na(mode) || !nzchar(mode)) mode <- "none"
  mode <- tolower(mode)
  if (identical(mode, "none")) return(off)
  if (identical(mode, "zero")) return(rep(0, length(off)))
  if (identical(mode, "fixed")) {
    val <- suppressWarnings(as.numeric(fixed_value))[1]
    if (!is.finite(val)) val <- 0
    return(rep(val, length(off)))
  }
  warning(sprintf("apply_coef_fc_lock: modo desconocido '%s'; se deja sin cambios.", mode))
  off
}


apply_inc_coef_fc_posthoc_lock <- function(grid_all,
                                           last_hist_year,
                                           mode = INC_COEF_FC_POSTHOC_LOCK_MODE,
                                           fixed_value = INC_COEF_FC_POSTHOC_LOCK_VALUE) {
  if (!is.data.frame(grid_all) || !nrow(grid_all) || !("coef_fc_offset_I" %in% names(grid_all))) {
    return(grid_all)
  }
  mode <- as.character(mode)[1]
  if (is.na(mode) || !nzchar(mode)) mode <- "none"
  mode <- tolower(mode)
  if (identical(mode, "none")) {
    grid_all$coef_fc_offset_I_raw <- suppressWarnings(as.numeric(grid_all$coef_fc_offset_I))
    grid_all$coef_fc_offset_I_effective <- suppressWarnings(as.numeric(grid_all$coef_fc_offset_I))
    grid_all$coef_fc_posthoc_adj <- 1
    grid_all$coef_fc_posthoc_lock_mode <- "none"
    return(grid_all)
  }
  out <- grid_all
  old <- suppressWarnings(as.numeric(out$coef_fc_offset_I))
  eff <- old
  idx_future <- is.finite(suppressWarnings(as.numeric(out$period))) & suppressWarnings(as.numeric(out$period)) > suppressWarnings(as.numeric(last_hist_year))[1]
  if (identical(mode, "zero")) {
    eff[idx_future] <- 0
  } else if (identical(mode, "fixed")) {
    val <- suppressWarnings(as.numeric(fixed_value))[1]
    if (!is.finite(val)) val <- 0
    eff[idx_future] <- val
  } else {
    warning(sprintf("apply_inc_coef_fc_posthoc_lock: modo desconocido '%s'; se deja sin cambios.", mode))
    mode <- "none"
  }
  adj <- exp(dplyr::coalesce(eff, 0) - dplyr::coalesce(old, 0))
  for (nm in c("mu_hat","mu_lwr","mu_upr","rate_hat","rate_lwr","rate_upr")) {
    if (nm %in% names(out)) out[[nm]] <- suppressWarnings(as.numeric(out[[nm]])) * adj
  }
  out$coef_fc_offset_I_raw <- old
  out$coef_fc_offset_I_effective <- eff
  out$coef_fc_offset_I <- eff
  out$coef_fc_posthoc_adj <- adj
  out$coef_fc_posthoc_lock_mode <- mode
  out
}



# ---- Diagnóstico explícito del borde en incidencia futura ----
make_incidence_border_diag <- function(df_future, lev_per, lev_coh, hist_df = NULL, last_hist_year = NULL) {
  if (!is.data.frame(df_future) || !nrow(df_future)) {
    return(tibble::tibble())
  }
  if (!length(lev_per) || !length(lev_coh)) {
    stop("make_incidence_border_diag: lev_per y lev_coh deben tener longitud positiva.")
  }

  per_min <- min(lev_per, na.rm = TRUE)
  per_max <- max(lev_per, na.rm = TRUE)
  coh_min <- min(lev_coh, na.rm = TRUE)
  coh_max <- max(lev_coh, na.rm = TRUE)
  last_hist_year <- suppressWarnings(as.integer(last_hist_year))[1]
  if (!is.finite(last_hist_year)) last_hist_year <- per_max

  out <- df_future %>%
    dplyr::mutate(
      period_raw = suppressWarnings(as.integer(period)),
      cohort_raw = suppressWarnings(as.integer(cohort)),
      mapped_period = pmin(pmax(period_raw, per_min), per_max),
      mapped_cohort = pmin(pmax(cohort_raw, coh_min), coh_max),
      period_is_clamped = is.finite(period_raw) & abs(period_raw - mapped_period) > 0L,
      cohort_is_edge = is.finite(cohort_raw) & abs(cohort_raw - mapped_cohort) > 0L,
      cohort_clamp_low = is.finite(cohort_raw) & cohort_raw < coh_min,
      cohort_clamp_high = is.finite(cohort_raw) & cohort_raw > coh_max,
      period_shift = suppressWarnings(as.integer(period_raw - mapped_period)),
      cohort_shift = suppressWarnings(as.integer(cohort_raw - mapped_cohort)),
      horizon = suppressWarnings(as.integer(period_raw - last_hist_year)),
      horizon_block = dplyr::case_when(
        !is.finite(horizon) ~ NA_character_,
        horizon <= 5L ~ "1_5",
        horizon <= 10L ~ "6_10",
        horizon <= 20L ~ "11_20",
        TRUE ~ "21p"
      )
    )

  if (is.data.frame(hist_df) && nrow(hist_df) && "cohort" %in% names(hist_df)) {
    hist_support <- hist_df %>%
      dplyr::mutate(cohort = suppressWarnings(as.integer(cohort))) %>%
      dplyr::filter(is.finite(cohort)) %>%
      dplyr::count(cohort, name = "support_n")

    max_support <- suppressWarnings(as.numeric(max(hist_support$support_n, na.rm = TRUE)))[1]
    if (!is.finite(max_support) || max_support <= 0) max_support <- 1

    out <- out %>%
      dplyr::left_join(hist_support %>% dplyr::rename(mapped_cohort_join = cohort),
                       by = c("mapped_cohort" = "mapped_cohort_join")) %>%
      dplyr::mutate(
        support_n = suppressWarnings(as.numeric(support_n)),
        support_frac = ifelse(is.finite(support_n) & max_support > 0, support_n / max_support, NA_real_)
      )
  } else {
    out <- out %>%
      dplyr::mutate(support_n = NA_real_, support_frac = NA_real_)
  }

  out %>%
    dplyr::select(dplyr::any_of(c(
      "sex", "age", "period", "cohort", "E", "exposure",
      "period_raw", "mapped_period", "period_is_clamped", "period_shift",
      "cohort_raw", "mapped_cohort", "cohort_is_edge", "cohort_clamp_low", "cohort_clamp_high", "cohort_shift",
      "support_n", "support_frac", "horizon", "horizon_block"
    )))
}

summarise_incidence_border_diag <- function(border_df, exposure_col = NULL) {
  if (!is.data.frame(border_df) || !nrow(border_df)) {
    return(tibble::tibble())
  }

  exposure_col <- exposure_col %||% intersect(c("exposure", "E"), names(border_df))[1]
  if (is.na(exposure_col) || !nzchar(exposure_col) || !(exposure_col %in% names(border_df))) {
    border_df$.w_border <- 1
  } else {
    border_df$.w_border <- suppressWarnings(as.numeric(border_df[[exposure_col]]))
    border_df$.w_border[!is.finite(border_df$.w_border) | border_df$.w_border <= 0] <- 1
  }

  .wmean_border <- function(x, w) {
    ok <- is.finite(x) & is.finite(w) & w > 0
    if (!any(ok)) return(NA_real_)
    stats::weighted.mean(x[ok], w[ok])
  }

  border_df %>%
    dplyr::mutate(
      sex = as.character(sex),
      horizon_block = factor(horizon_block, levels = c("1_5", "6_10", "11_20", "21p"))
    ) %>%
    dplyr::group_by(sex, horizon_block) %>%
    dplyr::summarise(
      n_cells = dplyr::n(),
      share_period_clamped = .wmean_border(as.numeric(period_is_clamped %in% TRUE), .w_border),
      share_cohort_edge = .wmean_border(as.numeric(cohort_is_edge %in% TRUE), .w_border),
      share_cohort_low = .wmean_border(as.numeric(cohort_clamp_low %in% TRUE), .w_border),
      share_cohort_high = .wmean_border(as.numeric(cohort_clamp_high %in% TRUE), .w_border),
      mean_abs_period_shift = .wmean_border(abs(suppressWarnings(as.numeric(period_shift))), .w_border),
      mean_abs_cohort_shift = .wmean_border(abs(suppressWarnings(as.numeric(cohort_shift))), .w_border),
      mean_support_frac = .wmean_border(suppressWarnings(as.numeric(support_frac)), .w_border),
      .groups = "drop"
    ) %>%
    dplyr::arrange(sex, horizon_block)
}

# --- ORTOGONALIDAD RW2(period) vs tendencia ----
make_slope_constr <- function(levels_vec, center) {
  t <- levels_vec - center
  list(A = matrix(t, nrow = 1), e = 0)   # impone sum(t_c * gamma_t) = 0
}


normalize_cause_id <- function(x) {
  k <- norm_txt(x)
  if (grepl("\bpulmon\b|\blung\b", k)) return("lung")
  if (grepl("\bpancreas\b", k)) return("pancreas")
  if (grepl("\brinon\b|\bkidney\b", k)) return("kidney")
  if (grepl("\bvejiga\b|\bbladder\b", k)) return("bladder")
  if (grepl("\blaringe\b|\blarynx\b", k)) return("larynx")
  if (grepl("\bestomago\b|\bstomach\b", k)) return("stomach")
  if (grepl("\besofago\b|\besophagus\b", k)) return("esophagus")
  if (grepl("cuello.*utero|cervix", k)) return("cervix")
  if (grepl("cavidad.*oral|faringe|oral.*phar|mouth|phary", k)) return("oralphar")
  return(k)
}

annualize_postdx_kernel <- function(p_0_1, p_1_3, p_3_5,
                                    mode = MORT_POSTDX_KERNEL_MODE) {
  p_0_1 <- as.numeric(p_0_1)[1]
  p_1_3 <- as.numeric(p_1_3)[1]
  p_3_5 <- as.numeric(p_3_5)[1]

  if (!all(is.finite(c(p_0_1, p_1_3, p_3_5)))) {
    stop("annualize_postdx_kernel(): probabilidades no finitas", call. = FALSE)
  }

  if (identical(mode, "midyear_uniform")) {
    w <- c(
      0.5 * p_0_1,
      0.5 * p_0_1 + 0.25 * p_1_3,
      0.5 * p_1_3,
      0.25 * p_1_3 + 0.25 * p_3_5,
      0.5 * p_3_5,
      0.25 * p_3_5
    )
  } else {
    stop("Modo de anualización post-diagnóstico no soportado: ", mode, call. = FALSE)
  }

  tibble::tibble(
    lag = 0:5,
    weight = pmax(as.numeric(w), 0),
    cum_weight = cumsum(pmax(as.numeric(w), 0))
  )
}

get_postdx_kernel <- function(cause_id, sex,
                              death_tbl = MORT_POSTDX_DEATH_TABLE,
                              mode = MORT_POSTDX_KERNEL_MODE) {
  cause_id <- normalize_cause_id(cause_id)
  sex <- toupper(substr(as.character(sex)[1], 1, 1))

  cause_id_key <- cause_id
  sex_key <- sex

  row <- death_tbl %>%
    dplyr::filter(.data$cause_id == .env$cause_id_key, .data$sex == .env$sex_key)

  if (nrow(row) != 1L) {
    stop(sprintf("No encuentro fila única en MORT_POSTDX_DEATH_TABLE para cause_id='%s', sex='%s'.", cause_id, sex),
         call. = FALSE)
  }

  annualize_postdx_kernel(
    p_0_1 = row$p_0_1[[1]],
    p_1_3 = row$p_1_3[[1]],
    p_3_5 = row$p_3_5[[1]],
    mode = mode
  ) %>%
    dplyr::mutate(
      cause_id = cause_id,
      sex = sex,
      p_0_1 = row$p_0_1[[1]],
      p_1_3 = row$p_1_3[[1]],
      p_3_5 = row$p_3_5[[1]],
      p_le_5 = row$p_le_5[[1]],
      kernel_mode = mode,
      .before = 1
    )
}

attach_external_mortality_offset <- function(mort_all, inc_rates_all, pop_all_tbl,
                                             cause_id, sex_sel,
                                             age_shift = MORT_POSTDX_USE_AGE_SHIFT,
                                             clip_min_age = MORT_POSTDX_CLIP_MIN_AGE,
                                             clip_min_year = MORT_POSTDX_CLIP_MIN_YEAR) {
  stopifnot(all(c("sex", "age", "period", "exposure") %in% names(mort_all)))
  stopifnot(all(c("sex", "age", "period", "rate_hat") %in% names(inc_rates_all)))
  stopifnot(all(c("sex", "age", "period", "exposure") %in% names(pop_all_tbl)))

  cause_id <- normalize_cause_id(cause_id)
  sex_sel <- toupper(substr(as.character(sex_sel)[1], 1, 1))
  kernel_tbl <- get_postdx_kernel(cause_id = cause_id, sex = sex_sel)

  inc_cases_all <- inc_rates_all %>%
    dplyr::left_join(
      pop_all_tbl %>% dplyr::select(.data$sex, .data$age, .data$period, .data$exposure),
      by = c("sex", "age", "period")
    ) %>%
    dplyr::mutate(
      exposure_inc = pmax(as.numeric(.data$exposure), 1e-12),
      cases_hat = pmax(as.numeric(.data$rate_hat), 1e-12) * exposure_inc
    ) %>%
    dplyr::select(.data$sex, .data$age, .data$period, .data$cases_hat)

  inc_cases_sel <- inc_cases_all %>% dplyr::filter(.data$sex == sex_sel)
  if (!nrow(inc_cases_sel)) {
    stop("attach_external_mortality_offset(): no hay incidencia disponible para el sexo seleccionado.", call. = FALSE)
  }

  age_min_inc <- suppressWarnings(min(inc_cases_sel$age, na.rm = TRUE))
  period_min_inc <- suppressWarnings(min(inc_cases_sel$period, na.rm = TRUE))

  base <- mort_all %>%
    dplyr::mutate(row_id = dplyr::row_number())

  contrib_list <- lapply(seq_len(nrow(kernel_tbl)), function(i) {
    lag_i <- as.integer(kernel_tbl$lag[[i]])
    w_i <- as.numeric(kernel_tbl$weight[[i]])

    age_diag_raw <- if (isTRUE(age_shift)) base$age - lag_i else base$age
    period_diag_raw <- base$period - lag_i

    age_diag <- if (isTRUE(clip_min_age)) pmax(age_diag_raw, age_min_inc) else age_diag_raw
    period_diag <- if (isTRUE(clip_min_year)) pmax(period_diag_raw, period_min_inc) else period_diag_raw

    base %>%
      dplyr::transmute(
        row_id = .data$row_id,
        sex = .data$sex,
        age_diag = age_diag,
        period_diag = period_diag,
        age_diag_raw = age_diag_raw,
        period_diag_raw = period_diag_raw
      ) %>%
      dplyr::left_join(
        inc_cases_sel,
        by = c("sex" = "sex", "age_diag" = "age", "period_diag" = "period")
      ) %>%
      dplyr::mutate(
        lag = lag_i,
        weight = w_i,
        cases_hat = dplyr::coalesce(.data$cases_hat, 0),
        deaths_ext_component = w_i * .data$cases_hat,
        clipped_age = .data$age_diag_raw != .data$age_diag,
        clipped_period = .data$period_diag_raw != .data$period_diag
      ) %>%
      dplyr::select(.data$row_id, .data$lag, .data$weight, .data$deaths_ext_component,
                    .data$clipped_age, .data$clipped_period)
  })

  contrib_tbl <- dplyr::bind_rows(contrib_list)

  contrib_summary <- contrib_tbl %>%
    dplyr::group_by(.data$lag) %>%
    dplyr::summarise(
      weight = dplyr::first(.data$weight),
      ext_deaths_total = sum(.data$deaths_ext_component, na.rm = TRUE),
      n_age_clipped = sum(.data$clipped_age, na.rm = TRUE),
      n_period_clipped = sum(.data$clipped_period, na.rm = TRUE),
      .groups = "drop"
    )

  offset_tbl <- contrib_tbl %>%
    dplyr::group_by(.data$row_id) %>%
    dplyr::summarise(
      mort_ext_deaths = sum(.data$deaths_ext_component, na.rm = TRUE),
      .groups = "drop"
    )


  out <- base %>%
    dplyr::left_join(offset_tbl, by = "row_id") %>%
    dplyr::mutate(
      E = pmax(as.numeric(.data$exposure), 1e-12),
      logE = log(.data$E),
      mort_ext_deaths = pmax(dplyr::coalesce(.data$mort_ext_deaths, 0), 1e-12),
      mort_ext_rate = .data$mort_ext_deaths / .data$E,
      log_mort_ext = log(.data$mort_ext_deaths),
      log_mort_ext_rate = log(pmax(.data$mort_ext_rate, 1e-12))
    ) %>%
    dplyr::select(-.data$row_id)

  attr(out, "kernel") <- kernel_tbl
  attr(out, "kernel_summary") <- contrib_summary
  out
}

# ---- Construir logI_lag a partir de cualquier tabla de tasas de Incidencia ----
attach_logI_to_mort <- function(mort_all, inc_rates_all,
                                L_I, Da_I, bridge_years, last_hist_year_m) {
  # mort_all: data.frame con sex, age, period, exposure, etc.
  # inc_rates_all: sex, age, period, rate_hat (hist+fut si aplica)
  stopifnot(all(c("sex","age","period") %in% names(mort_all)))
  stopifnot(all(c("sex","age","period","rate_hat") %in% names(inc_rates_all)))
  
  # 1) variables rezagadas para alinear Inc -> Mort
  df <- mort_all %>%
    dplyr::mutate(
      ageL = age - Da_I,
      perL = period - L_I
    )
  
  # 2) mapa de incidencia por (sex, ageL, perL)
  inc_map <- inc_rates_all %>%
    dplyr::transmute(
      sex,
      ageL = age - Da_I,
      perL = period,
      rate_I = pmax(rate_hat, 1e-12)
    )
  
  df <- df %>% dplyr::left_join(inc_map, by = c("sex","ageL","perL"))
  
  # 3) tasa "congelada" al último año histórico de incidencia (para el bridge)
  #    si no coincide exactamente con el último año histórico de mortalidad,
  #    tomamos el mayor período de incidencia <= last_hist_year_m; si no hay, tomamos el máximo global.
  last_hist_inc <- tryCatch(
    max(inc_rates_all$period[inc_rates_all$period <= last_hist_year_m], na.rm = TRUE),
    error = function(e) max(inc_rates_all$period, na.rm = TRUE)
  )
  
  inc_last <- inc_rates_all %>%
    dplyr::filter(period == last_hist_inc) %>%
    dplyr::transmute(
      sex,
      ageL = age - Da_I,
      rate_I0 = pmax(rate_hat, 1e-12)
    )
  
  df <- df %>% dplyr::left_join(inc_last, by = c("sex","ageL"))
  
  # 4) bridge: mezcla entre la tasa congelada (histórica) y la tasa rezagada (si existe)
  if (isTRUE(bridge_years > 0)) {
    df <- df %>%
      dplyr::mutate(
        w = pmax(0, pmin(1, (period - last_hist_year_m) / bridge_years)),
        rate_use = dplyr::case_when(
          period <= last_hist_year_m ~ dplyr::coalesce(rate_I, rate_I0),
          TRUE ~ (1 - w) * dplyr::coalesce(rate_I0, rate_I) + w * dplyr::coalesce(rate_I, rate_I0)
        ),
        logI_lag = log(pmax(rate_use, 1e-12))
      ) %>%
      dplyr::select(-w, -rate_use)
  } else {
    df <- df %>%
      dplyr::mutate(
        rate_use = dplyr::coalesce(rate_I, rate_I0),
        logI_lag = log(pmax(rate_use, 1e-12))
      ) %>%
      dplyr::select(-rate_use)
  }
  
  df$logI_lag[!is.finite(df$logI_lag)] <- log(1e-12)
  attr(df, "L_I") <- L_I
  return(df)
}

# Extrae efecto de período histórico (centrado) del fit APC de incidencia
extract_period_hist <- function(fit, lev_period, last_hist_year) {
  sp  <- fit$summary.random$period_id
  yrs <- lev_period[sp$ID]
  idx <- yrs <= last_hist_year
  tibble::tibble(period = yrs[idx],
                 per_eff = sp$mean[idx] - mean(sp$mean[idx]))
}


.rolling_mean_full <- function(x, k = 3L) {
  x <- as.numeric(x)
  k <- max(1L, as.integer(k))
  if (!length(x)) return(x)
  if (k <= 1L) return(x)
  if ((k %% 2L) == 0L) k <- k + 1L
  if (length(x) < k) return(rep(NA_real_, length(x)))

  h <- as.integer(k %/% 2L)
  out <- rep(NA_real_, length(x))
  for (i in seq_along(x)) {
    idx <- seq.int(i - h, i + h)
    if (all(idx >= 1L & idx <= length(x)) && all(is.finite(x[idx]))) {
      out[i] <- mean(x[idx])
    }
  }
  out
}


.expose_selected <- function(sex, A_I_used = NA, L_I_used = NA, bridge_years = NA) {
  if (!isTRUE(EXPOSE_SELECTED_LAGS_TO_ENV)) return(invisible())
  sx <- ifelse(identical(sex, "F"), "F", "M")
  if (!is.na(A_I_used))      assign(paste0("A_I_SELECTED_", sx), A_I_used, envir = .GlobalEnv)
  if (!is.na(L_I_used))      assign(paste0("L_I_SELECTED_", sx), L_I_used, envir = .GlobalEnv)
  if (!is.na(bridge_years))  assign(paste0("BRIDGE_SELECTED_", sx), bridge_years, envir = .GlobalEnv)
}

extract_beta_I <- function(fit, link_mode = MORT_I_LINK_MODE) {
  # En el nuevo camino externo no se estima beta: el enlace ya viene impuesto por el kernel.
  if (identical(link_mode, "external_kernel")) return(c(mean = 1, lwr = 1, upr = 1))
  if (identical(link_mode, "offset")) return(c(mean = 1, lwr = 1, upr = 1))
  if (!"logI_lag" %in% rownames(fit$summary.fixed)) return(c(mean = NA, lwr = NA, upr = NA))
  v <- fit$summary.fixed["logI_lag", c("mean","0.025quant","0.975quant")]
  names(v) <- c("mean","lwr","upr")
  as.numeric(v)
}

# =============================================================
# Helpers: scores (WAIC/DIC/LCPO) y export LaTeX
# =============================================================

get_fit_scores <- function(fit) {
  # Devuelve WAIC/DIC/LCPO para un objeto inla (o NA si no están disponibles)
  if (is.null(fit) || !inherits(fit, "inla")) {
    return(data.frame(WAIC = NA_real_, DIC = NA_real_, LCPO = NA_real_))
  }
  
  waic <- if (!is.null(fit$waic) && !is.null(fit$waic$waic)) fit$waic$waic else NA_real_
  dic  <- if (!is.null(fit$dic)  && !is.null(fit$dic$dic))   fit$dic$dic   else NA_real_
  
  cpo  <- if (!is.null(fit$cpo) && !is.null(fit$cpo$cpo)) fit$cpo$cpo else NA_real_
  fail <- if (!is.null(fit$cpo) && !is.null(fit$cpo$failure)) fit$cpo$failure else rep(FALSE, length(cpo))
  
  ok <- is.finite(cpo) & !is.na(cpo) & !fail
  lcpo <- if (any(ok)) sum(log(pmax(cpo[ok], 1e-15)), na.rm = TRUE) else NA_real_
  
  data.frame(WAIC = waic, DIC = dic, LCPO = lcpo)
}

# =============================================================
# Helpers: ΔLCPO con SE (pointwise) + backtest temporal Poisson
# =============================================================

.get_logcpo_vec <- function(fit, eps = 1e-15) {
  if (is.null(fit) || !inherits(fit, "inla") || is.null(fit$cpo) || is.null(fit$cpo$cpo)) {
    return(NULL)
  }
  cpo  <- fit$cpo$cpo
  fail <- fit$cpo$failure %||% rep(FALSE, length(cpo))
  out  <- rep(NA_real_, length(cpo))
  ok <- !fail & is.finite(cpo) & (cpo > 0)
  out[ok] <- log(pmax(cpo[ok], eps))
  out
}

delta_lcpo_se <- function(fit_alt, fit_ref, eps = 1e-15) {
  la <- .get_logcpo_vec(fit_alt, eps = eps)
  lr <- .get_logcpo_vec(fit_ref, eps = eps)
  if (is.null(la) || is.null(lr)) return(list(delta = NA_real_, se = NA_real_, n = 0L))
  ok <- is.finite(la) & is.finite(lr)
  d  <- la[ok] - lr[ok]
  if (!length(d)) return(list(delta = NA_real_, se = NA_real_, n = 0L))
  delta <- sum(d)
  se    <- sqrt(length(d) * stats::var(d))   # estilo loo: se(elpd_diff)
  list(delta = delta, se = se, n = as.integer(length(d)))
}

.log_pred_from_marg_mu <- function(y, marg_mu, eps_w = 1e-300) {
  # marg_mu: matriz (x, density) para mu = E[Y] (INLA fitted.values)
  if (is.null(marg_mu) || nrow(marg_mu) < 2) return(NA_real_)
  x <- marg_mu[, 1]
  d <- pmax(marg_mu[, 2], 0)
  dx <- diff(x)
  dx <- c(dx, tail(dx, 1))
  w <- d * dx
  s <- sum(w)
  if (!is.finite(s) || s <= 0) return(NA_real_)
  w <- w / s
  
  lt <- stats::dpois(y, lambda = pmax(x, 1e-12), log = TRUE) + log(pmax(w, eps_w))
  m  <- max(lt)
  m + log(sum(exp(lt - m)))
}

backtest_inla_poisson <- function(formula, data, train_end,
                                  tag,
                                  y_col = "y", period_col = "period",
                                  control.family = NULL,
                                  control.fixed = NULL,
                                  control.inla = NULL) {
  stopifnot(is.data.frame(data))
  dat <- data
  y   <- dat[[y_col]]
  per <- dat[[period_col]]
  hold <- !is.na(y) & per > train_end
  
  if (!any(hold)) {
    return(list(train_end = train_end, n = 0L, BT_LPD = NA_real_, BT_RMSE = NA_real_))
  }
  
  dat[[y_col]][hold] <- NA
  
  fit_bt <- inla_tag(tag, formula = formula, family = "poisson", data = dat,
                     control.family    = control.family,
                     control.fixed     = control.fixed,
                     control.predictor = list(compute = TRUE),
                     control.compute   = list(config = FALSE),
                     control.inla      = control.inla)
  
  mu_mean <- fit_bt$summary.fitted.values$mean[hold]
  BT_RMSE <- sqrt(mean((y[hold] - mu_mean)^2, na.rm = TRUE))
  
  idx <- which(hold)
  lpd_i <- vapply(idx, function(i) .log_pred_from_marg_mu(y[i], fit_bt$marginals.fitted.values[[i]]),
                  numeric(1))
  # --- BT LPD robusto ---
  # (A) Plug-in LPD usando el posterior mean de mu (siempre disponible)
  y_hold  <- y[hold]
  mu_hold <- pmax(mu_mean[hold], 1e-12)  # guardrail numérico
  BT_LPD_plugin <- sum(stats::dpois(y_hold, lambda = mu_hold, log = TRUE), na.rm = TRUE)
  
  # (B) Posterior-predictive LPD si existen marginals (opcional)
  BT_LPD_pp <- NA_real_
  if (!is.null(fit_bt$marginals.fitted.values)) {
    lpd_i <- vapply(which(hold), function(i) .log_pred_from_marg_mu(
      y = y[i],
      marg = fit_bt$marginals.fitted.values[[i]],
      w = if (!is.null(w)) w[i] else 1
    ), numeric(1))
    
    if (!all(is.na(lpd_i))) BT_LPD_pp <- sum(lpd_i, na.rm = TRUE)
  }
  
  # Prioridad: si pp existe, úsalo; si no, caé al plug-in
  BT_LPD <- if (is.finite(BT_LPD_pp)) BT_LPD_pp else BT_LPD_plugin
  
  if (is.null(fit_bt$marginals.fitted.values)) {
    .bapc_verbose("Backtest: marginals.fitted.values is NULL -> using plugin LPD")
  }
  
  list(train_end = train_end, n = as.integer(sum(hold)), BT_LPD = BT_LPD, BT_RMSE = BT_RMSE)
}

score_row <- function(model_label, fit) {
  s <- get_fit_scores(fit)
  data.frame(
    model = model_label,
    WAIC  = s$WAIC,
    DIC   = s$DIC,
    LCPO  = s$LCPO,
    stringsAsFactors = FALSE
  )
}

write_fit_scores_tex <- function(df, file, digits = 1) {
  stopifnot(is.data.frame(df))
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  
  keep <- intersect(c("sex","model","WAIC","DIC","LCPO","dLCPO","se_dLCPO","BT_LPD","BT_RMSE"), names(df))
  df2  <- df[, keep, drop = FALSE]
  
  # --- Bloques y etiquetas cortas ---
  df2 <- df2 %>%
    dplyr::mutate(
      block = dplyr::case_when(
        grepl("^Prevalence", model) ~ "Prevalence",
        grepl("^Incidence",  model) ~ "Incidence",
        grepl("^Mortality",  model) ~ "Mortality",
        TRUE ~ "Other"
      ),
      model_short = dplyr::case_when(
        model == "Prevalence BAPC"                 ~ "BAPC",
        model == "Incidence benchmark (APC)"       ~ "Benchmark (APC)",
        model == "Incidence prevalence-informed"   ~ "Prevalence-informed",
        model == "Mortality benchmark (APC)"       ~ "Benchmark (APC)",
        model == "Mortality anchored on I|P"       ~ "Anchored on I|P",
        model == "Mortality anchored on I only"    ~ "Anchored on I only",
        TRUE ~ model
      ),
      block = factor(block, levels = c("Prevalence","Incidence","Mortality","Other")),
      sex   = factor(sex, levels = c("F","M"))
    ) %>%
    dplyr::arrange(sex, block, model_short)
  
  fmt_num <- function(x) ifelse(is.na(x), "", formatC(x, format = "f", digits = digits))
  num_cols <- intersect(c("WAIC","DIC","LCPO","dLCPO","se_dLCPO","BT_LPD","BT_RMSE"), names(df2))
  for (nm in num_cols) df2[[nm]] <- fmt_num(df2[[nm]])
  
  # --- Cabecera (sin repetir "Prevalence/Incidence/Mortality" en cada fila) ---
  head_names <- c("Model","WAIC","DIC","LCPO","$\\Delta$LCPO","SE($\\Delta$LCPO)","BT LPD","BT RMSE")
  
  lines <- c(
    sprintf("%% Auto-generated by BAPC_PIM_8_loop.R on %s", Sys.Date()),
    "\\begin{tabular}{@{}lrrrrrrr@{}}",
    "\\toprule",
    paste(head_names, collapse = " & "), "\\\\",
    "\\midrule"
  )
  
  for (sx in levels(df2$sex)) {
    subx <- df2[df2$sex == sx, , drop = FALSE]
    if (nrow(subx) == 0) next
    
    # Sex header row (saves width vs a Sex column)
    sex_title <- if (sx == "F") "Females" else "Males"
    lines <- c(lines, sprintf("\\multicolumn{8}{@{}l}{\\textbf{%s}} \\\\", sex_title))
    lines <- c(lines, "\\addlinespace[0.25em]")
    
    for (blk in levels(subx$block)) {
      subb <- subx[subx$block == blk, , drop = FALSE]
      if (nrow(subb) == 0) next
      
      lines <- c(lines, sprintf("\\multicolumn{8}{@{}l}{\\textit{%s}} \\\\", blk))
      
      for (i in seq_len(nrow(subb))) {
        r <- subb[i, ]
        row <- c(
          r$model_short, r$WAIC, r$DIC, r$LCPO,
          if ("dLCPO" %in% names(r)) r$dLCPO else "",
          if ("se_dLCPO" %in% names(r)) r$se_dLCPO else "",
          if ("BT_LPD" %in% names(r)) r$BT_LPD else "",
          if ("BT_RMSE" %in% names(r)) r$BT_RMSE else ""
        )
        lines <- c(lines, paste(row, collapse = " & "), "\\\\")
      }
      
      lines <- c(lines, "\\addlinespace[0.35em]")
    }
    
    lines <- c(lines, "\\midrule")
  }
  
  if (tail(lines, 1) == "\\midrule") lines <- head(lines, -1)
  
  lines <- c(lines, "\\bottomrule", "\\end{tabular}")
  writeLines(lines, file, useBytes = TRUE)
}




# ---------- ARMAR EJES REALES DESDE LOS PARÁMETROS DEVUELTOS POR EL PIPELINE ----------
apc_axes_from_params <- function(res_sex, which = c("M","IP","I")) {
  which <- match.arg(which)
  p <- res_sex$params %||% list()
  
  if (which == "M") {
    a_min <- p$age_min_m %||% get0("AGE_M_MIN", ifnotfound = 35, envir = .GlobalEnv)
    a_max <- p$age_max_m %||% get0("AGE_M_MAX", ifnotfound = 89, envir = .GlobalEnv)
  } else {
    a_min <- p$age_min_i %||% get0("AGE_I_MIN", ifnotfound = 20, envir = .GlobalEnv)
    a_max <- p$age_max_i %||% get0("AGE_I_MAX", ifnotfound = 89, envir = .GlobalEnv)
  }
  t_min <- p$period_min_m %||% get0("PERIOD_M_MIN", ifnotfound = 1996, envir = .GlobalEnv)
  t_max <- p$period_max_m %||% get0("PERIOD_M_MAX", ifnotfound = 2022, envir = .GlobalEnv)
  
  age_vec    <- seq.int(as.integer(a_min), as.integer(a_max))
  period_vec <- seq.int(as.integer(t_min), as.integer(t_max))
  cohort_vec <- seq.int(min(period_vec) - max(age_vec), max(period_vec) - min(age_vec))
  
  list(age = age_vec, period = period_vec, cohort = cohort_vec)
}


.align_axis <- function(x_vals, n) {
  x_vals <- as.numeric(x_vals)
  if (length(x_vals) == n) return(x_vals)
  if (length(x_vals) >  n) return(tail(x_vals, n))  # si sobran, me quedo con el final
  if (length(x_vals) <  n) return(seq_len(n))       # fallback seguro
}

# Requiere patchwork
if (!"package:patchwork" %in% search()) {
  try({ library(patchwork) }, silent = TRUE)
}



## ===== Helpers diagnóstico PREV (tabla + export) =====

## ===== Helpers diagnóstico PREV =====
diag_prev_from_resboth <- function(rb) {
  suppressPackageStartupMessages({
    requireNamespace("dplyr"); requireNamespace("purrr"); requireNamespace("tibble")
  })
  src <- list(
    M = purrr::pluck(rb, "resM", "diag", "prev", .default = NULL),
    F = purrr::pluck(rb, "resF", "diag", "prev", .default = NULL)
  ) |> purrr::compact()
  if (length(src) == 0) return(tibble::tibble())
  src <- purrr::map(src, ~ tryCatch(tibble::as_tibble(.x), error = function(e) tibble::tibble()))
  out <- dplyr::bind_rows(src, .id = "sex")
  if (!"sex" %in% names(out)) out <- dplyr::mutate(out, sex = names(src)[1], .before = 1)
  out
}

diag_prev_from_scn <- function(scn) {
  purrr::imap_dfr(scn, function(x, nm) {
    tbl <- diag_prev_from_resboth(x)
    if (nrow(tbl) == 0) return(tibble::tibble())
    dplyr::mutate(tbl, scenario = nm, .before = 1)
  })
}

emit_diag_prev <- function(obj, out_dir = NULL, cause_label = NULL, scenario_label = NULL, also_print = TRUE, write_csv = TRUE) {
  # deps sin ensuciar el search path
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Falta 'dplyr'")
  if (!requireNamespace("readr", quietly = TRUE)) stop("Falta 'readr'")
  
  # carpeta de salida
  if (is.null(out_dir) || !length(out_dir) || !nzchar(out_dir)) {
    out_dir <- get0("BASE_RESULTS_DIR", ifnotfound = getwd(), inherits = TRUE)
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # construir tabla (acepta objeto 'res_both' o un elemento de 'scn')
  tbl <- if (!is.null(obj$resM) || !is.null(obj$resF)) {
    sc <- if (!is.null(scenario_label) && nzchar(as.character(scenario_label)[1])) as.character(scenario_label)[1] else "scenario"
    dplyr::mutate(diag_prev_from_resboth(obj), scenario = sc, .before = 1)
  } else {
    diag_prev_from_scn(obj)
  }
  
  # nada que reportar
  if (is.null(tbl) || !NROW(tbl)) {
    message(">> PREV diagnostic: no rows available (prev = NULL in both sexes).")
    return(invisible(NULL))
  }
  
  # etiqueta de causa
  if (is.null(cause_label) || !nzchar(cause_label)) {
    cause_label <- get0("CAUSE_ID", ifnotfound = NA_character_, inherits = TRUE)
    if (is.na(cause_label) || !nzchar(cause_label)) {
      cause_label <- if ("cause" %in% names(tbl)) {
        paste(unique(stats::na.omit(tbl$cause)), collapse = "_")
      } else "unknown_cause"
    }
  }
  
  # print compacto en consola (si se pide)
  if (isTRUE(also_print)) {
    if (!"sex" %in% names(tbl)) {
      tbl <- dplyr::mutate(tbl, sex = NA_character_, .before = 1)
    }
    proj_to <- get0("PROJ_TO", ifnotfound = NA_integer_, inherits = TRUE)
    zcol <- if (!is.na(proj_to)) paste0("z_", proj_to) else NA_character_
    if (!(zcol %in% names(tbl))) {
      zcand <- grep("^z_\\d+$", names(tbl), value = TRUE)
      zcol <- if (length(zcand)) tail(sort(zcand), 1L) else NULL
    }
    cols <- c("scenario","sex","cause","last_hist_year",
              "beta_mode","rr_inc","s_histP","z_base")
    print(dplyr::select(tbl, dplyr::any_of(cols)), n = Inf)
  }
  
  # guardar CSV (si se pide)
  if (isTRUE(write_csv)) {
    out_file <- file.path(out_dir, sprintf("diag_prev_%s.csv", cause_label))
    readr::write_csv(tbl, out_file)
    message(">> diag_prev saved to: ", out_file)
  }

  invisible(tbl)
}


# Calibración de QUIT_FLOOR_SD a partir de FAP (sin exfumadores)

.calib_floor_z_one_sex <- function(beta_pos, z_freeze_end,
                                   FAP_freeze_end,
                                   quit_mode = QUIT_MODE,
                                   HL        = QUIT_HALF_LIFE,
                                   rampY     = QUIT_RAMP_YEARS,
                                   proj_to  = PROJ_TO,
                                   base_year= PREV_BASE_YEAR) {
  if (!is.finite(beta_pos) || beta_pos <= 1e-8 || !is.finite(z_freeze_end)) return(NA_real_)
  Rtarget <- 1 - FAP_freeze_end
  Rtarget <- max(min(Rtarget, 1 - 1e-6), 1e-6)   # guardarraíl numérico
  z_quit_end <- z_freeze_end + (1/beta_pos) * log(Rtarget)
  
  if (identical(quit_mode, "floor")) {
    zf <- z_quit_end
  } else if (identical(quit_mode, "decay")) {
    t   <- proj_to - base_year
    lam <- log(2) / HL
    a   <- exp(-lam * t)
    zf  <- (z_quit_end - a * z_freeze_end) / (1 - a)
  } else { # "ramp"
    t <- proj_to - base_year
    w <- if (rampY > 0) min(1, t / rampY) else 1
    zf <- (z_quit_end - (1 - w) * z_freeze_end) / w
  }
  zf
}

# diag_freeze_tbl: salida de diag_prev_from_scn(...) filtrada a scenario=="freeze"
# Devuelve c(M=..., F=...) en **SD** del índice (lo que espera QUIT_FLOOR_SD)
calibrate_quit_floor_sd_by_sex <- function(diag_freeze_tbl, FAP_M, FAP_F, wI = W_I) {
  if (nrow(diag_freeze_tbl) == 0) return(c(M=NA_real_, F=NA_real_))
  zcol <- paste0("z_", PROJ_TO)
  if (!zcol %in% names(diag_freeze_tbl) && "z_proj_end" %in% names(diag_freeze_tbl)) zcol <- "z_proj_end"
  
  get_row <- function(sex_code) {
    dplyr::filter(diag_freeze_tbl, .data$sex == sex_code) %>% dplyr::slice(1)
  }
  rowM <- get_row("M"); rowF <- get_row("F")
  
  beta_pos_M <- if ("beta_P_pos" %in% names(rowM) && is.finite(rowM$beta_P_pos[1])) rowM$beta_P_pos[1] else 1
  beta_pos_F <- if ("beta_P_pos" %in% names(rowF) && is.finite(rowF$beta_P_pos[1])) rowF$beta_P_pos[1] else 1
  zM <- .calib_floor_z_one_sex(beta_pos = beta_pos_M, z_freeze_end = rowM[[zcol]][1], FAP_freeze_end = FAP_M)
  zF <- .calib_floor_z_one_sex(beta_pos = beta_pos_F, z_freeze_end = rowF[[zcol]][1], FAP_freeze_end = FAP_F)
  
  c(M = zM / wI, F = zF / wI)
}

# === DEBUG + WRAPPER LIMPIO ==================================
ok_modes <- c("logical","integer","double","numeric","complex","character","raw","list","expression")

.check_no_mode_in_f <- function(form) {
  txt <- paste(deparse(form), collapse = " ")
  if (grepl("f\\s*\\([^\\)]*\\bmode\\s*=", txt)) {
    stop("Se encontró 'mode=' dentro de f() (debe ser 'model='). Fórmula: ", txt, call. = FALSE)
  }
}

.inla_dump_on_error <- function(tag, e) {
  cat(sprintf("\n[INLA-FAIL:%s] %s\n", tag, conditionMessage(e)))
  dump.frames("dump_inla_err", TRUE)  # para traceback si querés
  stop(sprintf("[%s] %s", tag, conditionMessage(e)), call. = FALSE)
}

.edge_distance_to_border <- function(x) {
  lev <- sort(unique(as.numeric(x)))
  idx <- match(as.numeric(x), lev)
  pmin(idx - 1L, length(lev) - idx)
}

build_edge_weights_hist <- function(df,
                                    age_col = "age",
                                    period_col = "period",
                                    k_age = EDGE_WEIGHT_K_AGE,
                                    k_period = EDGE_WEIGHT_K_PERIOD,
                                    strength = EDGE_WEIGHT_STRENGTH,
                                    weight_min = EDGE_WEIGHT_MIN,
                                    geometry = EDGE_WEIGHT_GEOMETRY) {
  if (!is.data.frame(df) || nrow(df) == 0) return(tibble::tibble())
  
  d_age_edge <- .edge_distance_to_border(df[[age_col]])
  d_period_edge <- .edge_distance_to_border(df[[period_col]])
  
  k_age <- max(as.numeric(k_age), 1)
  k_period <- max(as.numeric(k_period), 1)
  
  e_age <- pmax(0, 1 - d_age_edge / k_age)
  e_period <- pmax(0, 1 - d_period_edge / k_period)
  
  edge_score <- switch(
    geometry,
    additive_mean = (e_age + e_period) / 2,
    stop("EDGE_WEIGHT_GEOMETRY no soportada: ", geometry)
  )
  
  edge_weight <- pmax(as.numeric(weight_min), 1 - as.numeric(strength) * edge_score)
  
  tibble::tibble(
    d_age_edge = as.integer(d_age_edge),
    d_period_edge = as.integer(d_period_edge),
    e_age = as.numeric(e_age),
    e_period = as.numeric(e_period),
    edge_score = as.numeric(edge_score),
    edge_weight = as.numeric(edge_weight)
  )
}

attach_edge_weights_hist <- function(df,
                                     stage,
                                     age_col = "age",
                                     period_col = "period",
                                     enabled = EDGE_WEIGHTING_ON,
                                     k_age = EDGE_WEIGHT_K_AGE,
                                     k_period = EDGE_WEIGHT_K_PERIOD,
                                     strength = EDGE_WEIGHT_STRENGTH,
                                     weight_min = EDGE_WEIGHT_MIN,
                                     geometry = EDGE_WEIGHT_GEOMETRY) {
  if (!is.data.frame(df) || nrow(df) == 0) return(df)
  
  if (!isTRUE(enabled)) {
    return(
      df %>%
        dplyr::mutate(
          edge_stage = stage,
          edge_geometry = geometry,
          d_age_edge = NA_integer_,
          d_period_edge = NA_integer_,
          e_age = 0,
          e_period = 0,
          edge_score = 0,
          edge_weight = 1
        )
    )
  }
  
  ew <- build_edge_weights_hist(
    df = df,
    age_col = age_col,
    period_col = period_col,
    k_age = k_age,
    k_period = k_period,
    strength = strength,
    weight_min = weight_min,
    geometry = geometry
  )
  
  dplyr::bind_cols(
    df,
    tibble::tibble(
      edge_stage = stage,
      edge_geometry = geometry
    ),
    ew
  )
}

summarise_edge_weights_hist <- function(df, stage = NA_character_) {
  if (!is.data.frame(df) || nrow(df) == 0 || !("edge_weight" %in% names(df))) {
    return(tibble::tibble())
  }
  
  out <- df %>%
    dplyr::group_by(dplyr::across(dplyr::any_of("sex"))) %>%
    dplyr::summarise(
      n_rows = dplyr::n(),
      min_weight = min(edge_weight, na.rm = TRUE),
      q10_weight = as.numeric(stats::quantile(edge_weight, 0.10, na.rm = TRUE, names = FALSE, type = 7)),
      mean_weight = mean(edge_weight, na.rm = TRUE),
      median_weight = stats::median(edge_weight, na.rm = TRUE),
      q90_weight = as.numeric(stats::quantile(edge_weight, 0.90, na.rm = TRUE, names = FALSE, type = 7)),
      max_weight = max(edge_weight, na.rm = TRUE),
      share_downweighted = mean(edge_weight < 0.999999, na.rm = TRUE),
      share_floor = mean(abs(edge_weight - min(edge_weight, na.rm = TRUE)) < 1e-12, na.rm = TRUE),
      mean_edge_score = mean(edge_score, na.rm = TRUE),
      .groups = "drop"
    )
  
  if (!is.na(stage)) out <- dplyr::mutate(out, stage = stage, .before = 1)
  out
}

inla_tag <- function(tag,
                     formula,
                     family,                 # <-- ahora explícito
                     data,
                     control.family = NULL,
                     control.fixed = NULL,
                     control.predictor = NULL,
                     control.compute = NULL,
                     control.inla = NULL,
                     ...) {                  # <-- ... al final (reenvío limpio)
  .check_no_mode_in_f(formula)
  
  dots_extra <- list(...)
  if (!is.null(dots_extra$weights)) {
    old_enable_weights <- tryCatch(INLA::inla.getOption("enable.inla.argument.weights"), error = function(e) NULL)
    tryCatch(INLA::inla.setOption(enable.inla.argument.weights = TRUE), error = function(e) NULL)
    
    if (!is.null(old_enable_weights)) {
      on.exit(
        tryCatch(INLA::inla.setOption(enable.inla.argument.weights = old_enable_weights), error = function(e) NULL),
        add = TRUE
      )
    }
  }
  
  # INLA ya no acepta `scale` dentro de `control.family` para Poisson en esta rama.
  # Conservamos la validación mínima y lo removemos antes del ajuste para evitar warnings
  # ruidosos y mantener la semántica del pipeline estable.
  if (!is.null(control.family) && !is.null(control.family$scale)) {
    sc <- control.family$scale
    if (!is.numeric(sc)) stop(sprintf("[%s] control.family$scale debe ser numérico", tag), call. = FALSE)
    if (NROW(sc) != NROW(data)) {
      stop(sprintf("[%s] largo de scale=%d ≠ nrow(data)=%d", tag, length(sc), NROW(data)), call. = FALSE)
    }
    control.family$scale <- NULL
    if (!length(control.family)) control.family <- NULL
  }
  
  .bapc_verbose(sprintf("[INLA-START:%s]", tag))
  out <- tryCatch(
    INLA::inla(
      formula = formula,
      family  = family,          # <-- reenviado explícitamente
      data    = data,
      control.family    = control.family,
      control.fixed     = control.fixed,
      control.predictor = control.predictor,
      control.compute   = control.compute,
      control.inla      = control.inla,
      ...               # Ntrials, E, etc.
    ),
    error = function(e) .inla_dump_on_error(tag, e)
  )
  .bapc_verbose(sprintf("[INLA-OK:%s]", tag))
  out
}

.chk <- function(tag) .bapc_verbose(sprintf("[CHK:%s]", tag))

.chk_vec <- function(tag, x) {
  if (inherits(x, "integer64")) x <- as.numeric(x)  # por si viniera de readr/bit64
  if (!is.numeric(x)) stop(sprintf("[%s] vector no numérico (mode=%s, class=%s)",
                                   tag, mode(x), paste(class(x), collapse=",")), call. = FALSE)
  if (any(!is.finite(x))) stop(sprintf("[%s] vector con valores no finitos", tag), call. = FALSE)
  x
}
# =============================================================





# --- Sanitización de objetos INLA para evitar referencias circulares y archivos pesados ---
clean_res_inla <- function(obj) {
  if (is.null(obj)) return(NULL)
  if (inherits(obj, "inla")) {
    obj$.args <- NULL
    obj$all.hyper <- NULL
    # obj$misc$configs <- NULL # opcional; a veces se usa para muestreo
    # obj$misc$lincomb.derived.correlation.matrix <- NULL
    # obj$misc$lincomb.derived.covariance.matrix <- NULL
    if (!is.null(obj$formula)) environment(obj$formula) <- .GlobalEnv
    if (!is.null(obj$call)) {
       # Limpiar el call para que no guarde el dataset entero en la cadena de texto
       # obj$call <- NULL 
    }
  }
  obj
}

# --- Sanitización recursiva de listas de resultados ---
sanitize_pipeline_output <- function(out) {
  if (is.list(out)) {
    # Evitar recursión infinita si ya hay un ciclo (aunque intentamos prevenirlos)
    # Por ahora solo limpiamos los niveles conocidos
    if (!is.null(out$fit_prev)) out$fit_prev <- clean_res_inla(out$fit_prev)
    if (!is.null(out$inc_fit$fit_inc)) out$inc_fit$fit_inc <- clean_res_inla(out$inc_fit$fit_inc)
    if (!is.null(out$fit_bapc)) out$fit_bapc <- clean_res_inla(out$fit_bapc)
    if (!is.null(out$fit_anchor_cond)) out$fit_anchor_cond <- clean_res_inla(out$fit_anchor_cond)
    if (!is.null(out$fit_anchor_noP)) out$fit_anchor_noP <- clean_res_inla(out$fit_anchor_noP)
  }
  out
}
