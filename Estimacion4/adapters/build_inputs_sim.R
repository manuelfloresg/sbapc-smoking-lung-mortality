# =============================================================
# Simulation input adapter
# =============================================================
# This module generates simulated data and standardizes it to
# the same contract expected by the core pipeline.

.make_rw1 <- function(n, sd) {
  x <- cumsum(stats::rnorm(n, 0, sd))
  x - mean(x)
}

.make_lag_kernel <- function(Lmax = 7L, center = 2.5, sd = 1.2) {
  Ls <- 0:Lmax
  w <- stats::dnorm(Ls, mean = center, sd = sd)
  w / sum(w)
}

.sim_delta_z_period <- function(grid, p0_by_sex, base_year, rate = 0.01, unit = 1.0) {
  k <- pmax(0L, grid$period - base_year)
  p0 <- p0_by_sex[as.character(grid$sex)]
  p_t <- p0 * (1 - rate)^k
  dz <- unit * (p_t - p0)
  dz[grid$period <= base_year] <- 0
  dz
}

simulate_prev_micro <- function(prev_truth, n_cell = 200L) {
  stopifnot(all(c("sex","age","period","cohort","p_true") %in% names(prev_truth)))
  idx <- rep(seq_len(nrow(prev_truth)), each = as.integer(n_cell))
  df <- prev_truth[idx, c("sex","age","period","cohort","p_true")]
  tibble::tibble(
    period = as.integer(df$period),
    age    = as.integer(df$age),
    cohort = as.integer(df$cohort),
    sex    = factor(as.character(df$sex), levels = c("M", "F")),
    fuma   = stats::rbinom(nrow(df), size = 1L, prob = pmax(pmin(df$p_true, 1 - 1e-8), 1e-8)),
    w      = 1,
    d_act  = 1L,
    d_12m  = 0L,
    d_30d  = 0L
  )
}

simulate_PIM_data <- function(cause_id = "lung",
                              seed = 1,
                              sexes = c("M", "F"),
                              age_min = 35, age_max = 89,
                              per_min = 1990, per_max = PROJ_TO,
                              last_hist = 2022,
                              N_prev_micro = 200L,
                              exposure = 100000,
                              A_true = 16L,
                              Lmax = 7L, L_center = 2.5, L_sd = 1.2,
                              theta_I = 0.35,
                              beta_mode = BETA_MODE,
                              rr_I = NULL,
                              theta_M = 1.0,
                              dgp = c("spec_linear", "misspec_tanh", "weak"),
                              kappa = 0.8,
                              scen_rate = 0.01,
                              scenario_name = c("freeze","up1pc","down1pc","down3pc","quit"),
                              prev_annual_rate_down3 = PREV_ANNUAL_RATE_DOWN3,
                              quit_half_life = QUIT_HALF_LIFE,
                              base_year = 2022,
                              prev_sex_shift = c(M = 0.15, F = -0.15),
                              inc_sex_shift  = c(M = 0.10, F = -0.10),
                              mort_sex_shift = c(M = 0.05, F = -0.05),
                               ...) {
  # Preserve RNG state to avoid side effects
  old_rng <- RNGkind()
  on.exit(do.call(RNGkind, as.list(old_rng)))
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  
  set.seed(seed)
  dgp <- match.arg(dgp)
  scenario_name <- match.arg(as.character(scenario_name)[1], c("freeze","up1pc","down1pc","down3pc","quit"))
  beta_mode <- match.arg(as.character(beta_mode)[1], c("estimate","prior_ols","offset","fixed_rr_offset"))
  sexes <- intersect(c("M", "F"), unique(as.character(sexes)))
  if (!length(sexes)) stop("sexes debe incluir M y/o F")
  args <- list(...)
  if (is.null(rr_I)) {
    if (!is.null(args$rr_inc)) {
      rr_I <- args$rr_inc
    } else {
      rr_I <- if (exists("INC_RR_TABLE")) INC_RR_TABLE$lung else 4.0
    }
  }
  if (isTRUE(BAPC_VERBOSE)) {
    message(sprintf("simulate_PIM_data: rr_I=%.2f, theta_I=%.2f", 
                    if(is.numeric(rr_I)) mean(rr_I) else -1, theta_I))
  }

  # Extraer RR real de la configuración del proyecto para la causa especificada
  rr_scalar_default <- get_inc_rr_by_cause(cause_id)
  rr_by_sex <- stats::setNames(rep(rr_scalar_default, length(sexes)), sexes)
  for (sx in sexes) {
    rr_by_sex[sx] <- get_inc_rr_by_cause_sex(cause_id, sx)
  }

  grid <- tidyr::expand_grid(
    sex = factor(sexes, levels = c("M", "F")),
    age = age_min:age_max,
    period = per_min:per_max
  ) %>% dplyr::mutate(cohort = period - age, grid_id = dplyr::row_number())

  lev_age <- sort(unique(grid$age))
  lev_per <- sort(unique(grid$period))
  lev_coh <- sort(unique(grid$cohort))

  # --------- PREVALENCIA (DGP binomial con APC) ----------
  aP_raw <- .make_rw1(length(lev_age), 0.08)
  aP_raw[lev_age > 65] <- aP_raw[max(which(lev_age <= 65))]
  aP <- aP_raw; names(aP) <- lev_age

  message(sprintf(">>> DGP SEED %s | aP[1:3]: %s", seed, paste(round(head(aP, 3), 4), collapse=" ")))
  pP <- .make_rw1(length(lev_per), 0.06); names(pP) <- lev_per
  cP <- .make_rw1(length(lev_coh), 0.10); names(cP) <- lev_coh
  sex_shift_prev <- prev_sex_shift[as.character(grid$sex)]
  sex_shift_prev[is.na(sex_shift_prev)] <- 0

  grid_age_clamped <- pmin(grid$age, 65)
  
  etaP <- -2.0 + sex_shift_prev +
    aP[as.character(grid_age_clamped)] + pP[as.character(grid$period)] + cP[as.character(grid$cohort)]
  p_true <- plogis(etaP)

  prev_truth <- grid %>%
    dplyr::transmute(sex, age, period, cohort, p_true = p_true)
  prev_micro <- simulate_prev_micro(prev_truth = prev_truth, n_cell = N_prev_micro)

  gammaP_true <- tibble::tibble(
    cohort = lev_coh,
    gammaP = as.numeric(cP)
  )

  p0_base <- prev_truth %>%
    dplyr::filter(period == base_year) %>%
    dplyr::group_by(sex) %>%
    dplyr::summarise(p0 = mean(p_true, na.rm = TRUE), .groups = "drop") %>%
    tibble::deframe()

  # --------- z_prev verdadero (baseline) ----------
  z_base <- if (identical(beta_mode, "fixed_rr_offset")) {
    rr_vec <- suppressWarnings(as.numeric(rr_by_sex[as.character(grid$sex)]))
    rr_vec[!is.finite(rr_vec) | rr_vec <= 1] <- rr_scalar_default
    
    # IMPORTANTE: El DGP debe usar la misma lógica que build_prev_rr_offset_stock_for_inc
    grid %>%
      dplyr::mutate(p_true = p_true) %>%
      dplyr::arrange(sex, cohort, age) %>%
      dplyr::group_by(sex, cohort) %>%
      dplyr::group_modify(function(df, key) {
        sx <- as.character(key$sex)
        L_q <- get_prev_inc_quit_horizon(cause_id = cause_id, sex = sx)
        w_vec <- .prev_rr_schedule_vec(cause_id = cause_id, sex = sx, L_q = L_q)
        rr_use <- as.numeric(rr_by_sex[sx])
        
        df$q_eff <- .calculate_stock_former_q_eff(as.numeric(df$p_true), w_vec, L_q)
        df$z_prev <- log1p(df$q_eff * (rr_use - 1))
        df
      }) %>%
      dplyr::ungroup()
  } else {
    build_prev_index_for_inc(df_inc_grid = grid, gammaP_all = gammaP_true, A_I = A_true, w_I = 1)
  } %>%
    dplyr::select(dplyr::any_of(c("sex", "age", "period", "cohort", "grid_id", "q_eff", "z_prev"))) %>%
    dplyr::right_join(grid %>% dplyr::select(sex, age, period, cohort, grid_id),
                      by = c("sex", "age", "period", "cohort", "grid_id")) %>%
    dplyr::arrange(grid_id)
  if (anyNA(z_base$z_prev)) stop("simulate_PIM_data: z_base quedó con NA tras alinear por grid_id.")

  s_hist <- stats::sd(z_base$z_prev[z_base$period <= base_year], na.rm = TRUE)
  unit <- 1 / pmax(s_hist, 1e-9)
  message(sprintf(">>> simulate_PIM_data: beta_mode='%s' scenario='%s'", beta_mode, scenario_name))
  if (identical(beta_mode, "fixed_rr_offset")) {
    message(">>> simulate_PIM_data: entró en bloque fixed_rr_offset")
    # Empezamos con la prevalencia corriente 'verdadera' (pP/p_true)
    z_scen <- z_base %>% dplyr::mutate(p_curr = p_true)
    k <- pmax(0L, z_scen$period - base_year)
    p0v <- p0_base[as.character(z_scen$sex)]
    idx_fut <- z_scen$period > base_year

    if (identical(scenario_name, "freeze")) {
      # No hay cambio en p_curr (se mantiene la tendencia proyectada de p_true)
    } else if (identical(scenario_name, "up1pc")) {
      fac <- (1 + scen_rate)^k
      z_scen$p_curr[idx_fut] <- pmin(pmax(z_scen$p_curr[idx_fut] * fac[idx_fut], 1e-8), 1 - 1e-8)
    } else if (identical(scenario_name, "down1pc")) {
      fac <- (1 - scen_rate)^k
      z_scen$p_curr[idx_fut] <- pmin(pmax(z_scen$p_curr[idx_fut] * fac[idx_fut], 1e-8), 1 - 1e-8)
    } else if (identical(scenario_name, "quit")) {
      # POLÍTICA: Todo el mundo deja de fumar de golpe (Instant Drop)
      z_scen$p_curr[idx_fut] <- 0
    }
    
    # RECALCULAR q_eff y z_prev a partir de la nueva trayectoria de p_curr
    # Usamos la misma maquinaria estructural que el estimador
    z_scen <- z_scen %>%
      dplyr::arrange(sex, cohort, age) %>%
      dplyr::group_by(sex, cohort) %>%
      dplyr::group_modify(function(df, key) {
        sx <- as.character(key$sex)
        L_q <- get_prev_inc_quit_horizon(cause_id = cause_id, sex = sx)
        w_vec <- .prev_rr_schedule_vec(cause_id = cause_id, sex = sx, L_q = L_q)
        rr_use <- as.numeric(rr_by_sex[sx])
        
        # El motor de stock convierte p_curr en q_eff usando el kernel w_vec
        df$q_eff <- .calculate_stock_former_q_eff(as.numeric(df$p_curr), w_vec, L_q)
        df$z_prev <- log1p(df$q_eff * (rr_use - 1))
        df
      }) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(grid_id)

  } else {
    dz <- .sim_delta_z_period(z_base, p0_by_sex = p0_base, base_year = base_year, rate = scen_rate, unit = unit)
    z_scen <- z_base %>% dplyr::mutate(z_prev = z_prev + dz) %>% dplyr::arrange(grid_id)
  }

  # Preserve the exact z regressors used in the DGP incidence equation
  zI_base_true_used <- z_base %>% dplyr::transmute(sex, age, period, cohort, grid_id, z_prev_used = z_prev)
  zI_scen_true_used <- z_scen %>% dplyr::transmute(sex, age, period, cohort, grid_id, z_prev_used = z_prev)

  # --------- INCIDENCIA (DGP Poisson con APC + g(z)) ----------
  aI <- .make_rw1(length(lev_age), 0.02); names(aI) <- lev_age
  # Estabilizar período en proyección para que los escenarios sean visibles
  pI_raw <- .make_rw1(length(lev_per), 0.015)
  idx_fut_per <- lev_per > last_hist
  if (any(idx_fut_per)) {
    pI_raw[idx_fut_per] <- pI_raw[max(which(!idx_fut_per))] 
  }
  pI <- pI_raw; names(pI) <- lev_per
  cI <- .make_rw1(length(lev_coh), 0.03); names(cI) <- lev_coh
  sex_shift_inc <- inc_sex_shift[as.character(grid$sex)]
  sex_shift_inc[is.na(sex_shift_inc)] <- 0

  # Sincronizar intercepto: 
  # Si el estimador usa offset(z - mean_z), el intercepto estimado es log_rate(mean_z).
  # El DGP usa exp(-8.5 + sex_shift + z). 
  # Para que coincidan, el DGP debería generar casos coherentes con un nivel base estable.
  etaI_apc <- -8.5 + sex_shift_inc +
    aI[as.character(grid$age)] + pI[as.character(grid$period)] + cI[as.character(grid$cohort)]

  z_base_used <- grid %>%
    dplyr::select(sex, age, period, cohort, grid_id) %>%
    dplyr::left_join(z_base, by = c("sex", "age", "period", "cohort", "grid_id")) %>%
    dplyr::arrange(grid_id)
    
  # IMPORTANTE: Asegurar que z_prev en el DGP sea exactamente lo que el estimador usará
  # pero con la no-linealidad si dgp == "misspec_tanh"
  g_z <- function(z, dgp_type) {
    if (dgp_type == "misspec_tanh") {
      # Aplicamos una saturación tipo tanh para romper la linealidad
      # Escalamos z por su desvío histórico para que el efecto sea visible
      return(2.0 * tanh(z / 1.5)) 
    }
    return(z) # spec_linear (default)
  }

  etaI_z_base <- g_z(z_base_used$z_prev, dgp)
  etaI_base <- etaI_apc + etaI_z_base
  rateI_base <- exp(etaI_base)
  
  z_scen_used <- grid %>%
    dplyr::select(sex, age, period, cohort, grid_id) %>%
    dplyr::left_join(z_scen, by = c("sex", "age", "period", "cohort", "grid_id")) %>%
    dplyr::arrange(grid_id)
  
  etaI_scen <- etaI_apc + g_z(z_scen_used$z_prev, dgp)
  rateI_scen <- exp(etaI_scen)

  idx_hist <- grid$period <= last_hist
  lamI_hist <- exposure * rateI_base[idx_hist]
  cases_hist <- stats::rpois(sum(idx_hist), lambda = lamI_hist)

  inc_hist <- grid[idx_hist, ] %>%
    dplyr::mutate(cases = cases_hist) %>%
    dplyr::select(sex, age, period, cohort, cases)

  pop_all <- grid %>%
    dplyr::transmute(sex, age, period, exposure = exposure)

  # --------- MORTALIDAD (DGP CONSISTENTE: INCIDENCIA ⊗ KERNEL) ----------
  # Extraemos el kernel de referencia para metadatos (usando M como default)
  wL <- get_postdx_kernel(cause_id = cause_id, sex = sexes[1])$weight

  build_mort_from_inc <- function(sex_lab) {
    k_tbl <- get_postdx_kernel(cause_id = cause_id, sex = sex_lab)
    idx_sex <- which(as.character(grid$sex) == sex_lab)

    inc_sex <- grid[idx_sex, ] %>%
      dplyr::mutate(
        inc_cases_base_true = pmax(rateI_base[idx_sex], 0) * exposure,
        inc_cases_scen_true = pmax(rateI_scen[idx_sex], 0) * exposure
      ) %>%
      dplyr::select(sex, age, period, inc_cases_base_true, inc_cases_scen_true)

    base_grid <- grid[idx_sex, ] %>%
      dplyr::mutate(row_id = dplyr::row_number())

    contrib <- lapply(seq_len(nrow(k_tbl)), function(i) {
      lag_i <- k_tbl$lag[i]
      w_i <- k_tbl$weight[i]

      base_grid %>%
        dplyr::mutate(
          age_diag_raw = age - lag_i,
          period_diag_raw = period - lag_i,
          age_diag = pmax(age_diag_raw, age_min),
          period_diag = pmax(period_diag_raw, per_min)
        ) %>%
        dplyr::left_join(
          inc_sex,
          by = c("sex", "age_diag" = "age", "period_diag" = "period")
        ) %>%
        dplyr::mutate(
          inc_cases_base_true = dplyr::coalesce(inc_cases_base_true, 0),
          inc_cases_scen_true = dplyr::coalesce(inc_cases_scen_true, 0),
          deaths_base = w_i * inc_cases_base_true,
          deaths_scen = w_i * inc_cases_scen_true
        ) %>%
        dplyr::select(row_id, deaths_base, deaths_scen)
    }) %>% dplyr::bind_rows()

    deaths_sum <- contrib %>%
      dplyr::group_by(row_id) %>%
      dplyr::summarise(
        mort_deaths_base_true = sum(deaths_base, na.rm = TRUE),
        mort_deaths_scen_true = sum(deaths_scen, na.rm = TRUE),
        .groups = "drop"
      )

    base_grid %>%
      dplyr::left_join(deaths_sum, by = "row_id") %>%
      dplyr::mutate(
        mort_deaths_base_true = pmax(dplyr::coalesce(mort_deaths_base_true, 0), 1e-12),
        mort_deaths_scen_true = pmax(dplyr::coalesce(mort_deaths_scen_true, 0), 1e-12),
        rateM_base_true = mort_deaths_base_true / exposure,
        rateM_scen_true = mort_deaths_scen_true / exposure,
        etaM_base_true = log(pmax(rateM_base_true, 1e-12)),
        etaM_scen_true = log(pmax(rateM_scen_true, 1e-12))
      ) %>%
      dplyr::select(
        sex, age, period, cohort, grid_id,
        mort_deaths_base_true, mort_deaths_scen_true,
        rateM_base_true, rateM_scen_true,
        etaM_base_true, etaM_scen_true
      )
  }

  mort_all <- dplyr::bind_rows(
    build_mort_from_inc("M"),
    build_mort_from_inc("F")
  ) %>%
    dplyr::arrange(grid_id)

  rateM_base <- mort_all$rateM_base_true
  rateM_scen <- mort_all$rateM_scen_true
  etaM_base  <- mort_all$etaM_base_true
  etaM_scen  <- mort_all$etaM_scen_true
  mort_deaths_base_true <- mort_all$mort_deaths_base_true

  lamM_hist <- mort_deaths_base_true[idx_hist]
  deaths_hist <- stats::rpois(sum(idx_hist), lambda = pmax(lamM_hist, 1e-12))

  mort_hist <- grid[idx_hist, ] %>%
    dplyr::mutate(
      deaths = deaths_hist,
      exposure = exposure,
      cause = "simulated"
    ) %>%
    dplyr::select(sex, age, period, cohort, deaths, exposure, cause)

  inc_truth_grid <- grid %>%
    dplyr::mutate(
      etaI_apc_true = etaI_apc,
      etaI_z_true = etaI_z_base,
      etaI_total_true = etaI_base,
      rateI_base_true = rateI_base,
      rateI_scen_true = rateI_scen,
      aI_true = aI[as.character(age)],
      pI_true = pI[as.character(period)],
      cI_true = cI[as.character(cohort)],
      zI_true_used = z_base_used$z_prev
    ) %>%
    dplyr::select(sex, age, period, cohort, grid_id,
                  etaI_apc_true, etaI_z_true, etaI_total_true,
                  rateI_base_true, rateI_scen_true,
                  aI_true, pI_true, cI_true, zI_true_used)

  mort_truth_grid <- grid %>%
    dplyr::left_join(
      mort_all %>%
        dplyr::select(sex, age, period, cohort, grid_id,
                      etaM_base_true, etaM_scen_true,
                      rateM_base_true, rateM_scen_true,
                      mort_deaths_base_true, mort_deaths_scen_true),
      by = c("sex", "age", "period", "cohort", "grid_id")
    ) %>%
    dplyr::mutate(
      logIker_base_true = NA_real_,
      logIker_scen_true = NA_real_,
      aM_true = NA_real_,
      pM_true = NA_real_,
      cM_true = NA_real_
    ) %>%
    dplyr::select(sex, age, period, cohort, grid_id,
                  etaM_base_true, etaM_scen_true,
                  rateM_base_true, rateM_scen_true,
                  logIker_base_true, logIker_scen_true,
                  aM_true, pM_true, cM_true,
                  mort_deaths_base_true, mort_deaths_scen_true)

  idx_fut <- grid$period > last_hist
  truth <- list(
    A_true = A_true,
    L_kernel = wL,
    L_center = L_center,
    base_year = base_year,
    scen_rate = scen_rate,
    theta_I = if (identical(beta_mode, "fixed_rr_offset")) 1 else theta_I,
    rr_I = rr_by_sex,
    prev_inc_mode = beta_mode,
    theta_M = NA_real_,
    prev_base_by_sex = p0_base,
    age_min = age_min,
    age_max = age_max,
    inc_delta_true = sum(exposure * (rateI_scen[idx_fut] - rateI_base[idx_fut])),
    mort_delta_true = sum(mort_all$mort_deaths_scen_true[idx_fut] - mort_all$mort_deaths_base_true[idx_fut], na.rm = TRUE)
  )

  list(
    prev_micro = prev_micro,
    prev_truth = prev_truth,
    gammaP_true = gammaP_true,
    z_base_true = z_base,
    z_scen_true = z_scen,
    zI_base_true_used = zI_base_true_used,
    zI_scen_true_used = zI_scen_true_used,
    inc_truth_grid = inc_truth_grid,
    mort_truth_grid = mort_truth_grid,
    inc_hist = inc_hist,
    pop_all = pop_all,
    mort_hist = mort_hist,
    truth = truth,
    meta = list(seed = seed, dgp = dgp, last_hist = last_hist, base_year = base_year,
                n_prev_micro = N_prev_micro, exposure = exposure,
                age_min = age_min, age_max = age_max)
  )
}

build_inputs_sim <- function(sim,
                             cause_id = "simulated",
                             label = "simulated") {
  stopifnot(is.list(sim))
  req <- c("prev_micro", "inc_hist", "pop_all", "mort_hist")
  miss <- setdiff(req, names(sim))
  if (length(miss)) stop("Faltan componentes en sim: ", paste(miss, collapse = ", "))
  make_bapc_inputs(
    mort_hist_tbl = sim$mort_hist,
    pop_all_tbl   = sim$pop_all,
    inc_hist_tbl  = sim$inc_hist,
    prev_path     = NULL,
    prev_data     = sim$prev_micro,
    metadata = c(
      list(cause_id = cause_id, label = label, source = "simulation"),
      if (is.null(sim$meta)) list() else sim$meta
    )
  )
}
# build_prev_index_for_inc: helper para simulaciones (usa gammaP_all, A_I, w_I) y devuelve z_prev
build_prev_index_for_inc <- function(df_inc_grid, gammaP_all, A_I, w_I = 1, base_year = 2022) {
  # Re-implementación simplificada para simulación
  grid <- df_inc_grid %>%
    dplyr::mutate(cohort_prev = period - A_I - age) %>%
    dplyr::left_join(gammaP_all, by = c("cohort_prev" = "cohort")) %>%
    dplyr::mutate(z_prev = w_I * dplyr::coalesce(gammaP, 0))
  
  # Calcular s_histP (desvío estándar histórico para normalización si fuera necesario)
  # Aquí lo hacemos simple: z_prev es directamente el efecto de cohorte rezagado
  attr(grid, "s_histP") <- 1
  grid
}

# build_prev_rr_offset_for_inc: helper para simulaciones en modo fixed_rr_offset
build_prev_rr_offset_for_inc <- function(df_inc_grid, gammaP_all, A_I, rr_inc, prev_base_prob, base_year = 2022) {
  req <- c("age", "period")
  miss <- setdiff(req, names(df_inc_grid))
  if (length(miss)) stop("build_prev_rr_offset_for_inc: faltan columnas: ", paste(miss, collapse = ", "))
  
  # 1) Índice base (efecto cohorte rezagado)
  idx <- build_prev_index_for_inc(df_inc_grid = df_inc_grid, gammaP_all = gammaP_all, A_I = A_I, w_I = 1, base_year = base_year)
  
  # 2) Anclaje a prevalencia base
  # Simplificación: asumimos que gammaP ya está centrada o que prev_base_prob captura el nivel medio
  out <- idx %>%
    dplyr::mutate(sex_chr = as.character(sex))
  
  if (length(prev_base_prob) == 1L) {
    out$base_p <- as.numeric(prev_base_prob)
  } else {
    out$base_p <- as.numeric(prev_base_prob[out$sex_chr])
  }
  
  if (length(rr_inc) == 1L) {
    out$rr_val <- as.numeric(rr_inc)
  } else {
    out$rr_val <- as.numeric(rr_inc[out$sex_chr])
  }
  
  out <- out %>%
    dplyr::mutate(
      q_eff = plogis(qlogis(pmin(pmax(base_p, 1e-6), 1-1e-6)) + z_prev),
      z_prev_offset = log1p(q_eff * (rr_val - 1))
    ) %>%
    dplyr::select(dplyr::all_of(names(df_inc_grid)), q_eff, z_prev = z_prev_offset)
  
  out
}
