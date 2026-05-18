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

.sim_cell_seed <- function(base_seed, sex, age, period, cohort) {
  modulus <- 2147483000
  sex_code <- ifelse(as.character(sex) == "F", 2, 1)
  raw <- as.numeric(base_seed) * 1000003 +
    as.numeric(sex_code) * 9176 +
    as.numeric(age) * 1009 +
    as.numeric(period) * 917 +
    as.numeric(cohort) * 101
  as.integer((raw %% modulus) + 1)
}

simulate_prev_micro <- function(prev_truth, n_cell = 200L, cell_seed = NULL) {
  stopifnot(all(c("sex","age","period","cohort","p_true") %in% names(prev_truth)))
  idx <- rep(seq_len(nrow(prev_truth)), each = as.integer(n_cell))
  df <- prev_truth[idx, c("sex","age","period","cohort","p_true")]
  prob <- pmax(pmin(df$p_true, 1 - 1e-8), 1e-8)
  fuma <- if (is.null(cell_seed)) {
    stats::rbinom(nrow(df), size = 1L, prob = prob)
  } else {
    draws <- vector("list", nrow(prev_truth))
    n_cell_int <- as.integer(n_cell)
    for (i in seq_len(nrow(prev_truth))) {
      set.seed(.sim_cell_seed(
        base_seed = cell_seed,
        sex = prev_truth$sex[i],
        age = prev_truth$age[i],
        period = prev_truth$period[i],
        cohort = prev_truth$cohort[i]
      ))
      draws[[i]] <- stats::rbinom(
        n_cell_int,
        size = 1L,
        prob = pmax(pmin(prev_truth$p_true[i], 1 - 1e-8), 1e-8)
      )
    }
    unlist(draws, use.names = FALSE)
  }
  tibble::tibble(
    period = as.integer(df$period),
    age    = as.integer(df$age),
    cohort = as.integer(df$cohort),
    sex    = factor(as.character(df$sex), levels = c("M", "F")),
    fuma   = fuma,
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
                              age_min_p = AGE_P_MIN, age_max_p = AGE_P_MAX,
                              per_min = 1990, per_max = PROJ_TO,
                              last_hist = 2022,
                              N_prev_micro = 200L,
                              exposure = 100000,
                              A_true = 16L,
                              Lmax = 7L, L_center = 2.5, L_sd = 1.2,
                              theta_I = 0.35,
                              rr_I = NULL,
                              theta_M = 1.0,
                              dgp = c("spec_linear", "misspec_tanh", "weak"),
                              kappa = 0.8,
                              scen_rate = 0.01,
                              scenario_name = c("freeze","up1pc","down1pc","down3pc","quit"),
                              prev_annual_rate_up = PREV_ANNUAL_RATE_UP,
                              prev_annual_rate_down = PREV_ANNUAL_RATE_DOWN,
                              prev_annual_rate_down3 = PREV_ANNUAL_RATE_DOWN3,
                              prev_obs_period_min = per_min,
                              prev_obs_period_max = per_max,
                              prev_obs_age_min = age_min_p,
                              prev_obs_age_max = age_max_p,
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
  set_stage_seed <- function(offset) {
    set.seed(as.integer(seed)[1] + as.integer(offset)[1])
  }
  dgp <- match.arg(dgp)
  scenario_name <- match.arg(as.character(scenario_name)[1], c("freeze","up1pc","down1pc","down3pc","quit"))
  beta_mode <- "fixed_rr_offset"
  sexes <- intersect(c("M", "F"), unique(as.character(sexes)))
  if (!length(sexes)) stop("sexes must include M and/or F")
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

  # Extract true RR from project configuration for the specified cause
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

  stock_per_min <- per_min - max(0L, as.integer(age_max) - as.integer(age_min_p))
  prev_stock_grid <- tidyr::expand_grid(
    sex = factor(sexes, levels = c("M", "F")),
    age = age_min_p:age_max,
    period = stock_per_min:per_max
  ) %>% dplyr::mutate(cohort = period - age)
  prev_obs_period_min <- max(stock_per_min, as.integer(prev_obs_period_min)[1])
  prev_obs_period_max <- min(per_max, as.integer(prev_obs_period_max)[1])
  prev_obs_age_min <- max(age_min_p, as.integer(prev_obs_age_min)[1])
  prev_obs_age_max <- min(age_max, as.integer(prev_obs_age_max)[1])
  prev_obs_grid <- prev_stock_grid %>%
    dplyr::filter(
      age >= prev_obs_age_min, age <= prev_obs_age_max,
      period >= prev_obs_period_min, period <= prev_obs_period_max
    )

  lev_age <- sort(unique(grid$age))
  lev_per <- sort(unique(grid$period))
  lev_coh <- sort(unique(grid$cohort))

  # --------- PREVALENCE (Binomial DGP with APC structure) ----------
  lev_age_p <- sort(unique(prev_stock_grid$age))
  lev_per_p <- sort(unique(prev_stock_grid$period))
  lev_coh_p <- sort(unique(prev_stock_grid$cohort))
  aP_raw <- .make_rw1(length(lev_age_p), 0.08)
  idx_age_cap <- max(which(lev_age_p <= age_max_p))
  aP_raw[lev_age_p > age_max_p] <- aP_raw[idx_age_cap]
  aP <- aP_raw; names(aP) <- lev_age_p

  message(sprintf(">>> DGP SEED %s | aP[1:3]: %s", seed, paste(round(head(aP, 3), 4), collapse=" ")))
  pP <- .make_rw1(length(lev_per_p), 0.06); names(pP) <- lev_per_p
  cP <- .make_rw1(length(lev_coh_p), 0.10); names(cP) <- lev_coh_p
  sex_shift_prev <- prev_sex_shift[as.character(prev_stock_grid$sex)]
  sex_shift_prev[is.na(sex_shift_prev)] <- 0

  stock_age_clamped <- pmin(pmax(prev_stock_grid$age, age_min_p), age_max_p)
  
  etaP_stock <- -2.0 + sex_shift_prev +
    aP[as.character(stock_age_clamped)] +
    pP[as.character(prev_stock_grid$period)] +
    cP[as.character(prev_stock_grid$cohort)]
  prev_stock <- prev_stock_grid %>%
    dplyr::mutate(p_true = plogis(etaP_stock))

  prev_truth <- prev_obs_grid %>%
    dplyr::left_join(prev_stock, by = c("sex", "age", "period", "cohort")) %>%
    dplyr::transmute(sex, age, period, cohort, p_true = p_true)
  set_stage_seed(100000L)
  prev_micro <- simulate_prev_micro(prev_truth = prev_truth, n_cell = N_prev_micro, cell_seed = seed)

  gammaP_true <- tibble::tibble(
    cohort = lev_coh_p,
    gammaP = as.numeric(cP)
  )

  p0_base <- prev_truth %>%
    dplyr::filter(period == base_year) %>%
    dplyr::group_by(sex) %>%
    dplyr::summarise(p0 = mean(p_true, na.rm = TRUE), .groups = "drop") %>%
    tibble::deframe()

  apply_prev_scenario <- function(prev_df, scenario) {
    out <- prev_df %>% dplyr::mutate(p_curr = p_true)
    k <- pmax(0L, out$period - base_year)
    idx_fut <- out$period > base_year
    ref_tbl <- out %>%
      dplyr::filter(period == base_year) %>%
      dplyr::select(sex, age, p_ref = p_true) %>%
      dplyr::distinct(sex, age, .keep_all = TRUE)
    out <- out %>% dplyr::left_join(ref_tbl, by = c("sex", "age"))
    idx_ref <- idx_fut & is.finite(out$p_ref)
    if (identical(scenario, "freeze")) {
      out$p_curr[idx_ref] <- out$p_ref[idx_ref]
    } else if (identical(scenario, "up1pc")) {
      fac <- (1 + abs(prev_annual_rate_up))^k
      out$p_curr[idx_ref] <- out$p_ref[idx_ref] * fac[idx_ref]
    } else if (identical(scenario, "down1pc")) {
      fac <- (1 - abs(prev_annual_rate_down))^k
      out$p_curr[idx_ref] <- out$p_ref[idx_ref] * fac[idx_ref]
    } else if (identical(scenario, "down3pc")) {
      fac <- (1 - abs(prev_annual_rate_down3))^k
      out$p_curr[idx_ref] <- out$p_ref[idx_ref] * fac[idx_ref]
    } else if (identical(scenario, "quit")) {
      out$p_curr[idx_fut] <- 0
    }
    out$p_curr <- pmin(pmax(out$p_curr, 1e-8), 1 - 1e-8)
    out %>% dplyr::select(-p_ref)
  }

  build_z_from_prev_stock <- function(prev_df) {
    rr_vec <- suppressWarnings(as.numeric(rr_by_sex[as.character(grid$sex)]))
    rr_vec[!is.finite(rr_vec) | rr_vec <= 1] <- rr_scalar_default
    
    z_stock <- prev_df %>%
      dplyr::arrange(sex, cohort, age) %>%
      dplyr::group_by(sex, cohort) %>%
      dplyr::group_modify(function(df, key) {
        sx <- as.character(key$sex)
        L_q <- get_prev_inc_quit_horizon(cause_id = cause_id, sex = sx)
        w_vec <- .prev_rr_schedule_vec(cause_id = cause_id, sex = sx, L_q = L_q)
        rr_use <- as.numeric(rr_by_sex[sx])
        
        df$q_eff <- .calculate_stock_former_q_eff(as.numeric(df$p_curr), w_vec, L_q)
        df$z_prev <- log1p(df$q_eff * (rr_use - 1))
        df
      }) %>%
      dplyr::ungroup()

    z_stock %>%
      dplyr::select(dplyr::any_of(c("sex", "age", "period", "cohort", "p_curr", "q_eff", "z_prev"))) %>%
      dplyr::inner_join(grid %>% dplyr::select(sex, age, period, cohort, grid_id),
                        by = c("sex", "age", "period", "cohort"))
  }

  # --------- True z_prev (baseline and scenario) ----------
  z_base <- if (identical(beta_mode, "fixed_rr_offset")) {
    build_z_from_prev_stock(apply_prev_scenario(prev_stock, "freeze"))
  } else {
    build_prev_index_for_inc(df_inc_grid = grid, gammaP_all = gammaP_true, A_I = A_true, w_I = 1)
  } %>%
    dplyr::select(dplyr::any_of(c("sex", "age", "period", "cohort", "grid_id", "p_curr", "q_eff", "z_prev"))) %>%
    dplyr::right_join(grid %>% dplyr::select(sex, age, period, cohort, grid_id),
                      by = c("sex", "age", "period", "cohort", "grid_id")) %>%
    dplyr::arrange(grid_id)
  if (anyNA(z_base$z_prev)) stop("simulate_PIM_data: z_base has NAs after aligning by grid_id.")

  s_hist <- stats::sd(z_base$z_prev[z_base$period <= base_year], na.rm = TRUE)
  unit <- 1 / pmax(s_hist, 1e-9)
  message(sprintf(">>> simulate_PIM_data: beta_mode='%s' scenario='%s'", beta_mode, scenario_name))
  if (identical(beta_mode, "fixed_rr_offset")) {
    message(">>> simulate_PIM_data: entered fixed_rr_offset block")
    z_scen <- build_z_from_prev_stock(apply_prev_scenario(prev_stock, scenario_name)) %>%
      dplyr::right_join(grid %>% dplyr::select(sex, age, period, cohort, grid_id),
                        by = c("sex", "age", "period", "cohort", "grid_id")) %>%
      dplyr::arrange(grid_id)
  } else {
    dz <- .sim_delta_z_period(z_base, p0_by_sex = p0_base, base_year = base_year, rate = scen_rate, unit = unit)
    z_scen <- z_base %>% dplyr::mutate(z_prev = z_prev + dz) %>% dplyr::arrange(grid_id)
  }

  # Preserve the exact z regressors used in the DGP incidence equation
  zI_base_true_used <- z_base %>% dplyr::transmute(sex, age, period, cohort, grid_id, z_prev_used = z_prev)
  zI_scen_true_used <- z_scen %>% dplyr::transmute(sex, age, period, cohort, grid_id, z_prev_used = z_prev)

  # --------- INCIDENCE (Poisson DGP with APC + g(z)) ----------
  set_stage_seed(200000L)
  aI <- .make_rw1(length(lev_age), 0.02); names(aI) <- lev_age
  # Stabilize projection period so that scenarios are visible
  pI_raw <- .make_rw1(length(lev_per), 0.015)
  # Do not flatten: let the RW1 trend continue into the future
  # so that Truth and Estimator share the same dynamic structure.
  pI <- pI_raw; names(pI) <- lev_per
  cI <- .make_rw1(length(lev_coh), 0.03); names(cI) <- lev_coh
  sex_shift_inc <- inc_sex_shift[as.character(grid$sex)]
  sex_shift_inc[is.na(sex_shift_inc)] <- 0

  # Synchronize intercept:
  # If the estimator uses offset(z - mean_z), the estimated intercept is log_rate(mean_z).
  # The DGP uses exp(-8.5 + sex_shift + z).
  # For them to align, the DGP should generate cases consistent with a stable baseline.
  etaI_apc <- -8.5 + sex_shift_inc +
    aI[as.character(grid$age)] + pI[as.character(grid$period)] + cI[as.character(grid$cohort)]

  z_base_used <- grid %>%
    dplyr::select(sex, age, period, cohort, grid_id) %>%
    dplyr::left_join(z_base, by = c("sex", "age", "period", "cohort", "grid_id")) %>%
    dplyr::arrange(grid_id)
    
  # IMPORTANT: Ensure z_prev in the DGP is exactly what the estimator will use
  # but with non-linearity if dgp == "misspec_tanh"
  g_z <- function(z, dgp_type) {
    if (dgp_type == "misspec_tanh") {
      # We apply a tanh-type saturation to break linearity
      # Scale z by its historical standard deviation so that the effect is visible
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
  set_stage_seed(300000L)
  cases_hist <- stats::rpois(sum(idx_hist), lambda = lamI_hist)

  inc_hist <- grid[idx_hist, ] %>%
    dplyr::mutate(cases = cases_hist) %>%
    dplyr::select(sex, age, period, cohort, cases)

  pop_all <- grid %>%
    dplyr::transmute(sex, age, period, exposure = exposure)

  # --------- MORTALITY (CONSISTENT DGP: INCIDENCE ⊗ KERNEL) ----------
  # Extract reference kernel for metadata (using M as default)
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
  set_stage_seed(400000L)
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
                age_min = age_min, age_max = age_max,
                age_min_p = age_min_p, age_max_p = age_max_p,
                prev_obs_period_min = prev_obs_period_min,
                prev_obs_period_max = prev_obs_period_max,
                prev_obs_age_min = prev_obs_age_min,
                prev_obs_age_max = prev_obs_age_max,
                stock_period_min = stock_per_min)
  )
}

build_inputs_sim <- function(sim,
                             cause_id = "simulated",
                             label = "simulated") {
  stopifnot(is.list(sim))
  req <- c("prev_micro", "inc_hist", "pop_all", "mort_hist")
  miss <- setdiff(req, names(sim))
  if (length(miss)) stop("Missing components in sim: ", paste(miss, collapse = ", "))
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
# build_prev_index_for_inc: helper for simulations (uses gammaP_all, A_I, w_I) and returns z_prev
build_prev_index_for_inc <- function(df_inc_grid, gammaP_all, A_I, w_I = 1, base_year = 2022) {
  # Simplified re-implementation for simulation
  grid <- df_inc_grid %>%
    dplyr::mutate(cohort_prev = period - A_I - age) %>%
    dplyr::left_join(gammaP_all, by = c("cohort_prev" = "cohort")) %>%
    dplyr::mutate(z_prev = w_I * dplyr::coalesce(gammaP, 0))
  
  # Calculate s_histP (historical standard deviation for normalization if needed)
  # Here we keep it simple: z_prev is directly the lagged cohort effect
  attr(grid, "s_histP") <- 1
  grid
}

# build_prev_rr_offset_for_inc: helper for simulations in fixed_rr_offset mode
build_prev_rr_offset_for_inc <- function(df_inc_grid, gammaP_all, A_I, rr_inc, prev_base_prob, base_year = 2022) {
  req <- c("age", "period")
  miss <- setdiff(req, names(df_inc_grid))
  if (length(miss)) stop("build_prev_rr_offset_for_inc: missing columns: ", paste(miss, collapse = ", "))
  
  # 1) Base index (lagged cohort effect)
  idx <- build_prev_index_for_inc(df_inc_grid = df_inc_grid, gammaP_all = gammaP_all, A_I = A_I, w_I = 1, base_year = base_year)
  
  # 2) Anchor to base prevalence
  # Simplification: we assume that gammaP is already centered or that prev_base_prob captures the mean level
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
