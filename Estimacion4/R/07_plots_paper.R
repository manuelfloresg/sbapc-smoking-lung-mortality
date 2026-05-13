# =============================================================
# Plot functions with publication-ready aesthetics (first cut)
# =============================================================
# 7) Gráficos
# =============================================================

plot_apc_panel <- function(sa, sp, sc, lev, title_txt) {
  df_age <- tibble(x = lev$age[sa$ID],    fit = sa$mean - mean(sa$mean), lwr = sa$`0.025quant`, upr = sa$`0.975quant`)
  df_per <- tibble(x = lev$period[sp$ID], fit = sp$mean - mean(sp$mean), lwr = sp$`0.025quant`, upr = sp$`0.975quant`)
  df_coh <- tibble(x = lev$cohort[sc$ID], fit = sc$mean - mean(sc$mean), lwr = sc$`0.025quant`, upr = sc$`0.975quant`)
  pA <- ggplot(df_age, aes(x, fit)) + geom_ribbon(aes(ymin=lwr,ymax=upr), alpha=.15) + geom_line() +
    labs(x="Age", y="Effect (link)") + theme_minimal(12)
  pP <- ggplot(df_per, aes(x, fit)) + geom_ribbon(aes(ymin=lwr,ymax=upr), alpha=.15) + geom_line() +
    labs(x="Period", y=NULL) + theme_minimal(12)
  pC <- ggplot(df_coh, aes(x, fit)) + geom_ribbon(aes(ymin=lwr,ymax=upr), alpha=.15) + geom_line() +
    labs(x="Cohort", y="Effect (link)") + theme_minimal(12)
  (pA + pP) / pC + plot_annotation(title = title_txt)
}

plot_apc_prev <- function(res_sex) {
  fit <- res_sex$fit_prev; lev <- res_sex$lev_prev
  plot_apc_panel(fit$summary.random$age_id, fit$summary.random$period_id, fit$summary.random$cohort_id,
                 lev, paste0("Prevalence APC (INLA, logit) — ", ifelse(res_sex$sex=="M","Males","Females")))
}

plot_apc_prev_both <- function(res_both,
                               base_size = 11,
                               male_lab = "Males",
                               female_lab = "Females") {
  stopifnot(!is.null(res_both$resM), !is.null(res_both$resF))
  make_df <- function(sumrnd, lev_vec, sex_lab) {
    mu <- mean(sumrnd$mean)
    tibble::tibble(
      x   = lev_vec[sumrnd$ID],
      fit = sumrnd$mean - mu,
      lwr = sumrnd$`0.025quant` - mu,
      upr = sumrnd$`0.975quant` - mu,
      sex = sex_lab
    )
  }
  fitM <- res_both$resM$fit_prev; levM <- res_both$resM$lev_prev
  fitF <- res_both$resF$fit_prev; levF <- res_both$resF$lev_prev
  if (is.null(fitM$summary.random$age_id) || is.null(fitF$summary.random$age_id)) {
    stop("No encuentro summary.random$age_id en fit_prev (¿corrió el modelo de prevalencia?).")
  }
  df_age <- dplyr::bind_rows(
    make_df(fitM$summary.random$age_id,    levM$age,    male_lab),
    make_df(fitF$summary.random$age_id,    levF$age,    female_lab)
  )
  df_per <- dplyr::bind_rows(
    make_df(fitM$summary.random$period_id, levM$period, male_lab),
    make_df(fitF$summary.random$period_id, levF$period, female_lab)
  )
  df_coh <- dplyr::bind_rows(
    make_df(fitM$summary.random$cohort_id, levM$cohort, male_lab),
    make_df(fitF$summary.random$cohort_id, levF$cohort, female_lab)
  )
  # Colores pedidos: hombres azul, mujeres verde
  cols <- c(setNames("#1f77b4", male_lab), setNames("#2ca02c", female_lab))
  theme_common <- ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5, size = base_size * 0.95, face = "plain"),
      legend.position = "bottom",
      legend.title = ggplot2::element_blank()
    )
  pA <- ggplot2::ggplot(df_age, ggplot2::aes(x = x, y = fit, color = sex, fill = sex)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lwr, ymax = upr), alpha = 0.20, colour = NA) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
    ggplot2::scale_color_manual(values = cols) +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::labs(title = "Age", x = NULL, y = NULL) +
    theme_common +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = base_size * 0.70))
  pP <- ggplot2::ggplot(df_per, ggplot2::aes(x = x, y = fit, color = sex, fill = sex)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lwr, ymax = upr), alpha = 0.20, colour = NA) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
    ggplot2::scale_color_manual(values = cols) +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::labs(title = "Period", x = NULL, y = NULL) +
    theme_common +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = base_size * 0.70))
  pC <- ggplot2::ggplot(df_coh, ggplot2::aes(x = x, y = fit, color = sex, fill = sex)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lwr, ymax = upr), alpha = 0.20, colour = NA) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
    ggplot2::scale_color_manual(values = cols) +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::labs(title = "Cohort", x = NULL, y = NULL) +
    ggplot2::scale_x_continuous(breaks = function(x) seq(floor(min(x)/10)*10, ceiling(max(x)/10)*10, by = 10)) +
    theme_common +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = base_size * 0.70))
  
  (pA + pP) / pC +
    patchwork::plot_layout(guides = "collect") &
    ggplot2::theme(legend.position = "bottom")
}
save_apc_prev_both_png <- function(res_both, file, width = 10, height = 6, dpi = 300, ...) {
  g <- plot_apc_prev_both(res_both, ...)
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  ggplot2::ggsave(filename = file, plot = g, width = width, height = height, dpi = dpi, bg = "white")
  invisible(file)
}

plot_ai_calibration_both <- function(res_both,
                                     base_size = 11,
                                     male_lab = "Males",
                                     female_lab = "Females") {
  
  tabM <- tryCatch(res_both$resM$diag$cal_AI_table, error = function(e) NULL)
  tabF <- tryCatch(res_both$resF$diag$cal_AI_table, error = function(e) NULL)
  
  AstarM <- tryCatch(res_both$resM$diag$A_I_star, error = function(e) NA_real_)
  AstarF <- tryCatch(res_both$resF$diag$A_I_star, error = function(e) NA_real_)
  
  if (is.null(tabM) && is.null(tabF)) {
    stop("No encuentro diag$cal_AI_table en res_both (ni M ni F).")
  }
  
  dfs <- list()
  cols <- c()
  
  if (!is.null(tabM)) {
    dfs[[male_lab]] <- dplyr::mutate(tabM, sex = male_lab)
    cols <- c(cols, setNames("#1f77b4", male_lab))
  }
  if (!is.null(tabF)) {
    dfs[[female_lab]] <- dplyr::mutate(tabF, sex = female_lab)
    cols <- c(cols, setNames("#2ca02c", female_lab))
  }
  
  df <- dplyr::bind_rows(dfs) %>%
    dplyr::filter(is.finite(A), is.finite(abs_r))
  
  df_star <- dplyr::bind_rows(
    if (!is.null(tabM)) dplyr::filter(df, sex == male_lab   & A == AstarM) else NULL,
    if (!is.null(tabF)) dplyr::filter(df, sex == female_lab & A == AstarF) else NULL
  )
  
  ggplot2::ggplot(df, ggplot2::aes(x = A, y = abs_r, color = sex)) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = 1.6) +
    ggplot2::geom_point(data = df_star, size = 3.2) +
    ggplot2::scale_color_manual(values = cols) +
    ggplot2::labs(
      x = expression("Exposure-history horizon " * A[I] * " (years)"),
      y = expression("| Cor(" * tilde(z)[t] * "(" * A[I] * "), " * hat(beta)[t]^I * ") |")
    ) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.title = ggplot2::element_blank()
    )
}


save_ai_calibration_png <- function(res_both, file, width = 10, height = 4, dpi = 300, ...) {
  g <- plot_ai_calibration_both(res_both, ...)
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  ggplot2::ggsave(filename = file, plot = g, width = width, height = height, dpi = dpi, bg = "white")
  invisible(file)
}

save_prev_apc_global <- function(res_any, out_base, file_name = "prev_apc_components_bothsex.png") {
  dir_global <- file.path(out_base, "plots")
  dir.create(dir_global, recursive = TRUE, showWarnings = FALSE)
  
  save_apc_prev_both_png(
    res_any,
    file = file.path(dir_global, file_name),
    width = 10, height = 6, dpi = 300
  )
}


compute_period_with_trend <- function(fit, lev_period, coef_prefix, degree, center) {
  cf <- fit$summary.fixed
  b1 <- if (paste0(coef_prefix,"trend_t") %in% rownames(cf))  cf[paste0(coef_prefix,"trend_t"),"mean"]  else 0
  b2 <- if (degree >= 2 && paste0(coef_prefix,"trend_t2") %in% rownames(cf)) cf[paste0(coef_prefix,"trend_t2"),"mean"] else 0
  t  <- if (degree >= 1) (lev_period - center) else 0
  t2 <- if (degree >= 2) (lev_period - center)^2 else 0
  b1 * t + b2 * t2
}

# === Helper: banda posterior para (período + tendencia) vía muestreo ===
.detect_degree <- function(fit, coef_prefix = "") {
  nmf <- rownames(fit$summary.fixed)
  has_t  <- paste0(coef_prefix, "trend_t")  %in% nmf
  has_t2 <- paste0(coef_prefix, "trend_t2") %in% nmf
  if (has_t2) 2 else if (has_t) 1 else 0
}

.make_idx_fixed <- function(nm, term) {
  idx <- which(nm == term)
  if (length(idx) == 0) {
    # Respaldo: a veces INLA anota "Beta for <term>" u otras variantes
    g <- grep(term, nm, fixed = TRUE)
    if (length(g)) idx <- g[1]
  }
  idx
}

make_period_band_sum <- function(fit, lev_period, center,
                                 nsamp = 1000, degree = NULL,
                                 coef_prefix = "") {
  stopifnot(!is.null(fit$summary.random$period_id))
  if (is.null(degree)) degree <- .detect_degree(fit, coef_prefix)
  
  sp <- fit$summary.random$period_id
  levP_hist <- lev_period[sp$ID]
  nP <- length(levP_hist)
  
  t  <- if (degree >= 1) (levP_hist - center) else rep(0, nP)
  t2 <- if (degree >= 2) (levP_hist - center)^2 else rep(0, nP)
  
  S  <- INLA::inla.posterior.sample(nsamp, fit)
  nm <- rownames(S[[1]]$latent)
  
  # índices del RW2 de período en el vector latente
  idx_per_all <- grep("^period_id", nm)
  
  # si agregaste iid de período, aparecerá aquí (si no, vector vacío)
  idx_iid_all <- grep("^period_iid", nm)
  
  # nombres de las betas fijas
  b1_name <- paste0(coef_prefix, "trend_t")
  b2_name <- paste0(coef_prefix, "trend_t2")
  
  # helper robusto: usa coincidencia exacta o "Beta for <term>" como fallback
  .idx_fixed <- function(term) {
    i <- which(nm == term)
    if (!length(i)) {
      # INLA suele anotar "Beta for term"
      i <- which(nm == paste0("Beta for ", term))
    }
    if (!length(i)) NA_integer_ else i[1]
  }
  
  draw_beta <- function(term) {
    i <- .idx_fixed(term)
    if (length(i) && is.finite(i)) {
      vapply(S, function(s) as.numeric(s$latent[i, 1]), numeric(1))
    } else {
      # respaldo: muestrear de la MVN de fixed si existe
      if (!is.null(fit$misc$cov.fixed) && term %in% rownames(fit$misc$cov.fixed)) {
        mu <- fit$summary.fixed[term, "mean", drop = TRUE]
        sg <- fit$misc$cov.fixed[term, term, drop = TRUE]
        rnorm(nsamp, mean = mu, sd = sqrt(max(sg, 0)))
      } else {
        rep(0, nsamp)
      }
    }
  }
  
  b1_vec <- if (degree >= 1) draw_beta(b1_name) else rep(0, nsamp)
  b2_vec <- if (degree >= 2) draw_beta(b2_name) else rep(0, nsamp)
  
  M <- matrix(NA_real_, nsamp, nP)
  for (k in seq_len(nsamp)) {
    per_k <- as.numeric(S[[k]]$latent[idx_per_all])[seq_len(nP)]
    iid_k <- if (length(idx_iid_all) >= nP) as.numeric(S[[k]]$latent[idx_iid_all])[seq_len(nP)] else rep(0, nP)
    M[k, ] <- per_k + iid_k + b1_vec[k] * t + b2_vec[k] * t2
  }
  
  mu  <- colMeans(M)
  cst <- mean(mu)
  
  trend_mean <- if (degree >= 1) mean(b1_vec) * t else rep(0, nP)
  if (degree >= 2) trend_mean <- trend_mean + mean(b2_vec) * t2
  trend_mean <- trend_mean - mean(trend_mean)
  
  tibble::tibble(
    x     = levP_hist,
    fit   = mu - cst,
    lwr   = apply(M, 2, stats::quantile, 0.025) - cst,
    upr   = apply(M, 2, stats::quantile, 0.975) - cst,
    trend = trend_mean
  )
}

plot_apc_inc <- function(res_sex) {
  fit    <- if (!is.null(res_sex$inc_fit) && !is.null(res_sex$inc_fit$fit)) res_sex$inc_fit$fit else res_sex$fit_inc
  lev    <- if (!is.null(res_sex$inc_fit) && !is.null(res_sex$inc_fit$lev_inc)) res_sex$inc_fit$lev_inc else res_sex$lev_inc
  deg    <- if (!is.null(res_sex$inc_fit$trend_meta$degree)) res_sex$inc_fit$trend_meta$degree else INC_TREND_DEGREE
  center <- if (!is.null(res_sex$inc_fit$trend_meta$center)) res_sex$inc_fit$trend_meta$center else mean(lev$period)
  
  sa <- fit$summary.random$age_id
  sp <- fit$summary.random$period_id
  sc <- fit$summary.random$cohort_id
  
  # --- EDAD (centramos fit y también lwr/upr con el mismo offset) ---
  mA <- mean(sa$mean)
  df_age <- tibble::tibble(
    x   = lev$age[sa$ID],
    fit = sa$mean - mA,
    lwr = sa$`0.025quant` - mA,
    upr = sa$`0.975quant` - mA
  )
  
  # --- COHORTE (centrado coherente) ---
  mC <- mean(sc$mean)
  df_coh <- tibble::tibble(
    x   = lev$cohort[sc$ID],
    fit = sc$mean - mC,
    lwr = sc$`0.025quant` - mC,
    upr = sc$`0.975quant` - mC
  )
  
  # --- PERÍODO: banda RW + (opcional) tendencia, solo aquí la línea roja ---
  band_per <- make_period_band_sum(
    fit, lev$period,
    center = center, nsamp = 1000, degree = deg, coef_prefix = "inc_"
  )
  
  # If compute_period_with_trend() exists, evaluate it on the same period support used by band_per
  trend_line <- try(compute_period_with_trend(fit, band_per$x,
                                              coef_prefix = "inc_",
                                              degree = deg, center = center),
                    silent = TRUE)
  if (inherits(trend_line, "try-error") || length(trend_line) != nrow(band_per) || all(!is.finite(trend_line))) {
    # fallback: use the embedded trend if available; otherwise use zeros with matching length
    trend_line <- if ("trend" %in% names(band_per) && length(band_per$trend) == nrow(band_per)) band_per$trend else rep(0, nrow(band_per))
  }
  trend_line <- trend_line - mean(trend_line)
  
  # Desplazamos la banda sumándole la tendencia (negra = RW + tendencia)
  band_per$fit <- band_per$fit + trend_line
  band_per$lwr <- band_per$lwr + trend_line
  band_per$upr <- band_per$upr + trend_line
  
  # ----- GGs -----
  theme_fix <- ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    plot.margin      = grid::unit(c(2, 2, 2, 2), "pt")
  )
  
  pA <- ggplot2::ggplot(df_age, ggplot2::aes(x, fit)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lwr, ymax = upr), alpha = .15) +
    ggplot2::geom_line() +
    ggplot2::labs(title = "Age", x = "Age", y = "Effect (link)") +
    ggplot2::theme_minimal(12) + theme_fix
  
  pP <- ggplot2::ggplot(band_per, ggplot2::aes(x, fit)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lwr, ymax = upr), alpha = .15, show.legend = FALSE) +
    ggplot2::geom_line() +
    ggplot2::geom_line(ggplot2::aes(y = trend_line), linetype = "dashed", color = "red", linewidth = 0.8) +
    ggplot2::labs(title = "Period", x = "Period", y = NULL) +
    ggplot2::theme_minimal(12) + theme_fix
  
  pC <- ggplot2::ggplot(df_coh, ggplot2::aes(x, fit)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lwr, ymax = upr), alpha = .15) +
    ggplot2::geom_line() +
    ggplot2::labs(title = "Cohort", x = "Cohort", y = NULL) +
    ggplot2::theme_minimal(12) + theme_fix
  
  # Layout: (Age | Period) / Cohort (Cohort ocupa todo el ancho)
  (pA + pP) / pC + patchwork::plot_annotation(
    title = paste0("Incidence APC (INLA, log) — ", ifelse(res_sex$sex == "M", "Males", "Females"))
  )
}


plot_apc_mort <- function(res_sex) {
  fit    <- res_sex$fit_bapc
  lev    <- res_sex$lev_mort
  deg    <- MORT_TREND_DEGREE
  center <- mean(lev$period)
  
  sa <- fit$summary.random$age_id
  sp <- fit$summary.random$period_id
  sc <- fit$summary.random$cohort_id
  
  # EDAD centrada
  mA <- mean(sa$mean)
  df_age <- tibble::tibble(
    x   = lev$age[sa$ID],
    fit = sa$mean - mA,
    lwr = sa$`0.025quant` - mA,
    upr = sa$`0.975quant` - mA
  )
  
  # COHORTE centrada
  mC <- mean(sc$mean)
  df_coh <- tibble::tibble(
    x   = lev$cohort[sc$ID],
    fit = sc$mean - mC,
    lwr = sc$`0.025quant` - mC,
    upr = sc$`0.975quant` - mC
  )
  
  # PERÍODO: banda RW + tendencia roja (solo aquí)
  band_per <- make_period_band_sum(
    fit, lev$period, center = center, nsamp = 1000, degree = deg, coef_prefix = "mort_"
  )
  
  trend_line <- try(compute_period_with_trend(
    fit, band_per$x, coef_prefix = "mort_", degree = deg, center = center
  ), silent = TRUE)
  if (inherits(trend_line, "try-error") || length(trend_line) != nrow(band_per) || all(!is.finite(trend_line))) {
    trend_line <- if ("trend" %in% names(band_per) && length(band_per$trend) == nrow(band_per)) band_per$trend else rep(0, nrow(band_per))
  }
  trend_line <- trend_line - mean(trend_line)
  
  band_per$fit <- band_per$fit + trend_line
  band_per$lwr <- band_per$lwr + trend_line
  band_per$upr <- band_per$upr + trend_line
  
  theme_fix <- ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    plot.margin      = grid::unit(c(2, 2, 2, 2), "pt")
  )
  
  pA <- ggplot2::ggplot(df_age, ggplot2::aes(x, fit)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lwr, ymax = upr), alpha = .15) + ggplot2::geom_line() +
    ggplot2::labs(title = "Age", x = "Age", y = "Effect (link)") +
    ggplot2::theme_minimal(12) + theme_fix
  
  pP <- ggplot2::ggplot(band_per, ggplot2::aes(x, fit)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lwr, ymax = upr), alpha = .15) + ggplot2::geom_line() +
    ggplot2::geom_line(ggplot2::aes(y = trend_line), linetype = "dashed", color = "red", linewidth = 0.8) +
    ggplot2::labs(title = "Period", x = "Period", y = NULL) +
    ggplot2::theme_minimal(12) + theme_fix
  
  pC <- ggplot2::ggplot(df_coh, ggplot2::aes(x, fit)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lwr, ymax = upr), alpha = .15) + ggplot2::geom_line() +
    ggplot2::labs(title = "Cohort", x = "Cohort", y = NULL) +
    ggplot2::theme_minimal(12) + theme_fix
  
  (pA + pP) / pC + patchwork::plot_annotation(
    title = paste0("Mortality APC (INLA, log) — ", ifelse(res_sex$sex == "M", "Males", "Females"))
  )
}

maybe_save_prev_apc_global <- function(res_both, out_base,
                                       file_name = "prev_apc_components_bothsex.png") {
  dir_global <- file.path(out_base, "plots")
  dir.create(dir_global, recursive = TRUE, showWarnings = FALSE)
  
  f <- file.path(dir_global, file_name)
  if (!file.exists(f)) {
    save_apc_prev_both_png(res_both, file = f, width = 10, height = 6, dpi = 300)
  }
  invisible(f)
}


plot_projections_mort <- function(res_sex,
                                  title_base = "Lung cancer mortality: observed and projected") {
  
  sex_tag   <- ifelse(res_sex$sex == "M", "Males", "Females")
  last_year <- res_sex$diag$last_hist_year
  obs_M     <- res_sex$obs_annual
  y0        <- obs_M$obs[obs_M$period == last_year]
  horizon_frontier <- tryCatch(res_sex$diag$projection_horizon_frontier, error = function(e) NULL)
  max_year <- projection_max_year_from_res_sex(res_sex, policy = "endogenous_max", default = max(c(last_year, PROJ_TO), na.rm = TRUE))
  end_credible <- projection_max_year_from_res_sex(res_sex, policy = "credible", default = max_year)
  end_caution  <- projection_max_year_from_res_sex(res_sex, policy = "caution", default = max_year)
  
  bapc_proj    <- res_sex$annual_bapc              %>% dplyr::filter(period >= last_year + 1) %>% clip_to_year(max_year = max_year) %>% dplyr::mutate(model = "BAPC M")
  anchor_InoP  <- res_sex$annual_anchor_noP        %>% dplyr::filter(period >= last_year + 1) %>% clip_to_year(max_year = max_year) %>% dplyr::mutate(model = "BAPC M | I")
  anchor_Icond <- res_sex$annual_anchor            %>% dplyr::filter(period >= last_year + 1) %>% clip_to_year(max_year = max_year) %>% dplyr::mutate(model = "BAPC M | I | P")
  
  add_2022 <- function(df, name)
    dplyr::bind_rows(tibble::tibble(period = last_year, deaths_hat = y0,
                                    deaths_lwr = NA_real_, deaths_upr = NA_real_, model = name), df)
  
  bapc_plot  <- add_2022(bapc_proj,   "BAPC M")
  ibapc_plot <- add_2022(anchor_InoP, "BAPC M | I")
  icond_plot <- add_2022(anchor_Icond,"BAPC M | I | P")
  
  ggplot2::ggplot() +
    ggplot2::geom_ribbon(data = bapc_plot,  ggplot2::aes(period, ymin = deaths_lwr, ymax = deaths_upr),
                         fill = "grey80", alpha = .35, show.legend = FALSE) +
    ggplot2::geom_ribbon(data = ibapc_plot, ggplot2::aes(period, ymin = deaths_lwr, ymax = deaths_upr),
                         fill = "#ff7f0e", alpha = .15, show.legend = FALSE) +
    ggplot2::geom_ribbon(data = icond_plot, ggplot2::aes(period, ymin = deaths_lwr, ymax = deaths_upr),
                         fill = "#2ca02c", alpha = .18, show.legend = FALSE) +
    ggplot2::geom_line(data = icond_plot, ggplot2::aes(period, deaths_hat, color = model, linetype = model), linewidth = 1) +
    ggplot2::geom_line(data = ibapc_plot, ggplot2::aes(period, deaths_hat, color = model, linetype = model), linewidth = 1) +
    ggplot2::geom_line(data = bapc_plot,  ggplot2::aes(period, deaths_hat, color = model, linetype = model), linewidth = 1) +
    ggplot2::geom_line(data = obs_M, ggplot2::aes(period, obs), color = "black", linewidth = 1) +
    ggplot2::geom_vline(xintercept = last_year + 0.5, linetype = 3) +
    ggplot2::scale_color_manual(values = c("BAPC M" = "black",
                                           "BAPC M | I" = "#ff7f0e",
                                           "BAPC M | I | P" = "#2ca02c")) +
    ggplot2::scale_linetype_manual(values = c("BAPC M" = "dashed",
                                              "BAPC M | I" = "solid",
                                              "BAPC M | I | P" = "solid")) +
    ggplot2::labs(
      x = "Year", y = "Deaths", color = "Projection", linetype = "Projection",
      title    = paste0(title_base, " — ", sex_tag),
      subtitle = paste0("Observed through ", last_year,
                        " | Endogenous projection window: ", last_year + 1, "–", max_year,
                        " | Credible through ", end_credible,
                        " | Caution through ", end_caution),
      caption  = "Bands: 95% CI (transformed INLA quantiles)"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::scale_y_continuous(limits = c(0, NA), expand = c(0, 0)) +
    ggplot2::theme(legend.position = "bottom")
}


plot_projections_mort_total <- function(res_both,
                                        title = "Mortality — Total",
                                        show_hist_points = TRUE) {
  stopifnot(!is.null(res_both))
  
  if (missing(title) || is.null(title)) {
    lab <- attr(res_both, "label")
    if (!is.null(lab)) title <- paste0("Mortality — ", lab, " (Total)")
  }
  
  # 1) Total proyectado (seguro, viene de run_pipeline_both$combined)
  tot <- res_both$combined$annual_anchor
  if (is.null(tot) || !all(c("period","deaths_hat","deaths_lwr","deaths_upr") %in% names(tot))) {
    # Fallback: intentar sumar manualmente si por alguna razón no está combined
    ab_m <- if (!is.null(res_both$resM$annual_anchor)) res_both$resM$annual_anchor else NULL
    ab_f <- if (!is.null(res_both$resF$annual_anchor)) res_both$resF$annual_anchor else NULL
    if (is.null(ab_m) && is.null(ab_f)) {
      stop("plot_projections_mort_total(): no hay anual_anchor disponible ni en combined ni por sexo.")
    }
    # normalizar columnas y sumar
    if (!is.null(ab_m)) ab_m <- dplyr::rename(ab_m, deaths_hat_M=deaths_hat, deaths_lwr_M=deaths_lwr, deaths_upr_M=deaths_upr)
    if (!is.null(ab_f)) ab_f <- dplyr::rename(ab_f, deaths_hat_F=deaths_hat, deaths_lwr_F=deaths_lwr, deaths_upr_F=deaths_upr)
    tot <- dplyr::full_join(ab_m, ab_f, by="period") %>%
      dplyr::mutate(
        deaths_hat = dplyr::coalesce(.data$deaths_hat_M,0) + dplyr::coalesce(.data$deaths_hat_F,0),
        deaths_lwr = dplyr::coalesce(.data$deaths_lwr_M,0) + dplyr::coalesce(.data$deaths_lwr_F,0),
        deaths_upr = dplyr::coalesce(.data$deaths_upr_M,0) + dplyr::coalesce(.data$deaths_upr_F,0)
      ) %>%
      dplyr::select(period, deaths_hat, deaths_lwr, deaths_upr) %>%
      dplyr::arrange(period)
  }
  
  # 2) Observado total
  obs_tot <- res_both$combined$obs_annual
  if (is.null(obs_tot) || !all(c("period","obs") %in% names(obs_tot))) {
    obs_m <- if (!is.null(res_both$resM$obs_annual)) res_both$resM$obs_annual else NULL
    obs_f <- if (!is.null(res_both$resF$obs_annual)) res_both$resF$obs_annual else NULL
    if (!is.null(obs_m) && !is.null(obs_f)) {
      obs_tot <- dplyr::full_join(obs_m, obs_f, by="period") %>%
        dplyr::mutate(obs = dplyr::coalesce(.data$obs.x,0) + dplyr::coalesce(.data$obs.y,0)) %>%
        dplyr::select(period, obs)
    } else if (!is.null(obs_m)) {
      obs_tot <- dplyr::select(obs_m, period, obs)
    } else if (!is.null(obs_f)) {
      obs_tot <- dplyr::select(obs_f, period, obs)
    } else {
      obs_tot <- tibble::tibble(period=integer(0), obs=double(0))
    }
  }
  
  # 3) Year bisagra (último histórico)
  lh <- res_both$combined$last_hist_year %||% NA_integer_
  
  max_year <- projection_common_max_year_from_res_both(res_both, policy = "endogenous_max", default = max(tot$period, na.rm = TRUE))
  if (is.finite(max_year)) {
    tot <- clip_to_year(tot, max_year = max_year)
  }

  # 4) Plot
  g <- ggplot2::ggplot(tot, ggplot2::aes(x = period)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = deaths_lwr, ymax = deaths_upr), alpha = 0.15) +
    ggplot2::geom_line(ggplot2::aes(y = deaths_hat), linewidth = 1) +
    ggplot2::labs(x = "Year", y = "Deaths", title = title,
                  subtitle = if (is.finite(max_year) && is.finite(lh)) paste0("Observed through ", lh, " | Endogenous projection window: ", lh + 1, "–", max_year) else NULL)
  
  if (show_hist_points && nrow(obs_tot) > 0) {
    g <- g + ggplot2::geom_point(data = obs_tot, ggplot2::aes(y = obs), size = 1.2)
  }
  if (is.finite(lh)) {
    g <- g + ggplot2::geom_vline(xintercept = lh, linetype = 3)
  }
  
  g + ggplot2::scale_y_continuous(limits = c(0, NA), expand = c(0, 0)) +
    ggplot2::theme_minimal(base_size = 12)
}



plot_incidence_proj_dual <- function(res_sex,
                                     title_base = "Lung cancer incidence: observed and projected") {
  stopifnot(!is.null(res_sex$inc_annual_cond),
            !is.null(res_sex$inc_obs_annual))
  inc_bapc <- tryCatch(res_sex$inc_annual_bapc, error = function(e) NULL) %||% tryCatch(res_sex$inc_annual_noP, error = function(e) NULL)
  stopifnot(!is.null(inc_bapc))

  sex_tag   <- ifelse(res_sex$sex == "M", "Males", "Females")
  last_year <- res_sex$diag$last_hist_year
  obs_I     <- res_sex$inc_obs_annual
  max_year <- projection_max_year_from_res_sex(res_sex, policy = "endogenous_max", default = max(c(last_year, PROJ_TO), na.rm = TRUE))
  end_credible <- projection_max_year_from_res_sex(res_sex, policy = "credible", default = max_year)
  end_caution  <- projection_max_year_from_res_sex(res_sex, policy = "caution", default = max_year)

  proj_bapc <- inc_bapc %>%
    dplyr::filter(period >= last_year + 1) %>%
    clip_to_year(max_year = max_year) %>%
    dplyr::mutate(model = "I")
  proj_cond <- res_sex$inc_annual_cond %>%
    dplyr::filter(period >= last_year + 1) %>%
    clip_to_year(max_year = max_year) %>%
    dplyr::mutate(model = "I | P")

  ggplot2::ggplot() +
    ggplot2::geom_ribbon(data = proj_bapc, ggplot2::aes(period, ymin = cases_lwr, ymax = cases_upr),
                         fill = "#ff7f0e", alpha = .15, show.legend = FALSE) +
    ggplot2::geom_ribbon(data = proj_cond, ggplot2::aes(period, ymin = cases_lwr, ymax = cases_upr),
                         fill = "#2ca02c", alpha = .18, show.legend = FALSE) +
    ggplot2::geom_line(data = proj_cond, ggplot2::aes(period, cases_hat, linetype = model, color = model), linewidth = 1) +
    ggplot2::geom_line(data = proj_bapc, ggplot2::aes(period, cases_hat, linetype = model, color = model), linewidth = 1) +
    ggplot2::geom_line(data = obs_I, ggplot2::aes(period, obs), color = "black", linewidth = 1) +
    ggplot2::geom_vline(xintercept = last_year + 0.5, linetype = 3) +
    ggplot2::scale_color_manual(values = c("I" = "#ff7f0e", "I | P" = "#2ca02c")) +
    ggplot2::scale_linetype_manual(values = c("I" = "solid", "I | P" = "solid")) +
    ggplot2::labs(x = "Year", y = "Cases", color = "Projection", linetype = "Projection",
                  title = paste0(title_base, " — ", sex_tag),
                  subtitle = paste0("Observed through ", last_year,
                                    " | Endogenous projection window: ", last_year + 1, "–", max_year,
                                    " | Credible through ", end_credible,
                                    " | Caution through ", end_caution),
                  caption  = "Bands: 95% CI (transformed INLA quantiles)") +
    ggplot2::scale_y_continuous(limits = c(0, NA), expand = c(0, 0)) +
    ggplot2::theme_minimal(base_size = 12) + ggplot2::theme(legend.position = "bottom")
}

plot_incidence_counterfactual_noP <- function(res_sex,
                                              title_base = "Lung cancer incidence counterfactual: no prevalence channel") {
  stopifnot(!is.null(res_sex$inc_annual_noP), !is.null(res_sex$inc_obs_annual))
  sex_tag   <- ifelse(res_sex$sex == "M", "Males", "Females")
  last_year <- res_sex$diag$last_hist_year
  obs_I     <- res_sex$inc_obs_annual
  max_year <- projection_max_year_from_res_sex(res_sex, policy = "endogenous_max", default = max(c(last_year, PROJ_TO), na.rm = TRUE))

  proj_noP <- res_sex$inc_annual_noP %>%
    dplyr::filter(period >= last_year + 1) %>%
    clip_to_year(max_year = max_year)

  ggplot2::ggplot() +
    ggplot2::geom_ribbon(data = proj_noP, ggplot2::aes(period, ymin = cases_lwr, ymax = cases_upr),
                         fill = "#ff7f0e", alpha = .15, show.legend = FALSE) +
    ggplot2::geom_line(data = proj_noP, ggplot2::aes(period, cases_hat), color = "#ff7f0e", linewidth = 1) +
    ggplot2::geom_line(data = obs_I, ggplot2::aes(period, obs), color = "black", linewidth = 1) +
    ggplot2::geom_vline(xintercept = last_year + 0.5, linetype = 3) +
    ggplot2::labs(x = "Year", y = "Cases",
                  title = paste0(title_base, " — ", sex_tag),
                  subtitle = paste0("Observed through ", last_year,
                                    " | Endogenous projection window: ", last_year + 1, "–", max_year),
                  caption = "Counterfactual series with prevalence channel removed ex post") +
    ggplot2::scale_y_continuous(limits = c(0, NA), expand = c(0, 0)) +
    ggplot2::theme_minimal(base_size = 12)
}

# === Etiquetas bonitas para la leyenda ===
.scn_label <- function(x) dplyr::recode(x,
                                        "down1pc" = "Down 1% per year",
                                        "freeze"  = "Frozen at 2022",
                                        "down3pc" = "Down 3% per year",
                                        "quit"    = "Full cessation",
                                        .default  = x
)

# === TOTAL (todas las causas) por escenarios de prevalencia ===
plot_total_mort_scenarios <- function(total_by_scn,
                                      obs_total,
                                      last_hist_year = PERIOD_M_MAX,
                                      title = "Mortality — Smoking prevalence scenarios") {
  stopifnot(all(c("scenario","period","deaths_hat","deaths_lwr","deaths_upr") %in% names(total_by_scn)))
  # Etiquetas y colores
  df <- total_by_scn %>%
    dplyr::mutate(
      scenario = factor(scenario, levels = c("freeze","down1pc","down3pc","quit")),
      esc = dplyr::recode(as.character(scenario),
                          "down1pc"="Down 1% per year",
                          "freeze" ="Frozen at 2022",
                          "down3pc"="Down 3% per year",
                          "quit"   ="Full cessation")
    ) %>%
    dplyr::arrange(scenario, period)
  
  cols <- c("Frozen at 2022"="#d62728", # rojo
            "Down 1% per year"    ="#ff7f0e", # orange
            "Down 3% per year"    ="#2ca02c", # green
            "Full cessation"="#1f77b4") # azul
  
  g <- ggplot2::ggplot(df, ggplot2::aes(x = period, group = esc)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = deaths_lwr, ymax = deaths_upr, fill = esc), alpha = 0.15) +
    ggplot2::geom_line(ggplot2::aes(y = deaths_hat, color = esc), linewidth = 1) +
    ggplot2::scale_color_manual(values = cols, name = "Smoking prevalence scenario") +
    ggplot2::scale_fill_manual(values = cols,  name = "Smoking prevalence scenario") +
    ggplot2::scale_y_continuous(labels = scales::label_comma(), limits = c(0, NA), expand = c(0, 0)) +
    ggplot2::labs(x = "Year", y = "Deaths", title = title) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom")
  
  # histórico en negro
  if (!missing(obs_total) && !is.null(obs_total) && nrow(obs_total) > 0) {
    g <- g + ggplot2::geom_line(data = obs_total, ggplot2::aes(x = period, y = obs),
                                inherit.aes = FALSE, color = "black", linewidth = 1)
  }
  # línea de corte
  if (is.finite(last_hist_year)) {
    g <- g + ggplot2::geom_vline(xintercept = last_hist_year + 0.5, linetype = 3)
  }
  g
}




# =============================================================

# =============================================================
# Additional plotting helpers extracted from the multi-cause section
# =============================================================
# ---------------------------------------------------------------------
# Plot TOTAL (both sexes aggregated) for one cause — robust to both call styles
# ---------------------------------------------------------------------

build_scn_list_from_csv_one_cause <- function(cause_id, agg_csv, plots_dir,
                                              baseline_year = 2022) {
  
  if (!file.exists(agg_csv)) stop("No encuentro agg_csv: ", agg_csv)
  
  # plots_dir = .../Resultados/cause_<id>/plots  -> root = .../Resultados/cause_<id>
  cause_root <- dirname(plots_dir)
  
  rds_file <- file.path(cause_root, paste0(cause_id, "_res_both_freeze.rds"))
  if (!file.exists(rds_file)) {
    # fallback por si el freeze se guardó con otro nombre
    rds_alt <- file.path(cause_root, paste0(cause_id, "_res_both.rds"))
    if (file.exists(rds_alt)) rds_file <- rds_alt
  }
  if (!file.exists(rds_file)) stop("No encuentro RDS para observados: ", rds_file)
  
  rb <- readRDS(rds_file)
  
  # Observed (black line)
  obs <- rb$combined$obs_annual %>%
    dplyr::transmute(
      year = period,
      mean = deaths,
      lwr  = deaths,
      upr  = deaths,
      scenario = "Observed"
    )
  
  # Scenario projections from the aggregated CSV
  by_cause_raw <- readr::read_csv(agg_csv, show_col_types = FALSE) %>%
    dplyr::filter(cause_id == !!cause_id)
  scn_codes <- unique(as.character(by_cause_raw$scenario))
  scn_levels <- unname(stats::na.omit(scenario_labels_en[scn_codes]))
  by_cause <- by_cause_raw %>%
    dplyr::mutate(
      esc = unname(scenario_labels_en[as.character(scenario)]),
      esc = factor(esc, levels = scn_levels)
    ) %>%
    dplyr::transmute(
      year = period,
      mean = deaths_hat,
      lwr  = deaths_lwr,
      upr  = deaths_upr,
      scenario = as.character(esc)
    ) %>%
    dplyr::filter(year >= baseline_year)
  
  # split scenarios + add observed
  scn_list <- split(by_cause, by_cause$scenario)
  scn_list[["Observed"]] <- obs
  scn_list
}


plot_and_save_scenarios_one_cause <- function(
    cause_id, cause_label, out_dir, agg_csv,
    lang = c("en","es"),
    file_stub = cause_id
) {
  lang <- match.arg(lang)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  agg <- readr::read_csv(agg_csv, show_col_types = FALSE)
  stopifnot(all(c("cause_id","scenario","period","deaths_hat","deaths_lwr","deaths_upr") %in% names(agg)))
  
  if (lang == "en") {
    xlab <- "Year"
    ylab <- "Deaths"
    leg  <- "Smoking prevalence scenario"
    esc_map    <- scenario_labels_en
    title_txt  <- cause_label
    out_png <- file.path(out_dir, sprintf("fig_%s_mort_scen_TOTAL.png", file_stub))
  } else {
    xlab <- "Year"
    ylab <- "Deaths"
    leg  <- "Smoking prevalence scenario"
    esc_map    <- scenario_labels
    title_txt  <- cause_label
    out_png <- file.path(out_dir, sprintf("fig_%s_mort_scen_TOTAL_ES.png", file_stub))
  }
  
  df <- agg %>%
    dplyr::filter(cause_id == !!cause_id) %>%
    dplyr::mutate(
      scenario = as.character(scenario),
      esc      = unname(esc_map[scenario])
    ) %>%
    dplyr::filter(!is.na(esc)) %>%
    dplyr::mutate(esc = factor(esc, levels = unique(unname(esc_map[unique(scenario)])))) %>%
    dplyr::arrange(scenario, period)
  esc_levels <- levels(df$esc)
  
  # Corte e histórico desde el RDS "freeze" de la causa
  dir_cause  <- file.path(BASE_RESULTS_DIR, paste0("cause_", cause_id))
  rds_freeze <- file.path(dir_cause, sprintf("%s_res_both_freeze.rds", cause_id))
  if (!file.exists(rds_freeze)) {
    stop("No encuentro RDS freeze para la causa ", cause_id, ": ", rds_freeze)
  }
  rb <- readRDS(rds_freeze)
  
  last_hist_year <- rb$combined$last_hist_year
  obs_hist <- rb$combined$obs_annual %>%
    dplyr::select(period, obs) %>%
    dplyr::arrange(period)
  
  df_future <- df %>% dplyr::filter(period > last_hist_year)
  
  col_map <- if (lang == "en") scenario_colors_en else scenario_colors
  
  p <- ggplot2::ggplot(df_future, ggplot2::aes(x = period, group = esc)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = deaths_lwr, ymax = deaths_upr, fill = esc), alpha = 0.15) +
    ggplot2::geom_line(ggplot2::aes(y = deaths_hat, color = esc), linewidth = 1) +
    ggplot2::geom_line(data = obs_hist, ggplot2::aes(x = period, y = obs),
                       inherit.aes = FALSE, color = "black", linewidth = 1) +
    ggplot2::geom_vline(xintercept = last_hist_year + 0.5, linetype = 3) +
    ggplot2::scale_color_manual(values = col_map, breaks = esc_levels, limits = esc_levels, name = leg) +
    ggplot2::scale_fill_manual(values  = col_map, breaks = esc_levels, limits = esc_levels, name = leg) +
    ggplot2::labs(x = xlab, y = ylab, title = title_txt) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::coord_cartesian(ylim = c(0, NA)) +
    ggplot2::theme_minimal(base_size = 16) +
    ggplot2::theme(
      legend.position = "bottom",
      plot.title  = ggplot2::element_text(hjust = 0.5, size = 20),
      axis.title  = ggplot2::element_text(size = 18),
      axis.text   = ggplot2::element_text(size = 14),
      legend.text = ggplot2::element_text(size = 13),
      legend.title= ggplot2::element_text(size = 14)
    )
  
  ggplot2::ggsave(out_png, p, width = 12, height = 7, bg = "white", dpi = 160)
  message("✔ Scenarios TOTAL (M+F) saved: ", out_png)
}



plot_and_save_scenarios_one_cause_bysex <- function(
    cause_id,
    out_dir,
    scenarios = scenario_levels,
    out_dir_results = BASE_RESULTS_DIR
) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Freeze RDS: provides historical series + last observed year
  rb_freeze <- get_res_both_from_rds(cause_id, scenario = "freeze", out_dir = out_dir_results, strict = TRUE)
  
  for (sx in c("M","F")) {
    
    rfreeze_sex <- if (sx == "M") rb_freeze$resM else rb_freeze$resF
    if (is.null(rfreeze_sex) || is.null(rfreeze_sex$obs_annual) || nrow(rfreeze_sex$obs_annual) == 0) {
      message("⚠ No historical series for ", cause_id, " sex=", sx, " (skipping).")
      next
    }
    
    last_hist_year <- rfreeze_sex$diag$last_hist_year %||% rb_freeze$combined$last_hist_year
    obs_hist <- rfreeze_sex$obs_annual %>%
      dplyr::select(period, obs) %>%
      dplyr::arrange(period)
    
    # Collect projections for each scenario from RDS
    proj_list <- list()
    for (scn in scenarios) {
      rb_scn <- get_res_both_from_rds(cause_id, scenario = scn, out_dir = out_dir_results, strict = FALSE)
      if (is.null(rb_scn)) next
      
      r_sx <- if (sx == "M") rb_scn$resM else rb_scn$resF
      if (is.null(r_sx) || is.null(r_sx$annual_anchor) || nrow(r_sx$annual_anchor) == 0) next
      
      proj_list[[length(proj_list) + 1]] <- r_sx$annual_anchor %>%
        dplyr::mutate(scenario = scn)
    }
    
    if (!length(proj_list)) {
      message("⚠ No projections found for ", cause_id, " sex=", sx, " (skipping).")
      next
    }
    
    esc_levels <- unname(scenario_labels_en[scenarios])
    df <- dplyr::bind_rows(proj_list) %>%
      dplyr::mutate(
        scenario = factor(scenario, levels = scenarios),
        esc      = factor(unname(scenario_labels_en[as.character(scenario)]), levels = esc_levels)
      ) %>%
      dplyr::arrange(scenario, period)
    
    df_future <- df %>% dplyr::filter(period > last_hist_year)
    
    p <- ggplot2::ggplot(df_future, ggplot2::aes(x = period, group = esc)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = deaths_lwr, ymax = deaths_upr, fill = esc), alpha = 0.15) +
      ggplot2::geom_line(ggplot2::aes(y = deaths_hat, color = esc), linewidth = 1) +
      ggplot2::geom_line(
        data = obs_hist,
        ggplot2::aes(x = period, y = obs),
        inherit.aes = FALSE, color = "black", linewidth = 1
      ) +
      ggplot2::geom_vline(xintercept = last_hist_year + 0.5, linetype = 3) +
      ggplot2::scale_color_manual(
        values = scenario_colors_en, breaks = esc_levels, limits = esc_levels,
        name = "Smoking prevalence scenario"
      ) +
      ggplot2::scale_fill_manual(
        values = scenario_colors_en, breaks = esc_levels, limits = esc_levels,
        name = "Smoking prevalence scenario"
      ) +
      ggplot2::labs(
        x = "Year", y = "Deaths",
        title = paste0(get_cause_label_en(cause_id), " — ", sex_labels_en[[sx]])
      ) +
      ggplot2::scale_y_continuous(expand = c(0, 0)) +
      ggplot2::coord_cartesian(ylim = c(0, NA)) +
      ggplot2::theme_minimal(base_size = 16) +
      ggplot2::theme(
        legend.position = "bottom",
        plot.title  = ggplot2::element_text(hjust = 0.5, size = 20),
        axis.title  = ggplot2::element_text(size = 18),
        axis.text   = ggplot2::element_text(size = 14),
        legend.text = ggplot2::element_text(size = 13),
        legend.title= ggplot2::element_text(size = 14)
      )
    
    out_png <- file.path(out_dir, sprintf("fig_%s_mort_scen_%s.png", cause_id, sx))
    ggplot2::ggsave(out_png, p, width = 12, height = 7, bg = "white", dpi = 160)
    message("✔ Scenarios by sex saved: ", out_png)
  }
  
  invisible(TRUE)
}



save_total_freeze_plots <- function(cause_ids = causes$cause_id,
                                    scenario = "freeze",
                                    out_dir = PLOTS_TOTAL_DIR,
                                    out_dir_results = BASE_RESULTS_DIR) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  tot_I <- build_total_incidence_from_rds(cause_ids = cause_ids, scenario = scenario, out_dir = out_dir_results)
  tot_M <- build_total_mortality_from_rds(cause_ids = cause_ids, scenario = scenario, out_dir = out_dir_results)

  pI <- plot_total_series(tot_I, "Incidence — sum across causes", paste0("scenario: ", scenario), "Cases")
  pM <- plot_total_series(tot_M, "Mortality — sum across causes", paste0("scenario: ", scenario), "Deaths")

  ggplot2::ggsave(file.path(out_dir, sprintf("total_incidence_%s.png", scenario)), pI, width = 10, height = 6, bg = "white", dpi = 160)
  ggplot2::ggsave(file.path(out_dir, sprintf("total_mortality_%s.png", scenario)), pM, width = 10, height = 6, bg = "white", dpi = 160)

  invisible(list(total_incidence = tot_I, total_mortality = tot_M, p_incidence = pI, p_mortality = pM))
}
# =============================================================
# Internal / diagnostic plots
# =============================================================
# First cut note:
# In this initial split, most plotting code with stable aesthetics was moved
# to R/07_plots_paper.R. Internal and lightweight diagnostic plots are still
# mixed into dev/legacy_debug_and_scenarios.R and will be separated in the
# next refactor pass.

# guarda los APC disponibles para un sexo si el objeto existe
save_apc_plots_for_sex <- function(res_sex, out_dir, label = NULL, sex_tag = "", overwrite = FALSE) {
  if (is.null(res_sex)) return(invisible(NULL))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  cid <- res_sex$cause_id
  if (is.null(cid) || !length(cid) || !nzchar(cid)) {
    cid <- attr(res_sex, "cause_id", exact = TRUE)
    if (is.null(cid) || !length(cid) || !nzchar(cid)) {
      cid <- get0("CAUSE_ID", ifnotfound = "cause", inherits = TRUE)
    }
  }

  g_apc_inc <- try(plot_apc_inc(res_sex), silent = TRUE)
  f_inc <- file.path(out_dir, sprintf("apc_inc_%s_%s.png", cid, sex_tag))
  if ((!file.exists(f_inc) || isTRUE(overwrite)) && !inherits(g_apc_inc, "try-error") && !is.null(g_apc_inc)) {
    ggplot2::ggsave(
      f_inc,
      g_apc_inc,
      width = 9.5, height = 7.5, dpi = 200, bg = "white"
    )
  }

  g_apc_mort <- try(plot_apc_mort(res_sex), silent = TRUE)
  f_mort <- file.path(out_dir, sprintf("apc_mort_%s_%s.png", cid, sex_tag))
  if ((!file.exists(f_mort) || isTRUE(overwrite)) && !inherits(g_apc_mort, "try-error") && !is.null(g_apc_mort)) {
    ggplot2::ggsave(
      f_mort,
      g_apc_mort,
      width = 9.5, height = 7.5, dpi = 200, bg = "white"
    )
  }

  invisible(NULL)
}
