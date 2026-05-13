# Common helpers for rebuilding prediction details consistently
safe_offset_numeric <- function(x) {
  dplyr::coalesce(as.numeric(x), 0)
}

build_bapc_prediction_detail <- function(df_all, lp_df, sex_sel) {
  stopifnot(is.data.frame(df_all), is.data.frame(lp_df), nrow(df_all) == nrow(lp_df))
  df_all %>%
    dplyr::transmute(
      sex = sex_sel,
      period,
      hist_flag = !is.na(y),
      mu_hat = exp(as.numeric(lp_df$mean)),
      mu_lwr = pmax(0, exp(as.numeric(lp_df$`0.025quant`))),
      mu_upr = pmax(0, exp(as.numeric(lp_df$`0.975quant`)))
    )
}

build_anchor_prediction_detail <- function(data_df, lp_df, sex_sel) {
  stopifnot(is.data.frame(data_df), is.data.frame(lp_df), nrow(data_df) == nrow(lp_df))
  data_df %>%
    dplyr::transmute(
      sex = sex_sel,
      period,
      hist_flag = !is.na(y),
      log_mort_ext = as.numeric(log_mort_ext),
      mort_anchor_tech_offset = safe_offset_numeric(mort_anchor_tech_offset),
      offset_total = log_mort_ext + mort_anchor_tech_offset,
      eta_lp_hat = as.numeric(lp_df$mean),
      eta_lp_lwr = as.numeric(lp_df$`0.025quant`),
      eta_lp_upr = as.numeric(lp_df$`0.975quant`),
      # FIX:
      # tratar SIEMPRE el predictor lineal como predictor total y recuperar
      # el residual restando offset_total tanto en histórico como en futuro.
      eta_resid_hat = eta_lp_hat - offset_total,
      eta_resid_lwr = eta_lp_lwr - offset_total,
      eta_resid_upr = eta_lp_upr - offset_total,
      mu_hat_legacy = exp(eta_lp_hat),
      mu_hat = exp(eta_resid_hat + offset_total),
      mu_lwr = pmax(0, exp(eta_resid_lwr + offset_total)),
      mu_upr = pmax(0, exp(eta_resid_upr + offset_total)),
      mort_ext_deaths = dplyr::coalesce(as.numeric(mort_ext_deaths), 0)
    )
}

summarise_annual_anchor_components <- function(pred_df) {
  if (!is.data.frame(pred_df) || !nrow(pred_df)) {
    return(tibble::tibble(
      sex = character(), period = numeric(), n_cells = integer(),
      deaths_ext = numeric(), log_mort_ext_mean = numeric(),
      offset_total_mean = numeric(), eta_resid_hat_mean = numeric(),
      mu_hat_sum = numeric()
    ))
  }
  pred_df %>%
    dplyr::mutate(.w_ext = pmax(dplyr::coalesce(as.numeric(mort_ext_deaths), 0), 1e-12)) %>%
    dplyr::group_by(sex, period) %>%
    dplyr::summarise(
      n_cells = dplyr::n(),
      deaths_ext = sum(mort_ext_deaths, na.rm = TRUE),
      log_mort_ext_mean = stats::weighted.mean(log_mort_ext, w = .w_ext, na.rm = TRUE),
      offset_total_mean = stats::weighted.mean(offset_total, w = .w_ext, na.rm = TRUE),
      eta_resid_hat_mean = mean(eta_resid_hat, na.rm = TRUE),
      mu_hat_sum = sum(mu_hat, na.rm = TRUE),
      .groups = "drop"
    )
}
