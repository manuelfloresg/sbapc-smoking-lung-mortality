
# Common helpers for engine consistency / annual aggregation
engine_require_stock_former_policy <- function(beta_mode,
                                               prev_cfg,
                                               prev_inc_channel_mode = PREV_INC_CHANNEL_MODE) {
  ok <- identical(beta_mode, "fixed_rr_offset") &&
    identical(prev_inc_channel_mode, "stock_former") &&
    identical(prev_cfg$axis, "period")
  if (!isTRUE(ok)) {
    stop(
      paste0(
        "Motor policy violation: la ruta operativa exige ",
        "beta_mode='fixed_rr_offset', PREV_INC_CHANNEL_MODE='stock_former', ",
        "PREV_INC_KEEP_LEGACY_AI=FALSE y prev_cfg$axis='period'."
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

summarise_annual_prediction <- function(pred_df) {
  if (!is.data.frame(pred_df) || !nrow(pred_df)) {
    return(tibble::tibble(
      sex = character(), period = numeric(),
      deaths_hat = numeric(), deaths_lwr = numeric(), deaths_upr = numeric()
    ))
  }
  pred_df %>%
    dplyr::group_by(sex, period) %>%
    dplyr::summarise(
      deaths_hat = sum(mu_hat, na.rm = TRUE),
      deaths_lwr = sum(mu_lwr, na.rm = TRUE),
      deaths_upr = sum(mu_upr, na.rm = TRUE),
      .groups = "drop"
    )
}

summarise_annual_external <- function(pred_df) {
  if (!is.data.frame(pred_df) || !nrow(pred_df)) {
    return(tibble::tibble(
      sex = character(), period = numeric(),
      deaths_ext = numeric(), log_mort_ext_mean = numeric(), offset_total_mean = numeric()
    ))
  }
  pred_df %>%
    dplyr::mutate(.w_ext = pmax(dplyr::coalesce(as.numeric(mort_ext_deaths), 0), 1e-12)) %>%
    dplyr::group_by(sex, period) %>%
    dplyr::summarise(
      deaths_ext = sum(mort_ext_deaths, na.rm = TRUE),
      log_mort_ext_mean = stats::weighted.mean(log_mort_ext, w = .w_ext, na.rm = TRUE),
      offset_total_mean = stats::weighted.mean(offset_total, w = .w_ext, na.rm = TRUE),
      .groups = "drop"
    )
}
