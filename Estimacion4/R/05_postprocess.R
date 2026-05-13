# Safe serializer: muffles the known non-fatal warning produced by namespace refs in INLA objects
safe_saveRDS <- function(object, file, ...) {
  withCallingHandlers(
    saveRDS(object, file = file, ...),
    warning = function(w) {
      msg <- conditionMessage(w)
      if (grepl("package:stats.*proceso de carga|package:stats.*loading process", msg, ignore.case = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )
  invisible(file)
}

.save_if_missing_plot <- function(file, expr, overwrite = FALSE) {
  if (!file.exists(file) || isTRUE(overwrite)) {
    force(expr)
  }
  invisible(file)
}

# =============================================================
# 9) Helpers multi-causa: guardado estándar + empaquetado
# =============================================================

save_all_outputs <- function(res_both,
                             cause_id = NULL,
                             label = NULL,
                             out_base = BASE_RESULTS_DIR,
                             static_out_base = NULL,
                             overwrite_static = FALSE,
                             flatten_single_cause = FALSE) {
  .write_diag_csv <- function(x, file) {
    if (is.null(x)) return(invisible(FALSE))
    if (!is.data.frame(x) || !nrow(x)) return(invisible(FALSE))
    readr::write_csv(tibble::as_tibble(x), file)
    invisible(TRUE)
  }

  .export_external_kernel_diags <- function(res_sex, sex_tag, out_dir) {
    if (is.null(res_sex)) return(invisible(NULL))
    diag <- tryCatch(res_sex$diag, error = function(e) NULL)
    if (is.null(diag)) return(invisible(NULL))

    .write_diag_csv(tryCatch(diag$mort_kernel, error = function(e) NULL),
                    file.path(out_dir, sprintf("mort_kernel_%s.csv", sex_tag)))
    .write_diag_csv(tryCatch(diag$mort_kernel_summary_cond, error = function(e) NULL),
                    file.path(out_dir, sprintf("mort_kernel_summary_cond_%s.csv", sex_tag)))

    invisible(NULL)
  }
  # fallbacks para cause_id/label
  if (is.null(cause_id) || !length(cause_id) || !nzchar(cause_id)) {
    cause_id <- attr(res_both, "cause_id", exact = TRUE)
    if (is.null(cause_id) || !length(cause_id) || !nzchar(cause_id)) {
      cause_id <- get0("CAUSE_ID", ifnotfound = "cause", inherits = TRUE)
    }
  }
  if (is.null(label) || !length(label) || !nzchar(label)) {
    label <- attr(res_both, "label", exact = TRUE)
    if (is.null(label) || !length(label) || !nzchar(label)) label <- cause_id
  }

  scenario_dir <- if (isTRUE(flatten_single_cause)) out_base else file.path(out_base, paste0("cause_", cause_id))
  scenario_plots_dir <- file.path(scenario_dir, "plots")
  scenario_diag_dir <- file.path(scenario_dir, "diagnostics")
  dir.create(scenario_plots_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(scenario_diag_dir, recursive = TRUE, showWarnings = FALSE)

  if (is.null(static_out_base)) static_out_base <- out_base
  static_dir <- if (isTRUE(flatten_single_cause)) static_out_base else file.path(static_out_base, paste0("cause_", cause_id))
  static_plots_dir <- file.path(static_dir, "plots")
  dir.create(static_plots_dir, recursive = TRUE, showWarnings = FALSE)

  # Full result object (scenario-specific because projections differ)
  safe_saveRDS(res_both, file.path(scenario_dir, sprintf("%s_res_both.rds", cause_id)))

  # --- Scenario-specific diagnostics useful for downstream tables/figures ---
  .export_external_kernel_diags(res_both$resM, "M", scenario_diag_dir)
  .export_external_kernel_diags(res_both$resF, "F", scenario_diag_dir)

  # --- Static outputs shared across scenarios but still cause-specific ---
  calib_file <- file.path(static_plots_dir, sprintf("%s_Astar_calibration_bothsex.png", cause_id))
  if (!file.exists(calib_file) || isTRUE(overwrite_static)) {
    try(save_ai_calibration_png(res_both, file = calib_file, width = 10, height = 4, dpi = 300), silent = TRUE)
  }

  if (!is.null(res_both$resM)) save_apc_plots_for_sex(res_both$resM, file.path(static_plots_dir, "APC_M"), label, "Males", overwrite = overwrite_static)
  if (!is.null(res_both$resF)) save_apc_plots_for_sex(res_both$resF, file.path(static_plots_dir, "APC_F"), label, "Females", overwrite = overwrite_static)

  # --- Scenario-dependent outputs ---
  if (!is.null(res_both$resM)) {
    ggplot2::ggsave(file.path(scenario_plots_dir, sprintf("%s_mort_M.png", cause_id)),
                    plot_projections_mort(res_both$resM, paste0("Cancer mortality for ", label)),
                    width = 9, height = 5.2, bg = "white", dpi = 200)
    ggplot2::ggsave(file.path(scenario_plots_dir, sprintf("%s_inc_M.png", cause_id)),
                    plot_incidence_proj_dual(res_both$resM, paste0("Cancer incidence for ", label)),
                    width = 9, height = 5.2, bg = "white", dpi = 200)
  }
  if (!is.null(res_both$resF)) {
    ggplot2::ggsave(file.path(scenario_plots_dir, sprintf("%s_mort_F.png", cause_id)),
                    plot_projections_mort(res_both$resF, paste0("Cancer mortality for ", label)),
                    width = 9, height = 5.2, bg = "white", dpi = 200)
    ggplot2::ggsave(file.path(scenario_plots_dir, sprintf("%s_inc_F.png", cause_id)),
                    plot_incidence_proj_dual(res_both$resF, paste0("Cancer incidence for ", label)),
                    width = 9, height = 5.2, bg = "white", dpi = 200)
  }

  ggplot2::ggsave(file.path(scenario_plots_dir, sprintf("%s_mort_total.png", cause_id)),
                  plot_projections_mort_total(res_both, paste0("Cancer mortality for ", label, " — Total")),
                  width = 9, height = 5.2, bg = "white", dpi = 200)

  invisible(list(
    scenario_dir = scenario_dir,
    scenario_plots_dir = scenario_plots_dir,
    scenario_diag_dir = scenario_diag_dir,
    static_dir = static_dir,
    static_plots_dir = static_plots_dir
  ))
}


pack_params <- function(res_both, cause_id, label) {
  out <- list()
  if (!is.null(res_both$resM)) out[[length(out)+1]] <- tibble::tibble(
    cause_id = cause_id, label = label, sex = "M",
    A_I_star   = res_both$resM$diag$A_I_star,
    prev_sign  = res_both$resM$diag$prev_sign,
    beta_mode  = tryCatch(res_both$resM$diag$prev$beta_mode, error = function(e) BETA_MODE),
    rr_inc     = tryCatch(res_both$resM$diag$prev$rr_inc, error = function(e) NA_real_),
    projection_end_credible = tryCatch(projection_max_year_from_frontier(res_both$resM$diag$projection_horizon_frontier, policy = "credible"), error = function(e) NA_integer_),
    projection_end_caution  = tryCatch(projection_max_year_from_frontier(res_both$resM$diag$projection_horizon_frontier, policy = "caution"), error = function(e) NA_integer_),
    projection_end_risky    = tryCatch(projection_max_year_from_frontier(res_both$resM$diag$projection_horizon_frontier, policy = "risky"), error = function(e) NA_integer_),
    projection_end_max      = tryCatch(projection_max_year_from_frontier(res_both$resM$diag$projection_horizon_frontier, policy = "endogenous_max"), error = function(e) NA_integer_),
    beta_I_mean= res_both$resM$params$beta_I_mean,
    beta_I_lwr = res_both$resM$params$beta_I_lwr,
    beta_I_upr = res_both$resM$params$beta_I_upr,
    L_I        = if (identical(tryCatch(res_both$resM$params$mort_link_mode, error = function(e) NA_character_), "external_kernel")) NA_integer_ else res_both$resM$params$L_I,
    bridge_years = res_both$resM$params$bridge_years,
    mort_link_mode = tryCatch(res_both$resM$params$mort_link_mode, error = function(e) NA_character_),
    mort_kernel_max_lag = tryCatch(res_both$resM$params$mort_kernel_max_lag, error = function(e) NA_integer_),
    mort_kernel_total_mass = tryCatch(res_both$resM$params$mort_kernel_total_mass, error = function(e) NA_real_),
    mort_bapc_trend_scenario = tryCatch(res_both$resM$diag$mort_bapc_trend_scenario, error = function(e) NA_character_),
    mort_bapc_future_mode = tryCatch(res_both$resM$diag$mort_bapc_future_mode, error = function(e) NA_character_)
  )
  if (!is.null(res_both$resF)) out[[length(out)+1]] <- tibble::tibble(
    cause_id = cause_id, label = label, sex = "F",
    A_I_star   = res_both$resF$diag$A_I_star,
    prev_sign  = res_both$resF$diag$prev_sign,
    beta_mode  = tryCatch(res_both$resF$diag$prev$beta_mode, error = function(e) BETA_MODE),
    rr_inc     = tryCatch(res_both$resF$diag$prev$rr_inc, error = function(e) NA_real_),
    projection_end_credible = tryCatch(projection_max_year_from_frontier(res_both$resF$diag$projection_horizon_frontier, policy = "credible"), error = function(e) NA_integer_),
    projection_end_caution  = tryCatch(projection_max_year_from_frontier(res_both$resF$diag$projection_horizon_frontier, policy = "caution"), error = function(e) NA_integer_),
    projection_end_risky    = tryCatch(projection_max_year_from_frontier(res_both$resF$diag$projection_horizon_frontier, policy = "risky"), error = function(e) NA_integer_),
    projection_end_max      = tryCatch(projection_max_year_from_frontier(res_both$resF$diag$projection_horizon_frontier, policy = "endogenous_max"), error = function(e) NA_integer_),
    beta_I_mean= res_both$resF$params$beta_I_mean,
    beta_I_lwr = res_both$resF$params$beta_I_lwr,
    beta_I_upr = res_both$resF$params$beta_I_upr,
    L_I        = if (identical(tryCatch(res_both$resF$params$mort_link_mode, error = function(e) NA_character_), "external_kernel")) NA_integer_ else res_both$resF$params$L_I,
    bridge_years = res_both$resF$params$bridge_years,
    mort_link_mode = tryCatch(res_both$resF$params$mort_link_mode, error = function(e) NA_character_),
    mort_kernel_max_lag = tryCatch(res_both$resF$params$mort_kernel_max_lag, error = function(e) NA_integer_),
    mort_kernel_total_mass = tryCatch(res_both$resF$params$mort_kernel_total_mass, error = function(e) NA_real_),
    mort_bapc_trend_scenario = tryCatch(res_both$resF$diag$mort_bapc_trend_scenario, error = function(e) NA_character_),
    mort_bapc_future_mode = tryCatch(res_both$resF$diag$mort_bapc_future_mode, error = function(e) NA_character_)
  )
  dplyr::bind_rows(out)
}


.common_horizon_row <- function(h_rows, cause_id, label) {
  h_rows <- tibble::as_tibble(h_rows)
  if (!nrow(h_rows)) {
    return(tibble::tibble(
      cause_id = cause_id,
      label = label,
      sex = "T",
      last_hist_year = NA_integer_,
      max_horizon_available = NA_real_,
      first_horizon_caution = NA_real_,
      first_horizon_risky = NA_real_,
      first_horizon_beyond_max = NA_real_,
      end_horizon_credible = NA_real_,
      end_horizon_caution = NA_real_,
      end_horizon_risky = NA_real_,
      end_year_credible = NA_real_,
      end_year_caution = NA_real_,
      end_year_risky = NA_real_,
      max_projection_year_endogenous = NA_integer_,
      support_floor_credible = NA_real_,
      support_floor_caution = NA_real_,
      support_floor_max = NA_real_,
      edge_share_credible = NA_real_,
      edge_share_caution = NA_real_,
      edge_share_max = NA_real_
    ))
  }

  first_non_na <- function(x) {
    x <- x[is.finite(x)]
    if (!length(x)) return(NA_real_)
    x[[1]]
  }
  min_finite <- function(x) {
    x <- x[is.finite(x)]
    if (!length(x)) return(NA_real_)
    min(x)
  }
  max_finite <- function(x) {
    x <- x[is.finite(x)]
    if (!length(x)) return(NA_real_)
    max(x)
  }

  tibble::tibble(
    cause_id = cause_id,
    label = label,
    sex = "T",
    last_hist_year = suppressWarnings(as.integer(max_finite(h_rows$last_hist_year))),
    max_horizon_available = min_finite(h_rows$max_horizon_available),
    first_horizon_caution = min_finite(h_rows$first_horizon_caution),
    first_horizon_risky = min_finite(h_rows$first_horizon_risky),
    first_horizon_beyond_max = min_finite(h_rows$first_horizon_beyond_max),
    end_horizon_credible = min_finite(h_rows$end_horizon_credible),
    end_horizon_caution = min_finite(h_rows$end_horizon_caution),
    end_horizon_risky = min_finite(h_rows$end_horizon_risky),
    end_year_credible = min_finite(h_rows$end_year_credible),
    end_year_caution = min_finite(h_rows$end_year_caution),
    end_year_risky = min_finite(h_rows$end_year_risky),
    max_projection_year_endogenous = suppressWarnings(as.integer(min_finite(h_rows$max_projection_year_endogenous))),
    support_floor_credible = first_non_na(h_rows$support_floor_credible),
    support_floor_caution = first_non_na(h_rows$support_floor_caution),
    support_floor_max = first_non_na(h_rows$support_floor_max),
    edge_share_credible = first_non_na(h_rows$edge_share_credible),
    edge_share_caution = first_non_na(h_rows$edge_share_caution),
    edge_share_max = first_non_na(h_rows$edge_share_max)
  )
}

pack_horizon <- function(res_both, cause_id, label) {
  out <- list()
  if (!is.null(res_both$resM)) {
    hM <- tryCatch(res_both$resM$diag$projection_horizon_frontier, error = function(e) NULL)
    if (is.data.frame(hM) && nrow(hM)) out[[length(out) + 1L]] <- hM %>% dplyr::mutate(cause_id = cause_id, label = label, .before = 1)
  }
  if (!is.null(res_both$resF)) {
    hF <- tryCatch(res_both$resF$diag$projection_horizon_frontier, error = function(e) NULL)
    if (is.data.frame(hF) && nrow(hF)) out[[length(out) + 1L]] <- hF %>% dplyr::mutate(cause_id = cause_id, label = label, .before = 1)
  }
  h_all <- dplyr::bind_rows(out)
  common_tbl <- .common_horizon_row(h_all, cause_id = cause_id, label = label)
  dplyr::bind_rows(h_all, common_tbl)
}

pack_proj <- function(res_both, cause_id, label) {
  out <- list()
  if (!is.null(res_both$resM)) {
    inc_bapc_M <- tryCatch(res_both$resM$inc_annual_bapc, error = function(e) NULL) %||% tryCatch(res_both$resM$inc_annual_noP, error = function(e) NULL)
    out[[length(out)+1]] <- dplyr::bind_rows(
      inc_bapc_M %>% dplyr::mutate(sex="M", cause_id, label, metric="incidence", series="I") %>%
        dplyr::rename(mean=cases_hat,lwr=cases_lwr,upr=cases_upr),
      res_both$resM$inc_annual_cond %>% dplyr::mutate(sex="M", cause_id, label, metric="incidence", series="I|P") %>%
        dplyr::rename(mean=cases_hat,lwr=cases_lwr,upr=cases_upr),
      res_both$resM$inc_annual_noP  %>% dplyr::mutate(sex="M", cause_id, label, metric="incidence", series="I_noP_cf") %>%
        dplyr::rename(mean=cases_hat,lwr=cases_lwr,upr=cases_upr),
      res_both$resM$annual_anchor     %>% dplyr::mutate(sex="M", cause_id, label, metric="mortality", series="M|I|P") %>%
        dplyr::rename(mean=deaths_hat,lwr=deaths_lwr,upr=deaths_upr),
      res_both$resM$annual_anchor_noP %>% dplyr::mutate(sex="M", cause_id, label, metric="mortality", series="M|I") %>%
        dplyr::rename(mean=deaths_hat,lwr=deaths_lwr,upr=deaths_upr)
    )
  }
  if (!is.null(res_both$resF)) {
    inc_bapc_F <- tryCatch(res_both$resF$inc_annual_bapc, error = function(e) NULL) %||% tryCatch(res_both$resF$inc_annual_noP, error = function(e) NULL)
    out[[length(out)+1]] <- dplyr::bind_rows(
      inc_bapc_F %>% dplyr::mutate(sex="F", cause_id, label, metric="incidence", series="I") %>%
        dplyr::rename(mean=cases_hat,lwr=cases_lwr,upr=cases_upr),
      res_both$resF$inc_annual_cond %>% dplyr::mutate(sex="F", cause_id, label, metric="incidence", series="I|P") %>%
        dplyr::rename(mean=cases_hat,lwr=cases_lwr,upr=cases_upr),
      res_both$resF$inc_annual_noP  %>% dplyr::mutate(sex="F", cause_id, label, metric="incidence", series="I_noP_cf") %>%
        dplyr::rename(mean=cases_hat,lwr=cases_lwr,upr=cases_upr),
      res_both$resF$annual_anchor     %>% dplyr::mutate(sex="F", cause_id, label, metric="mortality", series="M|I|P") %>%
        dplyr::rename(mean=deaths_hat,lwr=deaths_lwr,upr=deaths_upr),
      res_both$resF$annual_anchor_noP %>% dplyr::mutate(sex="F", cause_id, label, metric="mortality", series="M|I") %>%
        dplyr::rename(mean=deaths_hat,lwr=deaths_lwr,upr=deaths_upr)
    )
  }
  proj_tbl <- dplyr::bind_rows(out) %>% dplyr::select(cause_id,label,sex,metric,series,period,mean,lwr,upr)

  horizon_year <- tryCatch(projection_common_max_year_from_res_both(res_both, policy = "endogenous_max"), error = function(e) NA_integer_)
  if (is.finite(horizon_year)) {
    proj_tbl <- clip_to_year(proj_tbl, max_year = horizon_year, year_var = "period")
  }

  horizon_zone <- dplyr::bind_rows(
    tryCatch(res_both$resM$diag$projection_horizon_year, error = function(e) NULL),
    tryCatch(res_both$resF$diag$projection_horizon_year, error = function(e) NULL)
  )
  if (is.data.frame(horizon_zone) && nrow(horizon_zone)) {
    horizon_zone <- horizon_zone %>%
      dplyr::transmute(sex = as.character(sex), period = suppressWarnings(as.integer(period)), projection_zone = as.character(projection_zone))
    proj_tbl <- proj_tbl %>% dplyr::left_join(horizon_zone, by = c("sex", "period"))
  } else {
    proj_tbl$projection_zone <- NA_character_
  }

  hist_cutoff <- tryCatch(suppressWarnings(as.integer(res_both$combined$last_hist_year)[1]), error = function(e) NA_integer_)
  if (is.finite(hist_cutoff)) {
    proj_tbl <- proj_tbl %>%
      dplyr::mutate(
        projection_zone = dplyr::case_when(
          period <= hist_cutoff ~ "historical",
          TRUE ~ dplyr::coalesce(projection_zone, "beyond_max")
        ),
        projection_zone = factor(projection_zone, levels = c("historical", "credible", "caution", "risky", "beyond_max"))
      )
  }

  proj_tbl
}


# Suma TOTAL de causas: método normal-approx para IC
agregar_todas_causas <- function(proj_tbl, method = c("normal","conservador")) {
  method <- match.arg(method)

  expected_causes <- length(unique(stats::na.omit(proj_tbl$cause_id)))
  if (is.finite(expected_causes) && expected_causes > 0 && "cause_id" %in% names(proj_tbl)) {
    proj_tbl <- proj_tbl %>%
      dplyr::group_by(sex, metric, series, period) %>%
      dplyr::mutate(.n_causes_present = dplyr::n_distinct(cause_id)) %>%
      dplyr::ungroup() %>%
      dplyr::filter(.n_causes_present >= expected_causes) %>%
      dplyr::select(-.n_causes_present)
  }
  
  if (method == "normal") {
    by_sex <- proj_tbl %>%
      dplyr::group_by(sex, metric, series, period) %>%
      dplyr::summarise(
        mean = sum(mean, na.rm=TRUE),
        se2  = sum( ((upr - lwr)/(2*1.96))^2, na.rm=TRUE ),
        .groups="drop"
      ) %>%
      dplyr::mutate(
        lwr = pmax(0, mean - 1.96*sqrt(se2)),
        upr = mean + 1.96*sqrt(se2)
      ) %>% dplyr::select(-se2)
  } else {
    by_sex <- proj_tbl %>%
      dplyr::group_by(sex, metric, series, period) %>%
      dplyr::summarise(
        mean = sum(mean, na.rm=TRUE),
        lwr  = sum(lwr,  na.rm=TRUE),
        upr  = sum(upr,  na.rm=TRUE),
        .groups="drop"
      )
  }
  
  total <- by_sex %>%
    dplyr::group_by(metric, series, period) %>%
    dplyr::summarise(
      mean = sum(mean, na.rm=TRUE),
      lwr  = sum(lwr,  na.rm=TRUE),
      upr  = sum(upr,  na.rm=TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(sex = "T") %>%
    dplyr::select(sex, dplyr::everything())
  
  dplyr::bind_rows(by_sex, total) %>%
    dplyr::arrange(metric, series, sex, period)
}


# =============================================================


map_to_cause_id <- function(x) {
  k <- norm_txt(x)
  if (grepl("\bpulmon\b|\blung\b", k))                    return("lung")
  if (grepl("\bpancreas\b", k))                             return("pancreas")
  if (grepl("\brinon\b|\bkidney\b", k))                return("kidney")
  if (grepl("\bvejiga\b|\bbladder\b", k))              return("bladder")
  if (grepl("\blaringe\b|\blarynx\b", k))              return("larynx")
  if (grepl("\bestomago\b|\bstomach\b", k))            return("stomach")
  if (grepl("\besofago\b|\besophagus\b", k))           return("esophagus")
  if (grepl("cuello.*utero|cervix", k))                           return("cervix")
  if (grepl("cavidad.*oral|faringe|oral.*phar|mouth|phary", k))   return("oralphar")
  return(k)
}

list_available_cause_ids <- function(out_dir = BASE_RESULTS_DIR) {
  dirs <- list.dirs(out_dir, recursive = FALSE, full.names = FALSE)
  ids  <- sub("^cause_", "", dirs[grepl("^cause_", dirs)])
  tibble::tibble(
    cause_id = ids,
    rds_path = file.path(out_dir, paste0("cause_", ids), paste0(ids, "_res_both.rds")),
    exists   = file.exists(file.path(out_dir, paste0("cause_", ids), paste0(ids, "_res_both.rds")))
  )
}

get_res_both_from_rds <- function(cause,
                                  scenario = "freeze",
                                  out_dir  = BASE_RESULTS_DIR,
                                  strict   = TRUE) {
  id <- map_to_cause_id(cause)
  dir_cause <- file.path(out_dir, paste0("cause_", id))
  if (!dir.exists(dir_cause)) {
    if (isTRUE(strict)) stop("No existe la carpeta: ", dir_cause, call. = FALSE) else return(NULL)
  }
  avail <- list.files(dir_cause, pattern = sprintf("^%s_res_both_.*\\.rds$", id), full.names = TRUE)
  target <- list.files(dir_cause, pattern = sprintf("^%s_res_both_%s\\.rds$", id, scenario), full.names = TRUE)
  if (length(target) == 0L && identical(scenario, "freeze")) {
    fallback <- file.path(dir_cause, sprintf("%s_res_both.rds", id))
    if (file.exists(fallback)) target <- fallback
  }
  if (length(target) == 0L) {
    if (isTRUE(strict)) {
      disp <- if (length(avail) == 0L) "(ninguno)" else paste(basename(avail), collapse = ", ")
      stop(sprintf("No encuentro RDS para '%s' (id='%s') con scenario='%s'. Disponibles: %s", cause, id, scenario, disp), call. = FALSE)
    } else {
      return(NULL)
    }
  }
  readRDS(target[[1]])
}

build_total_incidence_from_rds <- function(cause_ids = causes$cause_id,
                                           scenario = "freeze",
                                           out_dir = BASE_RESULTS_DIR) {
  rows <- list()
  add_if_ok <- function(tbl, sx) {
    if (!is.null(tbl) && all(c("period","cases_hat","cases_lwr","cases_upr") %in% names(tbl))) {
      rows[[length(rows) + 1]] <<- dplyr::transmute(
        tbl,
        sex    = sx,
        period = as.integer(period),
        mean   = as.numeric(cases_hat),
        lwr    = as.numeric(cases_lwr),
        upr    = as.numeric(cases_upr)
      )
    }
  }
  for (id in cause_ids) {
    rb <- try(get_res_both_from_rds(id, scenario = scenario, out_dir = out_dir, strict = FALSE), silent = TRUE)
    if (inherits(rb, "try-error") || is.null(rb)) next
    if (!is.null(rb$resM)) add_if_ok(rb$resM$inc_annual_cond, "M")
    if (!is.null(rb$resF)) add_if_ok(rb$resF$inc_annual_cond, "F")
  }
  if (length(rows) == 0) {
    return(tibble::tibble(sex = character(), period = integer(), mean = double(), lwr = double(), upr = double()))
  }
  by_sex <- dplyr::bind_rows(rows)
  tot <- by_sex %>%
    dplyr::group_by(sex, period) %>%
    dplyr::summarise(mean = sum(mean, na.rm = TRUE), lwr = sum(lwr, na.rm = TRUE), upr = sum(upr, na.rm = TRUE), .groups = "drop")
  tot_T <- tot %>%
    dplyr::group_by(period) %>%
    dplyr::summarise(mean = sum(mean), lwr = sum(lwr), upr = sum(upr), .groups = "drop") %>%
    dplyr::mutate(sex = "T")
  dplyr::bind_rows(tot, tot_T) %>% dplyr::arrange(sex, period)
}

build_total_mortality_from_rds <- function(cause_ids = causes$cause_id,
                                           scenario = "freeze",
                                           out_dir = BASE_RESULTS_DIR) {
  rows <- list()
  for (id in cause_ids) {
    rb <- try(get_res_both_from_rds(id, scenario = scenario, out_dir = out_dir, strict = FALSE), silent = TRUE)
    if (inherits(rb, "try-error") || is.null(rb)) next
    anch <- rb$combined$annual_anchor
    if (!is.null(anch) && all(c("period","deaths_hat","deaths_lwr","deaths_upr") %in% names(anch))) {
      rows[[length(rows) + 1]] <- dplyr::transmute(
        anch,
        period = as.integer(period), mean = as.numeric(deaths_hat), lwr = as.numeric(deaths_lwr), upr = as.numeric(deaths_upr)
      )
    }
  }
  if (length(rows) == 0) {
    return(tibble::tibble(period = integer(), mean = double(), lwr = double(), upr = double(), sex = character()))
  }
  dplyr::bind_rows(rows) %>%
    dplyr::group_by(period) %>%
    dplyr::summarise(mean = sum(mean, na.rm = TRUE), lwr = sum(lwr, na.rm = TRUE), upr = sum(upr, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(sex = "T")
}
