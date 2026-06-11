# =============================================================
# Trace and Transmission Diagnostics Export
# =============================================================

export_scenario_transmission_diagnostics <- function(run_list, out_dir, ref_scenario = "freeze", scenarios = NULL, stem = "transmission_trace") {
  if (is.null(scenarios)) scenarios <- names(run_list)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  diag_prev_all <- list()
  inc_rates_all <- list()
  mort_pred_all <- list()
  mort_kernel_all <- list()
  
  for (scn in scenarios) {
    if (!scn %in% names(run_list)) next
    res <- run_list[[scn]]$res
    
    for (sex in c("M", "F")) {
      sex_res <- if (sex == "M") res$resM else res$resF
      if (is.null(sex_res)) next
      
      # 1. Prevalence to Incidence transmission metrics
      if (!is.null(sex_res$diag$prev) && nrow(sex_res$diag$prev)) {
        dp <- sex_res$diag$prev %>% dplyr::mutate(scenario = scn)
        diag_prev_all[[paste(scn, sex)]] <- dp
      }
      
      # 2. Incidence detailed TRACE (rates, offsets)
      if (!is.null(sex_res$inc_fit$rates_all_full) && nrow(sex_res$inc_fit$rates_all_full)) {
        ir <- sex_res$inc_fit$rates_all_full %>% dplyr::mutate(scenario = scn)
        inc_rates_all[[paste(scn, sex)]] <- ir
      }
      
      # 3. Incidence to Mortality transmission metrics (Kernel)
      if (!is.null(sex_res$diag$mort_kernel) && nrow(sex_res$diag$mort_kernel)) {
        mk <- sex_res$diag$mort_kernel %>% dplyr::mutate(scenario = scn, sex = sex)
        mort_kernel_all[[paste(scn, sex)]] <- mk
      }
      
      # 4. Mortality detailed TRACE (anchored predictions, external offsets)
      if (!is.null(sex_res$mort_anchor_pred_detail) && nrow(sex_res$mort_anchor_pred_detail)) {
        mp <- sex_res$mort_anchor_pred_detail %>% dplyr::mutate(scenario = scn)
        mort_pred_all[[paste(scn, sex)]] <- mp
      }
    }
  }
  
  files_saved <- character()
  
  if (length(diag_prev_all)) {
    out_prev <- dplyr::bind_rows(diag_prev_all)
    f_prev <- file.path(out_dir, paste0(stem, "_prev_to_inc_summary.csv"))
    readr::write_csv(out_prev, f_prev)
    files_saved <- c(files_saved, f_prev)
  }
  
  if (length(inc_rates_all)) {
    out_inc <- dplyr::bind_rows(inc_rates_all)
    f_inc <- file.path(out_dir, paste0(stem, "_incidence_detailed_trace.csv"))
    readr::write_csv(out_inc, f_inc)
    files_saved <- c(files_saved, f_inc)
  }
  
  if (length(mort_kernel_all)) {
    out_mk <- dplyr::bind_rows(mort_kernel_all)
    f_mk <- file.path(out_dir, paste0(stem, "_mortality_kernel.csv"))
    readr::write_csv(out_mk, f_mk)
    files_saved <- c(files_saved, f_mk)
  }
  
  if (length(mort_pred_all)) {
    out_mort <- dplyr::bind_rows(mort_pred_all)
    f_mort <- file.path(out_dir, paste0(stem, "_mortality_detailed_trace.csv"))
    readr::write_csv(out_mort, f_mort)
    files_saved <- c(files_saved, f_mort)
  }
  
  invisible(files_saved)
}
