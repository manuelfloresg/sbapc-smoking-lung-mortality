# runs/replication_diagnostics.R
# =============================================================
# REPLICATION HUB: Section 4 (Simulations) & Appendix C
# =============================================================

source("runs/_runtime_setup.R")
source("runs/_source_all.R")
source("adapters/build_inputs_sim.R")
source("R/09_figures_maintext.R")

library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
library(patchwork)
library(future)
library(future.apply)

# --- Canonical Configuration for Paper ---
CANONICAL_SEEDS <- 1:50
CANONICAL_DGPS  <- c("spec_linear")
CANONICAL_SCENS <- c("freeze", "up1pc", "down3pc", "quit")
CAUSE_ID        <- "lung"

# Output Directories
OUT_BASE    <- "results/20260515_FINAL_PROD"
OUT_SEC4    <- file.path(OUT_BASE, "section4")
OUT_APPENDIX <- file.path(OUT_BASE, "appendixC")
OUT_RAW     <- file.path(OUT_BASE, "raw_data")

dir.create(OUT_SEC4, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_APPENDIX, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_RAW, recursive = TRUE, showWarnings = FALSE)

# =============================================================
# 1. SIMULATION RUNNER
# =============================================================

run_single_seed_replication <- function(seed, dgp, scens = CANONICAL_SCENS, force_rerun = FALSE, ...) {
  # Check if all scenarios for this seed/dgp exist
  all_exist <- all(vapply(scens, function(sc) {
    file.exists(file.path(OUT_RAW, sprintf("res_%s_s%d_%s.rds", dgp, seed, sc)))
  }, logical(1)))
  
  if (all_exist && !force_rerun) {
    return(NULL)
  }
  
  message("\n>>> Processing SEED: ", seed, " | DGP: ", dgp)
  sim_base <- simulate_PIM_data(cause_id = CAUSE_ID, seed = seed, dgp = dgp, scenario_name = "freeze", ...)
  inputs <- build_inputs_sim(sim_base, cause_id = CAUSE_ID)
  
  cfg_row <- tibble::tibble(
    cause_id = CAUSE_ID,
    AGE_M_MIN = sim_base$meta$age_min, AGE_M_MAX = sim_base$meta$age_max,
    AGE_P_MIN = sim_base$meta$age_min, AGE_P_MAX = sim_base$meta$age_max,
    AGE_I_MIN = sim_base$meta$age_min, AGE_I_MAX = sim_base$meta$age_max,
    L_I_MAX_YEARS = 3L,
    MORT_SHOCK_YEARS = list(integer(0)),
    DOWNWEIGHT_F = list(integer(0))
  )
  
  # 1) Always run FREEZE first as the benchmark
  message("  Simulating base: freeze")
  sim_freeze <- simulate_PIM_data(cause_id = CAUSE_ID, seed = seed, dgp = dgp, scenario_name = "freeze", ...)
  prev_cfg_freeze <- get_prev_config(scenario = "freeze")
  res_freeze <- run_pipeline_both_from_inputs(inputs = inputs, cfg_row = cfg_row, prev_cfg = prev_cfg_freeze, ...)
  # SANITIZAR INMEDIATAMENTE para evitar errores de memoria (buffer overflow en serialize)
  if (!is.null(res_freeze$resM)) res_freeze$resM <- sanitize_pipeline_output(res_freeze$resM)
  if (!is.null(res_freeze$resF)) res_freeze$resF <- sanitize_pipeline_output(res_freeze$resF)
  res_freeze$meta  <- list(seed = seed, dgp = dgp, scenario = "freeze", args = list(...))
  # No adjuntamos la verdad todavía para mantener res_freeze liviano para el rebuilder
  
  # 2) Rebuild other scenarios from the freeze benchmark
  other_scens <- setdiff(scens, "freeze")
  gc() # Limpiar memoria antes del bucle pesado
  for (scen in other_scens) {
    message("  Rebuilding: ", scen)
    prev_cfg_scen <- get_prev_config(scenario = scen)
    
    out_rebuild <- .rebuild_scenario_freeze_benchmark(
      res_base = res_freeze, 
      inputs = inputs,
      cfg_row = cfg_row,
      prev_cfg_scen = prev_cfg_scen
    )
    
    # Capture truth for this specific scenario
    sim_scen <- simulate_PIM_data(cause_id = CAUSE_ID, seed = seed, dgp = dgp, scenario_name = scen, ...)
    
    res_scen <- out_rebuild$res_scen
    res_scen$meta  <- list(seed = seed, dgp = dgp, scenario = scen, args = list(...))
    res_scen$truth <- sim_scen$truth
    res_scen$inc_truth_grid <- sim_scen$inc_truth_grid
    res_scen$mort_truth_grid <- sim_scen$mort_truth_grid
    # res_scen$pop_all <- sim_scen$pop_all # Demasiado pesado
    
    # Now sanitize and save the scenario result
    if (exists("sanitize_pipeline_output", inherits = TRUE)) {
      if (!is.null(res_scen$resM)) res_scen$resM <- sanitize_pipeline_output(res_scen$resM)
      if (!is.null(res_scen$resF)) res_scen$resF <- sanitize_pipeline_output(res_scen$resF)
    }
    saveRDS(res_scen, file.path(OUT_RAW, sprintf("res_%s_s%d_%s.rds", dgp, seed, scen)))
  }
  
  # 3) Finally attach truth to freeze and save
  res_freeze$truth <- sim_freeze$truth
  res_freeze$inc_truth_grid <- sim_freeze$inc_truth_grid
  res_freeze$mort_truth_grid <- sim_freeze$mort_truth_grid
  # res_freeze$pop_all <- sim_freeze$pop_all 
  
  if (exists("sanitize_pipeline_output", inherits = TRUE)) {
    if (!is.null(res_freeze$resM)) res_freeze$resM <- sanitize_pipeline_output(res_freeze$resM)
    if (!is.null(res_freeze$resF)) res_freeze$resF <- sanitize_pipeline_output(res_freeze$resF)
  }
  saveRDS(res_freeze, file.path(OUT_RAW, sprintf("res_%s_s%d_freeze.rds", dgp, seed)))
}

run_simulation_replication <- function(seeds = CANONICAL_SEEDS, 
                                       dgps = CANONICAL_DGPS, 
                                       scens = CANONICAL_SCENS,
                                       force_rerun = FALSE,
                                       n_cores = NULL, ...) {
  
  if (!is.null(n_cores) && n_cores > 1) {
    plan(multisession, workers = n_cores)
    message("Using parallel execution with ", n_cores, " cores.")
  } else {
    plan(sequential)
  }
  
  tasks <- expand.grid(seed = seeds, dgp = dgps, stringsAsFactors = FALSE)
  
  future_lapply(seq_len(nrow(tasks)), function(i) {
    seed <- tasks$seed[i]
    dgp <- tasks$dgp[i]
    run_single_seed_replication(seed = seed, dgp = dgp, scens = scens, force_rerun = force_rerun, ...)
  }, future.seed = TRUE)
  
  plan(sequential) # Reset
}

# =============================================================
# 2. DATA EXTRACTION
# =============================================================

# --- Helpers ---
read_rds_safe <- function(file) {
  rb <- try(readRDS(file), silent = TRUE)
  if (inherits(rb, "try-error")) return(rb)
  
  # Fix-up if the RDS was saved with the buggy JSON roundtrip (lists instead of data frames)
  for (sx in c("resM", "resF")) {
    if (!is.null(rb[[sx]]$diag)) {
      for (df_name in c("z_prev_hist", "z_prev_future")) {
        if (is.list(rb[[sx]]$diag[[df_name]]) && !is.data.frame(rb[[sx]]$diag[[df_name]])) {
          rb[[sx]]$diag[[df_name]] <- as.data.frame(rb[[sx]]$diag[[df_name]])
        }
      }
    }
  }
  return(rb)
}

extract_all_metrics <- function(seeds = CANONICAL_SEEDS, dgps = CANONICAL_DGPS, scens = CANONICAL_SCENS) {
  metrics_list <- list()
  deltas_list  <- list()
  support_list <- list()
  inc_list     <- list()
  mort_list    <- list()
  
  for (seed in seeds) {
    for (dgp in dgps) {
      # First, get freeze mort for deltas
      freeze_rds <- file.path(OUT_RAW, sprintf("res_%s_s%d_freeze.rds", dgp, seed))
      if (!file.exists(freeze_rds)) next
      rb_freeze <- readRDS(freeze_rds)
      if (is.null(rb_freeze$inc_truth_grid)) {
        message("  WARNING: Skipping freeze file without embedded truth: ", freeze_rds)
        next
      }
      
      # We need truth for freeze to get diag_res$mort
      diag_freeze <- compare_pipeline_to_truth(rb_freeze, rb_freeze, out_dir = NULL)
      
      freeze_mort <- diag_freeze$mort %>% dplyr::select(period, sex, deaths_freeze = deaths_hat)
      
      for (scen in scens) {
        rds_file <- file.path(OUT_RAW, sprintf("res_%s_s%d_%s.rds", dgp, seed, scen))
        if (!file.exists(rds_file)) next
        # Read and process with robustness
        rb <- read_rds_safe(rds_file)
        if (inherits(rb, "try-error") || is.null(rb$inc_truth_grid)) {
          message("  WARNING: Skipping corrupt or incomplete file: ", rds_file)
          next
        }
        
        diag_res <- compare_pipeline_to_truth(rb, rb, out_dir = NULL)
        
        # Metrics
        m <- diag_res$metrics %>% mutate(seed = seed, dgp = dgp, scenario = scen)
        metrics_list[[paste(seed, dgp, scen, sep="_")]] <- m
        
        # Detailed data for plotting
        inc_list[[paste(seed, dgp, scen, sep="_")]] <- diag_res$inc %>% mutate(seed = seed, dgp = dgp, scenario = scen)
        mort_list[[paste(seed, dgp, scen, sep="_")]] <- diag_res$mort %>% mutate(seed = seed, dgp = dgp, scenario = scen)
        
        # Support
        if (!is.null(diag_res$support)) {
          support_list[[paste(seed, dgp, scen, sep="_")]] <- diag_res$support %>% mutate(seed = seed, dgp = dgp, scenario = scen)
        }
        
        # Deltas
        if (scen != "freeze") {
          current_mort <- diag_res$mort %>% dplyr::select(period, sex, deaths_hat)
          delta_df <- current_mort %>%
            dplyr::left_join(freeze_mort, by = c("period", "sex")) %>%
            dplyr::mutate(
              delta_deaths = deaths_hat - deaths_freeze,
              rel_delta = delta_deaths / pmax(deaths_freeze, 1e-12),
              seed = seed, dgp = dgp, scenario = scen
            )
          deltas_list[[paste(seed, dgp, scen, sep="_")]] <- delta_df
        }
      }
    }
  }
  
  return(list(
    metrics = bind_rows(metrics_list),
    deltas  = bind_rows(deltas_list),
    support = bind_rows(support_list),
    inc     = bind_rows(inc_list),
    mort    = bind_rows(mort_list)
  ))
}

# =============================================================
# 3. PLOTTING FUNCTIONS
# =============================================================

plot_deconstruction_figure <- function(seed = 4, dgp = "spec_linear", scen = "quit") {
  rds_file <- file.path(OUT_RAW, sprintf("res_%s_s%d_%s.rds", dgp, seed, scen))
  if (!file.exists(rds_file)) stop("RDS not found for deconstruction.")
  rb <- read_rds_safe(rds_file)
  if (inherits(rb, "try-error")) stop("Corrupt RDS for deconstruction: ", rds_file)
  
  diag_res <- compare_pipeline_to_truth(rb, rb, out_dir = NULL)
  
  df_mort <- diag_res$mort
  
  # Preparar para plot comparativo con nomenclatura del usuario
  plot_df <- df_mort %>%
    dplyr::select(
      period, sex, 
      Truth = deaths_true, 
      `Informed SBAPC (M | I | P)` = deaths_hat, 
      `Uninformed SBAPC (M | I)` = deaths_noP, 
      `Pure BAPC (M)` = deaths_bapc
    ) %>%
    tidyr::pivot_longer(
      cols = c(Truth, `Informed SBAPC (M | I | P)`, `Uninformed SBAPC (M | I)`, `Pure BAPC (M)`), 
      names_to = "Series", values_to = "Deaths"
    )
  
  plot_df$Series <- factor(plot_df$Series, levels = c("Truth", "Informed SBAPC (M | I | P)", "Uninformed SBAPC (M | I)", "Pure BAPC (M)"))
  
  last_hist <- rb$meta$last_hist %||% 2022
  
  g <- ggplot(plot_df, aes(x = period, y = Deaths, color = Series, linetype = Series)) +
    geom_line(linewidth = 1) +
    geom_vline(xintercept = last_hist, linetype = "dotted", color = "gray50") +
    facet_wrap(~sex, scales = "free_y") +
    scale_color_manual(values = c(
      "Truth" = "black", 
      "Informed SBAPC (M | I | P)" = "#CD5C5C", 
      "Uninformed SBAPC (M | I)" = "#ff7f0e", 
      "Pure BAPC (M)" = "#4682B4"
    )) +
    scale_linetype_manual(values = c(
      "Truth" = "dashed", 
      "Informed SBAPC (M | I | P)" = "solid", 
      "Uninformed SBAPC (M | I)" = "dotdash", 
      "Pure BAPC (M)" = "dotted"
    )) +
    labs(title = "Information Gain Deconstruction",
         subtitle = sprintf("Seed %d | DGP: %s | Scenario: %s", seed, dgp, scen),
         y = "Annual Deaths", x = "Year") +
    theme_paper_main(base_size = 11) +
    theme(legend.position = "bottom")
  
  return(g)
}

plot_scenario_sensitivity_informed <- function(seed = 4, dgp = "spec_linear") {
  scens <- c("freeze", "up1pc", "down1pc", "quit")
  data_list <- list()
  
  for (sc in scens) {
    rds_file <- file.path(OUT_RAW, sprintf("res_%s_s%d_%s.rds", dgp, seed, sc))
    if (!file.exists(rds_file)) next
    rb <- read_rds_safe(rds_file)
    if (inherits(rb, "try-error")) next
    
    # Extract Informed SBAPC (annual_anchor)
    # resM and resF
    df_m <- rb$resM$annual_anchor %>% mutate(sex = "M", scenario = sc)
    df_f <- rb$resF$annual_anchor %>% mutate(sex = "F", scenario = sc)
    data_list[[sc]] <- bind_rows(df_m, df_f)
  }
  
  if (length(data_list) == 0) stop("No RDS found for sensitivity plot.")
  plot_df <- bind_rows(data_list)
  
  # Factor scenarios for legend
  plot_df$scenario <- factor(plot_df$scenario, levels = scens)
  
  g <- ggplot(plot_df, aes(x = period, y = deaths_hat, color = scenario, group = scenario)) +
    geom_line(linewidth = 1) +
    geom_vline(xintercept = 2022, linetype = "dotted", color = "gray50") +
    facet_wrap(~sex, scales = "free_y") +
    scale_color_manual(values = c(
      "freeze" = "gray30", 
      "up1pc" = "#CD5C5C", 
      "down1pc" = "#4682B4", 
      "quit" = "#228B22"
    )) +
    labs(title = sprintf("Scenario Sensitivity: Informed SBAPC (Seed %d)", seed),
         subtitle = "Comparing Informed SBAPC Projections across Policy Scenarios",
         y = "Projected Deaths", x = "Year", color = "Scenario") +
    theme_paper_main(base_size = 11)
  
  return(g)
}

plot_transmission_waterfall <- function(seed = 4, dgp = "spec_linear", scen = "quit") {
  # Panel A: Delta Smoking Stock (1 - p_never)
  # Panel B: Delta Incidence
  # Panel C: Delta Mortality
  
  rds_file <- file.path(OUT_RAW, sprintf("res_%s_s%d_%s.rds", dgp, seed, scen))
  rds_freeze <- file.path(OUT_RAW, sprintf("res_%s_s%d_freeze.rds", dgp, seed))
  if (!file.exists(rds_file) || !file.exists(rds_freeze)) stop("RDS not found for waterfall.")
  
  rb_scen <- read_rds_safe(rds_file)
  rb_frz  <- read_rds_safe(rds_freeze)
  if (inherits(rb_scen, "try-error") || inherits(rb_frz, "try-error")) {
    stop("Corrupt RDS for waterfall: ", rds_file, " or ", rds_freeze)
  }
  
  sim_args <- rb_scen$meta$args %||% list()
  sim_scen <- rb_scen
  sim_frz  <- rb_frz
  
  get_stock <- function(rb) {
    # Combinar z_prev_hist y z_prev_future para ambos sexos
    data_list <- list()
    for (sx in c("M", "F")) {
      sex_res <- if (sx == "M") rb$resM else rb$resF
      if (is.null(sex_res)) next
      
      # Rates all (for incidence rate)
      r_all <- sex_res$inc_fit$rates_all %>% mutate(sex = sx)
      
      # Prevalence stock (from diag)
      # p_cur is the current smoker proportion
      z_h <- sex_res$diag$z_prev_hist
      z_f <- sex_res$diag$z_prev_future
      z_all <- bind_rows(z_h, z_f) %>% 
        group_by(period) %>% 
        summarise(current_prev = mean(as.numeric(p_cur), na.rm = TRUE), .groups = "drop") %>%
        mutate(sex = sx)
      
      # Join
      res_sx <- r_all %>%
        group_by(period, sex) %>%
        summarise(inc_rate = mean(rate_hat, na.rm = TRUE), .groups = "drop") %>%
        left_join(z_all, by = c("period", "sex"))
      
      data_list[[sx]] <- res_sx
    }
    bind_rows(data_list)
  }
  
  stock_scen <- get_stock(rb_scen)
  stock_frz  <- get_stock(rb_frz)
  
  diff_stock <- stock_scen %>%
    left_join(stock_frz, by = c("period", "sex"), suffix = c("_scen", "_frz")) %>%
    mutate(
      diff_pcur = current_prev_scen - current_prev_frz,
      diff_inc  = (inc_rate_scen - inc_rate_frz) / pmax(inc_rate_frz, 1e-12) * 100
    )
  
  # Panel A: Current Prevalence change
  pA <- ggplot(diff_stock, aes(x = period, y = diff_pcur, color = sex)) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    geom_line(linewidth = 1) +
    scale_y_continuous(limits = c(-0.4, 0.4)) +
    labs(title = "Stage 1: Change in Current Prevalence", y = "Delta Prevalence", x = NULL) +
    theme_paper_main(base_size = 10)
  
  # Panel B: Incidence Rate % change
  pB <- ggplot(diff_stock, aes(x = period, y = diff_inc, color = sex)) +
    geom_line(linewidth = 1) +
    labs(title = "Stage 2: Change in Incidence Rate", y = "% Delta Rate", x = NULL) +
    theme_paper_main(base_size = 10)
  
  # Panel C: Mortality % change (Stratified by sex)
  get_mort_diff <- function(sex_res_scen, sex_res_frz, sex_lab) {
    m_scen <- sex_res_scen$annual_anchor %>% select(period, deaths_scen = deaths_hat)
    m_frz  <- sex_res_frz$annual_anchor %>% select(period, deaths_frz = deaths_hat)
    m_scen %>%
      left_join(m_frz, by = "period") %>%
      mutate(diff_pct = (deaths_scen - deaths_frz) / pmax(deaths_frz, 1e-12) * 100, sex = sex_lab)
  }
  
  diff_mort <- bind_rows(
    get_mort_diff(rb_scen$resM, rb_frz$resM, "M"),
    get_mort_diff(rb_scen$resF, rb_frz$resF, "F")
  )
  
  pC <- ggplot(diff_mort, aes(x = period, y = diff_pct, color = sex)) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    geom_line(linewidth = 1) +
    labs(title = "Stage 3: Change in Total Deaths", y = "% Delta Deaths", x = "Year") +
    theme_paper_main(base_size = 10)
  
  (pA / pB / pC) + plot_annotation(title = sprintf("Transmission Waterfall: %s Scenario (Seed %d)", scen, seed)) &
    theme_paper_main(base_size = 11)
}

plot_reliability_calibration <- function(data) {
  # data$inc contains the errors per period/sex/seed/dgp
  # but we need to join with support info to get 'horizon' and 'support_frac'
  
  # For the calibration plot, we aggregate all seeds and DGPs
  # to show how error grows with horizon.
  
  df <- data$inc %>%
    dplyr::filter(period > 2022) %>%
    dplyr::mutate(horizon = period - 2022)
  
  # Join with support info if available, otherwise just use horizon
  g <- ggplot(df, aes(x = horizon, y = abs(rel_error) * 100)) +
    stat_summary(fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.2, fill = "blue") +
    stat_summary(fun = "mean", geom = "line", linewidth = 1, color = "blue") +
    geom_vline(xintercept = c(5, 10, 20), linetype = "dashed", color = "gray60") +
    annotate("text", x = 2.5, y = Inf, label = "Credible", vjust = 1.5, size = 3.5, family = "serif") +
    annotate("text", x = 7.5, y = Inf, label = "Caution", vjust = 1.5, size = 3.5, family = "serif") +
    labs(title = "Reliability Calibration: Error vs. Horizon",
         y = "Mean Absolute Relative Error (%)", x = "Projection Horizon (Years)") +
    theme_paper_main(base_size = 11)
  
  return(g)
}


# =============================================================
# 4. MAIN OUTPUT GENERATOR
# =============================================================

replicate_main_paper <- function() {
  # 1. Deconstruction
  g1 <- plot_deconstruction_figure(seed = 4, dgp = "spec_linear", scen = "quit")
  ggsave(file.path(OUT_SEC4, "fig_deconstruction_seed4.png"), g1, width = 10, height = 6, bg = "white")
  
  # 2. Waterfall
  g2 <- plot_transmission_waterfall(seed = 4, dgp = "spec_linear", scen = "quit")
  ggsave(file.path(OUT_SEC4, "fig_waterfall_seed4.png"), g2, width = 8, height = 10, bg = "white")
  
  # 3. Scenario Sensitivity (New)
  g_sens <- plot_scenario_sensitivity_informed(seed = 4, dgp = "spec_linear")
  ggsave(file.path(OUT_SEC4, "fig_sensitivity_seed4.png"), g_sens, width = 10, height = 6, bg = "white")
  
  # 4. Bias Table
  data <- extract_all_metrics()
  export_latex_bias_summary(data$metrics, file.path(OUT_SEC4, "tab_bias_summary.tex"))
  write_csv(data$metrics, file.path(OUT_RAW, "all_metrics.csv"))
  
  # 5. Reliability Plot
  g3 <- plot_reliability_calibration(data)
  ggsave(file.path(OUT_SEC4, "fig_reliability_calibration.png"), g3, width = 10, height = 6, bg = "white")
  
  # 6. Support Summary
  write_csv(data$support, file.path(OUT_SEC4, "support_summary.csv"))
  
  message("\nMain paper replication files generated in: ", OUT_SEC4)
}

# =============================================================
# 5. APPENDIX C GENERATOR
# =============================================================

generate_appendix_c <- function() {
  data <- extract_all_metrics()
  
  # 1. Distribution of Bias (Boxplots)
  g_bias <- ggplot(data$metrics, aes(x = scenario, y = proj_bias, fill = scenario)) +
    geom_boxplot(alpha = 0.7) +
    facet_wrap(~dgp + sex) +
    labs(title = "Appendix C: Distribution of Projection Bias across Seeds",
         y = "Projection Bias (%)", x = "Scenario") +
    theme_minimal() + theme(legend.position = "none")
  ggsave(file.path(OUT_APPENDIX, "fig_bias_distributions.png"), g_bias, width = 10, height = 7, bg = "white")
  
  # 2. Case Studies Selection (Based on Projection Bias in 'quit' scenario)
  # We'll pick Best (min abs bias), Median, and Worst (max abs bias)
  case_seeds <- data$metrics %>%
    dplyr::filter(scenario == "quit", dgp == "spec_linear", sex == "M") %>%
    dplyr::arrange(abs(proj_bias)) %>%
    dplyr::slice(c(1, n() %/% 2, n())) %>%
    dplyr::mutate(label = c("Best", "Median", "Worst"))
  
  for (i in 1:nrow(case_seeds)) {
    s <- case_seeds$seed[i]
    lbl <- case_seeds$label[i]
    g_case <- plot_deconstruction_figure(seed = s, dgp = "spec_linear", scen = "quit") +
      labs(title = sprintf("Appendix C Case Study: %s Performance (Seed %d)", lbl, s))
    ggsave(file.path(OUT_APPENDIX, sprintf("fig_case_study_%s_s%d.png", tolower(lbl), s)), g_case, width = 10, height = 6, bg = "white")
  }
  
  # 3. Full Detailed Table (CSV)
  write_csv(data$metrics, file.path(OUT_APPENDIX, "full_simulation_matrix.csv"))
  
  message("\nAppendix C replication files generated in: ", OUT_APPENDIX)
}

# =============================================================
# 6. ORCHESTRATOR
# =============================================================

replicate_all_simulations <- function() {
  message("STARTING FULL REPLICATION WORKFLOW...")
  
  # 1. Run Simulations
  # Using 6 cores to avoid memory allocation errors (INLA is memory intensive)
  run_simulation_replication(n_cores = 6)
  
  # 2. Section 4
  replicate_main_paper()
  
  # 3. Appendix C
  generate_appendix_c()
  
  message("\nALL REPLICATION TASKS COMPLETED SUCCESSFULLY.")
}
