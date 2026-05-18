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
# ------------------------------------------------------------
# PROJECT PARAMETERS
# ------------------------------------------------------------
CANONICAL_SEEDS <- 1:50
CANONICAL_DGPS  <- c("spec_linear")
CANONICAL_SCENS <- c("freeze", "up1pc", "down1pc", "quit")
CAUSE_ID        <- "lung"

# Scenario Palette
SCEN_COLORS <- c(
  "up1pc"   = "#D73027", # Red
  "freeze"  = "#FDAE61", # Orange
  "down1pc" = "#00796B",
  "quit"    = "#512DA8"
)

SCEN_LABELS <- c(
  "up1pc"   = "\u2191 1% p.a.",
  "freeze"  = "Freeze (2022)",
  "down1pc" = "\u2193 1% p.a.",
  "quit"    = "Quit"
)

# Output Directories
.out_base_env <- Sys.getenv("BAPC_OUT_BASE", unset = "")
OUT_BASE    <- if (nzchar(.out_base_env)) .out_base_env else file.path(BAPC_PATHS$results, "20260515_FINAL_PROD")
OUT_BASE    <- normalizePath(OUT_BASE, winslash = "/", mustWork = FALSE)
OUT_SEC4    <- file.path(OUT_BASE, "section4")
OUT_APPENDIX <- file.path(OUT_BASE, "appendixC")
OUT_RAW     <- file.path(OUT_BASE, "raw_data")
OUT_RAW_ORACLE <- file.path(OUT_BASE, "raw_data_oracle")

dir.create(OUT_SEC4, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_APPENDIX, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_RAW, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_RAW_ORACLE, recursive = TRUE, showWarnings = FALSE)

FIG_FORMAT <- {
  val <- Sys.getenv("BAPC_FIG_FORMAT", unset = "")
  if (!nzchar(val)) val <- getOption("BAPC_FIG_FORMAT", BAPC_FIG_FORMAT %||% "svg")
  match.arg(as.character(val)[1], c("svg", "png", "both"))
}

figure_exts <- function(format = FIG_FORMAT) {
  if (identical(format, "both")) c("svg", "png") else format
}

save_paper_plot <- function(plot, path_no_ext, width, height, bg = "white", format = FIG_FORMAT, ...) {
  for (ext in figure_exts(format)) {
    ggplot2::ggsave(
      filename = paste0(path_no_ext, ".", ext),
      plot = plot,
      width = width,
      height = height,
      bg = bg,
      ...
    )
  }
  invisible(path_no_ext)
}

# =============================================================
# 1. SIMULATION RUNNER
# =============================================================

.simulation_information_args <- function(information_set = c("realistic", "oracle")) {
  information_set <- match.arg(as.character(information_set)[1], c("realistic", "oracle"))
  if (identical(information_set, "oracle")) {
    return(list(
      prev_obs_period_min = -9999L,
      prev_obs_period_max = 2022L,
      prev_obs_age_min = AGE_P_MIN,
      prev_obs_age_max = 9999L
    ))
  }
  list()
}

run_single_seed_replication <- function(seed, dgp, scens = CANONICAL_SCENS, force_rerun = FALSE,
                                        information_set = c("realistic", "oracle"),
                                        raw_dir = NULL, ...) {
  information_set <- match.arg(as.character(information_set)[1], c("realistic", "oracle"))
  if (is.null(raw_dir)) raw_dir <- if (identical(information_set, "oracle")) OUT_RAW_ORACLE else OUT_RAW
  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
  sim_info_args <- .simulation_information_args(information_set)
  sim_call <- function(scenario) {
    do.call(simulate_PIM_data, c(
      list(cause_id = CAUSE_ID, seed = seed, dgp = dgp, scenario_name = scenario),
      sim_info_args,
      list(...)
    ))
  }

  # Check if all scenarios for this seed/dgp exist
  all_exist <- all(vapply(scens, function(sc) {
    file.exists(file.path(raw_dir, sprintf("res_%s_s%d_%s.rds", dgp, seed, sc)))
  }, logical(1)))
  
  if (all_exist && !force_rerun) {
    return(NULL)
  }
  
  message("\n>>> Processing SEED: ", seed, " | DGP: ", dgp, " | information_set: ", information_set)
  sim_base <- sim_call("freeze")
  inputs <- build_inputs_sim(sim_base, cause_id = CAUSE_ID)
  
  # Use level trend and freeze for model alignment.
  extra_args <- list(trend_type = "level", gammaP_method = "freeze", sd_theta_IP = 2.0)
  
  cfg_row <- tibble::tibble(
    cause_id = CAUSE_ID,
    AGE_M_MIN = sim_base$meta$age_min, AGE_M_MAX = sim_base$meta$age_max,
    PERIOD_P_MIN = if (identical(information_set, "oracle")) sim_base$meta$prev_obs_period_min else PERIOD_M_MIN,
    PERIOD_P_MAX = if (identical(information_set, "oracle")) sim_base$meta$prev_obs_period_max else PERIOD_M_MAX,
    AGE_P_MIN = sim_base$meta$prev_obs_age_min %||% sim_base$meta$age_min_p %||% sim_base$meta$age_min,
    AGE_P_MAX = if (identical(information_set, "oracle")) {
      sim_base$meta$prev_obs_age_max %||% sim_base$meta$age_max
    } else {
      sim_base$meta$age_max_p %||% sim_base$meta$age_max
    },
    AGE_I_MIN = sim_base$meta$age_min, AGE_I_MAX = sim_base$meta$age_max,
    L_I_MAX_YEARS = 3L,
    MORT_SHOCK_YEARS = list(integer(0)),
    DOWNWEIGHT_F = list(integer(0))
  )
  
  # 1) Always run FREEZE first as the benchmark
  message("  Simulating base: freeze")
  sim_freeze <- sim_call("freeze")
  prev_cfg_freeze <- get_prev_config(scenario = "freeze")
  res_freeze <- do.call(run_pipeline_both_from_inputs, c(list(inputs = inputs, cfg_row = cfg_row, prev_cfg = prev_cfg_freeze), extra_args))
  # SANITIZE IMMEDIATELY to avoid serialization memory limits/buffer overflows
  if (!is.null(res_freeze$resM)) res_freeze$resM <- sanitize_pipeline_output(res_freeze$resM, keep_rebuilder = TRUE)
  if (!is.null(res_freeze$resF)) res_freeze$resF <- sanitize_pipeline_output(res_freeze$resF, keep_rebuilder = TRUE)
  res_freeze$meta  <- list(seed = seed, dgp = dgp, scenario = "freeze", information_set = information_set,
                           args = c(sim_info_args, list(...)))
  # Do not attach truth yet to keep res_freeze lightweight for the scenario rebuilder
  
  # 2) Rebuild other scenarios from the freeze benchmark
  other_scens <- setdiff(scens, "freeze")
  gc() # Free memory before launching parallel loops
  for (scen in other_scens) {
    message("  Rebuilding: ", scen)
    prev_cfg_scen <- get_prev_config(scenario = scen)
    
    out_rebuild <- do.call(.rebuild_scenario_freeze_benchmark, c(list(
      res_base = res_freeze, 
      inputs = inputs,
      cfg_row = cfg_row,
      prev_cfg_scen = prev_cfg_scen
    ), extra_args))
    
    # Capture truth for this specific scenario
    sim_scen <- sim_call(scen)
    
    res_scen <- out_rebuild$res_scen
    res_scen$meta  <- list(seed = seed, dgp = dgp, scenario = scen, information_set = information_set,
                           args = c(sim_info_args, list(...)))
    res_scen$truth <- sim_scen$truth
    res_scen$inc_truth_grid <- sim_scen$inc_truth_grid
    res_scen$mort_truth_grid <- sim_scen$mort_truth_grid
    # res_scen$pop_all <- sim_scen$pop_all # Too heavy
    
    # Now sanitize and save the scenario result
    if (exists("sanitize_pipeline_output", inherits = TRUE)) {
      if (!is.null(res_scen$resM)) res_scen$resM <- sanitize_pipeline_output(res_scen$resM, keep_rebuilder = FALSE)
      if (!is.null(res_scen$resF)) res_scen$resF <- sanitize_pipeline_output(res_scen$resF, keep_rebuilder = FALSE)
    }
    saveRDS(res_scen, file.path(raw_dir, sprintf("res_%s_s%d_%s.rds", dgp, seed, scen)))
  }
  
  # 3) Finally attach truth to freeze and save
  res_freeze$truth <- sim_freeze$truth
  res_freeze$inc_truth_grid <- sim_freeze$inc_truth_grid
  res_freeze$mort_truth_grid <- sim_freeze$mort_truth_grid
  # res_freeze$pop_all <- sim_freeze$pop_all 
  
  if (exists("sanitize_pipeline_output", inherits = TRUE)) {
    if (!is.null(res_freeze$resM)) res_freeze$resM <- sanitize_pipeline_output(res_freeze$resM, keep_rebuilder = FALSE)
    if (!is.null(res_freeze$resF)) res_freeze$resF <- sanitize_pipeline_output(res_freeze$resF, keep_rebuilder = FALSE)
  }
  saveRDS(res_freeze, file.path(raw_dir, sprintf("res_%s_s%d_freeze.rds", dgp, seed)))
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

compact_replication_rds <- function(raw_dir = OUT_RAW,
                                    compact_dir = file.path(dirname(OUT_RAW), "raw_data_compact"),
                                    pattern = "^res_.*\\.rds$",
                                    overwrite = FALSE) {
  dir.create(compact_dir, recursive = TRUE, showWarnings = FALSE)
  files <- list.files(raw_dir, pattern = pattern, full.names = TRUE)
  manifest <- vector("list", length(files))

  for (i in seq_along(files)) {
    src <- files[[i]]
    dst <- file.path(compact_dir, basename(src))
    if (file.exists(dst) && !isTRUE(overwrite)) {
      manifest[[i]] <- tibble::tibble(
        file = basename(src),
        status = "skipped_exists",
        source_bytes = file.info(src)$size,
        compact_bytes = file.info(dst)$size
      )
      next
    }
    rb <- read_rds_safe(src)
    if (inherits(rb, "try-error")) {
      manifest[[i]] <- tibble::tibble(
        file = basename(src),
        status = "read_error",
        source_bytes = file.info(src)$size,
        compact_bytes = NA_real_
      )
      next
    }
    if (!is.null(rb$resM)) rb$resM <- sanitize_pipeline_output(rb$resM, keep_rebuilder = FALSE)
    if (!is.null(rb$resF)) rb$resF <- sanitize_pipeline_output(rb$resF, keep_rebuilder = FALSE)
    saveRDS(rb, dst)
    manifest[[i]] <- tibble::tibble(
      file = basename(src),
      status = "compacted",
      source_bytes = file.info(src)$size,
      compact_bytes = file.info(dst)$size
    )
    rm(rb)
    if (i %% 10L == 0L) gc()
  }

  out <- dplyr::bind_rows(manifest) %>%
    dplyr::mutate(ratio = compact_bytes / pmax(source_bytes, 1))
  readr::write_csv(out, file.path(compact_dir, "compact_manifest.csv"))
  out
}

extract_all_metrics <- function(seeds = CANONICAL_SEEDS, dgps = CANONICAL_DGPS, scens = CANONICAL_SCENS, force_refresh = FALSE) {
  cache_file <- file.path(OUT_RAW, "all_extracted_data.rds")
  if (file.exists(cache_file) && !isTRUE(force_refresh)) {
    message("Loading cached metrics from: ", cache_file)
    cached <- readRDS(cache_file)
    cached_scens <- sort(unique(as.character(cached$metrics$scenario %||% character(0))))
    requested_scens <- sort(unique(as.character(scens)))
    if (identical(cached_scens, requested_scens)) {
      return(cached)
    }
    message("Cached metrics use a different scenario set; rebuilding extraction cache.")
  }
  
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
  
  res <- list(
    metrics = bind_rows(metrics_list),
    deltas  = bind_rows(deltas_list),
    support = bind_rows(support_list),
    inc     = bind_rows(inc_list),
    mort    = bind_rows(mort_list)
  )
  saveRDS(res, cache_file)
  return(res)
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
  
  # Prepare comparative plotting data using user-friendly series names
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

plot_scenario_atlas_by_sex <- function(seed = 4, sex_lab = "M") {
  # Scenarios to include in the requested order
  scens <- c("up1pc", "freeze", "down1pc", "quit")
  all_data <- list()
  
  for (scen in scens) {
    f_path <- file.path(OUT_RAW, sprintf("res_%s_s%d_%s.rds", "spec_linear", seed, scen))
    if (!file.exists(f_path)) next
    
    rb <- readRDS(f_path)
    sex_res <- if (sex_lab == "M") rb$resM else rb$resF
    
    # 1. Truth
    df_truth <- rb$mort_truth_grid[rb$mort_truth_grid$sex == sex_lab, ]
    df_truth <- aggregate(mort_deaths_scen_true ~ period, data = df_truth, FUN = sum, na.rm = TRUE)
    colnames(df_truth)[2] <- "deaths"
    df_truth$model <- "Truth"
    
    # 2. Informed SBAPC
    df_inf <- data.frame(period = sex_res$annual_anchor$period, deaths = sex_res$annual_anchor$deaths_hat, model = "Informed SBAPC")
    
    # 3. Uninformed SBAPC (No-Prev)
    df_uninf <- data.frame(period = sex_res$annual_anchor_noP$period, deaths = sex_res$annual_anchor_noP$deaths_hat, model = "Uninformed SBAPC")
    
    # 4. Pure BAPC
    df_bapc <- data.frame(period = sex_res$annual_bapc$period, deaths = sex_res$annual_bapc$deaths_hat, model = "Pure BAPC")
    
    combined <- bind_rows(df_truth, df_inf, df_uninf, df_bapc) %>%
      mutate(scenario = SCEN_LABELS[[scen]])
    
    all_data[[scen]] <- combined
  }
  
  df_final <- bind_rows(all_data) %>% filter(period > 1995)
  df_final$scenario <- factor(df_final$scenario, levels = SCEN_LABELS[scens])
  model_levels <- c("Truth", "Informed SBAPC", "Uninformed SBAPC", "Pure BAPC")
  df_final$model <- factor(df_final$model, levels = model_levels)
  
  # Colors and Types
  pal <- c("Truth" = "black", "Informed SBAPC" = "#D32F2F", "Uninformed SBAPC" = "#EF6C00", "Pure BAPC" = "#1976D2")
  types <- c("Truth" = "dashed", "Informed SBAPC" = "solid", "Uninformed SBAPC" = "longdash", "Pure BAPC" = "dotted")
  
  g <- ggplot(df_final, aes(x = period, y = deaths, color = model, linetype = model)) +
    facet_wrap(~scenario, nrow = 1, scales = "fixed") +
    geom_vline(xintercept = 2022, linetype = "dotted", color = "gray60") +
    geom_line(linewidth = 0.9) +
    scale_color_manual(values = pal, breaks = model_levels) +
    scale_linetype_manual(values = types, breaks = model_levels) +
    labs(title = sprintf("SBAPC Performance Atlas: %s (Seed %d)", ifelse(sex_lab=="M", "Males", "Females"), seed),
         subtitle = "Comparison across scenarios with shared Y-axis",
         y = "Annual Deaths", x = "Year", color = "Model", linetype = "Model") +
    theme_paper_main(base_size = 12) +
    theme(legend.position = "bottom", strip.background = element_rect(fill = "gray95"))
  
  return(g)
}

# =============================================================

plot_scenario_sensitivity_informed <- function(seed = 4, dgp = "spec_linear") {
  scens_to_plot <- c("up1pc", "freeze", "down1pc", "quit")
  
  data_list <- list()
  for (scen in scens_to_plot) {
    f_path <- file.path(OUT_RAW, sprintf("res_%s_s%d_%s.rds", dgp, seed, scen))
    if (!file.exists(f_path)) {
      message("Sensitivity informed: File not found ", f_path)
      next
    }
    res_both <- readRDS(f_path)
    
    # Extract aggregated total mortality (both sexes combined)
    m_m <- res_both$resM$annual_anchor %>% select(period, deaths_hat) %>% mutate(sex = "M")
    m_f <- res_both$resF$annual_anchor %>% select(period, deaths_hat) %>% mutate(sex = "F")
    
    data_list[[scen]] <- bind_rows(m_m, m_f) %>%
      group_by(period) %>%
      summarise(deaths = sum(deaths_hat, na.rm = TRUE), .groups = "drop") %>%
      mutate(scenario = scen)
  }
  
  df_all <- bind_rows(data_list)
  df_all$scenario <- factor(df_all$scenario, levels = scens_to_plot)
  last_hist <- 2022
  
  ggplot(df_all, aes(x = period, y = deaths, color = scenario)) +
    geom_vline(xintercept = last_hist, linetype = "dotted") +
    geom_line(linewidth = 1.2) +
    scale_color_manual(values = SCEN_COLORS, labels = SCEN_LABELS) +
    labs(title = "Policy Sensitivity: Projected Deaths by Scenario",
         subtitle = sprintf("Seed %d | DGP: %s | (Total Population)", seed, dgp),
         y = "Total Annual Deaths", x = "Year", color = "Scenario") +
    theme_paper_main()
}

plot_transmission_waterfall <- function(seed = 4, dgp = "spec_linear", scen = "quit") {
  # Panel A: Current Prevalence Level
  # Panel B: Effective Exposure (The 'Slide')
  # Panel C: Incidence Rate % change
  # Panel D: Mortality % change
  
  rds_file <- file.path(OUT_RAW, sprintf("res_%s_s%d_%s.rds", dgp, seed, scen))
  rds_freeze <- file.path(OUT_RAW, sprintf("res_%s_s%d_freeze.rds", dgp, seed))
  if (!file.exists(rds_file) || !file.exists(rds_freeze)) stop("RDS not found for waterfall.")
  
  rb_scen <- read_rds_safe(rds_file)
  rb_frz  <- read_rds_safe(rds_freeze)
  if (inherits(rb_scen, "try-error") || inherits(rb_frz, "try-error")) {
    stop("Corrupt RDS for waterfall: ", rds_file, " or ", rds_freeze)
  }
  
  get_stock <- function(rb) {
    data_list <- list()
    for (sx in c("M", "F")) {
      sex_res <- if (sx == "M") rb$resM else rb$resF
      r_all <- sex_res$inc_fit$rates_all %>% mutate(sex = sx)
      
      # Prevalence stock (from diag)
      z_h <- sex_res$diag$z_prev_hist
      z_f <- sex_res$diag$z_prev_future
      z_all <- bind_rows(z_h, z_f) %>% 
        group_by(period) %>% 
        summarise(
          current_prev  = mean(as.numeric(p_cur), na.rm = TRUE),
          current_q_eff = mean(as.numeric(q_eff), na.rm = TRUE), 
          .groups = "drop"
        ) %>%
        mutate(sex = sx)
      
      res_sx <- sex_res$annual_anchor %>% 
        select(period, deaths_hat) %>%
        mutate(sex = sx) %>%
        left_join(
          r_all %>% group_by(period, sex) %>% summarise(inc_rate = mean(rate_hat, na.rm = TRUE), .groups = "drop"),
          by = c("period", "sex")
        ) %>%
        left_join(z_all, by = c("period", "sex"))
      
      data_list[[sx]] <- res_sx
    }
    bind_rows(data_list)
  }
  
  stock_scen <- get_stock(rb_scen)
  stock_frz  <- get_stock(rb_frz)
  
  # Panel A: Current Prevalence Level
  pA <- ggplot(stock_scen, aes(x = period, y = current_prev * 100, color = sex)) +
    geom_line(linewidth = 1) +
    geom_line(data = stock_frz, aes(y = current_prev * 100), linetype = "dotted", alpha = 0.7) +
    labs(title = "Smoking Prevalence Level (p_curr)", y = "Rate (%)", x = NULL) +
    theme_paper_main(base_size = 9)
    
  # Panel B: Effective Exposure
  pB <- ggplot(stock_scen, aes(x = period, y = current_q_eff * 100, color = sex)) +
    geom_line(linewidth = 1) +
    geom_line(data = stock_frz, aes(y = current_q_eff * 100), linetype = "dotted", alpha = 0.7) +
    labs(title = "Effective Exposure Level (q_eff)", y = "Stock (%)", x = NULL) +
    theme_paper_main(base_size = 9)
  
  # Panel C: Incidence Rate Levels
  pC <- ggplot(stock_scen, aes(x = period, y = inc_rate * 100000, color = sex)) +
    geom_line(linewidth = 1) +
    geom_line(data = stock_frz, aes(y = inc_rate * 100000), linetype = "dotted", alpha = 0.7) +
    labs(title = "Incidence Rate Levels", y = "Rate per 100k", x = NULL) +
    theme_paper_main(base_size = 9) + theme(legend.position = "none")
  
  # Panel D: Total Deaths Levels
  pD <- ggplot(stock_scen, aes(x = period, y = deaths_hat, color = sex)) +
    geom_line(linewidth = 1) +
    geom_line(data = stock_frz, aes(y = deaths_hat), linetype = "dotted", alpha = 0.7) +
    labs(title = "Total Deaths Levels", y = "Total Deaths", x = "Year") +
    theme_paper_main(base_size = 9) + theme(legend.position = "none")
  
  (pA | pB) / (pC | pD) + plot_annotation(title = sprintf("Transmission Waterfall: %s Scenario (Seed %d)", scen, seed)) &
    theme_paper_main(base_size = 10)
}

plot_transmission_map <- function(seed = 4,
                                  dgp = "spec_linear",
                                  sex_lab = "M",
                                  scens = c("up1pc", "freeze", "down1pc", "quit"),
                                  raw_dir = OUT_RAW,
                                  title_suffix = NULL) {
  sex_lab <- match.arg(as.character(sex_lab)[1], c("M", "F"))
  series_levels <- c("Truth", "Informed SBAPC")
  metric_levels <- c(
    "Current smoking prevalence",
    "Effective smoking exposure",
    "Annual incident cases",
    "Annual deaths"
  )

  data_list <- list()
  freeze_hist_inc_by_seed <- NULL
  for (scen in scens) {
    rds_file <- file.path(raw_dir, sprintf("res_%s_s%d_%s.rds", dgp, seed, scen))
    if (!file.exists(rds_file)) {
      message("Transmission map: RDS not found ", rds_file)
      next
    }
    rb <- read_rds_safe(rds_file)
    if (inherits(rb, "try-error")) {
      message("Transmission map: corrupt RDS ", rds_file)
      next
    }
    sex_res <- if (identical(sex_lab, "M")) rb$resM else rb$resF

    sim_truth <- simulate_PIM_data(
      cause_id = CAUSE_ID,
      seed = seed,
      dgp = dgp,
      scenario_name = scen,
      sexes = sex_lab
    )

    truth_stock <- sim_truth$z_scen_true %>%
      dplyr::filter(as.character(sex) == sex_lab) %>%
      dplyr::group_by(period) %>%
      dplyr::summarise(
        `Current smoking prevalence` = mean(as.numeric(p_curr), na.rm = TRUE) * 100,
        `Effective smoking exposure` = mean(as.numeric(q_eff), na.rm = TRUE) * 100,
        .groups = "drop"
      ) %>%
      tidyr::pivot_longer(-period, names_to = "metric", values_to = "value") %>%
      dplyr::mutate(series = "Truth")

    est_stock <- dplyr::bind_rows(
      tibble::as_tibble(sex_res$diag$z_prev_hist %||% tibble::tibble()),
      tibble::as_tibble(sex_res$diag$z_prev_future %||% tibble::tibble())
    ) %>%
      dplyr::filter(as.character(sex) == sex_lab) %>%
      dplyr::group_by(period) %>%
      dplyr::summarise(
        `Current smoking prevalence` = mean(as.numeric(p_cur), na.rm = TRUE) * 100,
        `Effective smoking exposure` = mean(as.numeric(q_eff), na.rm = TRUE) * 100,
        .groups = "drop"
      ) %>%
      tidyr::pivot_longer(-period, names_to = "metric", values_to = "value") %>%
      dplyr::mutate(series = "Informed SBAPC")

    exposure <- suppressWarnings(as.numeric(rb$meta$args$exposure %||% sim_truth$meta$exposure %||% 100000))[1]
    if (!is.finite(exposure) || exposure <= 0) exposure <- 100000

    truth_inc <- rb$inc_truth_grid %>%
      dplyr::filter(as.character(sex) == sex_lab) %>%
      dplyr::group_by(period) %>%
      dplyr::summarise(value = sum(as.numeric(rateI_scen_true) * exposure, na.rm = TRUE), .groups = "drop") %>%
      dplyr::mutate(metric = "Annual incident cases", series = "Truth")

    est_inc <- sex_res$inc_annual_cond %>%
      dplyr::transmute(period = as.integer(period), value = as.numeric(cases_hat),
                       metric = "Annual incident cases", series = "Informed SBAPC")
    if (!any(est_inc$period <= 2022, na.rm = TRUE) && !identical(scen, "freeze")) {
      if (is.null(freeze_hist_inc_by_seed)) {
        freeze_file <- file.path(raw_dir, sprintf("res_%s_s%d_freeze.rds", dgp, seed))
        rb_freeze <- read_rds_safe(freeze_file)
        if (!inherits(rb_freeze, "try-error")) {
          sex_freeze <- if (identical(sex_lab, "M")) rb_freeze$resM else rb_freeze$resF
          freeze_hist_inc_by_seed <- sex_freeze$inc_annual_cond %>%
            dplyr::filter(as.integer(period) <= 2022) %>%
            dplyr::transmute(period = as.integer(period), value = as.numeric(cases_hat),
                             metric = "Annual incident cases", series = "Informed SBAPC")
        } else {
          freeze_hist_inc_by_seed <- tibble::tibble()
        }
      }
      est_inc <- dplyr::bind_rows(freeze_hist_inc_by_seed, est_inc) %>%
        dplyr::arrange(period)
    }

    truth_mort <- rb$mort_truth_grid %>%
      dplyr::filter(as.character(sex) == sex_lab) %>%
      dplyr::group_by(period) %>%
      dplyr::summarise(value = sum(as.numeric(mort_deaths_scen_true), na.rm = TRUE), .groups = "drop") %>%
      dplyr::mutate(metric = "Annual deaths", series = "Truth")

    est_mort <- sex_res$annual_anchor %>%
      dplyr::transmute(period = as.integer(period), value = as.numeric(deaths_hat),
                       metric = "Annual deaths", series = "Informed SBAPC")

    data_list[[scen]] <- dplyr::bind_rows(truth_stock, est_stock, truth_inc, est_inc, truth_mort, est_mort) %>%
      dplyr::mutate(
        scenario = SCEN_LABELS[[scen]],
        metric = factor(metric, levels = metric_levels),
        series = factor(series, levels = series_levels)
      )
  }

  df <- dplyr::bind_rows(data_list) %>%
    dplyr::filter(period >= 1998) %>%
    dplyr::mutate(scenario = factor(scenario, levels = SCEN_LABELS[scens]))

  ggplot(df, aes(x = period, y = value, color = series, linetype = series)) +
    facet_grid(metric ~ scenario, scales = "free_y") +
    geom_vline(xintercept = 2022, linetype = "dotted", color = "gray60", linewidth = 0.35) +
    geom_line(linewidth = 0.85, na.rm = TRUE) +
    scale_color_manual(values = c("Truth" = "black", "Informed SBAPC" = "#D32F2F"), breaks = series_levels) +
    scale_linetype_manual(values = c("Truth" = "dashed", "Informed SBAPC" = "solid"), breaks = series_levels) +
    labs(
      title = paste0(
        sprintf("Smoking-to-Mortality Transmission Map: %s (Seed %d)",
                ifelse(sex_lab == "M", "Males", "Females"), seed),
        if (!is.null(title_suffix) && nzchar(title_suffix)) paste0(" - ", title_suffix) else ""
      ),
      subtitle = "Data-based pathway comparison from prevalence state to projected deaths",
      x = "Year", y = NULL, color = "Series", linetype = "Series"
    ) +
    theme_paper_main(base_size = 10) +
    theme(
      legend.position = "bottom",
      strip.background = element_rect(fill = "gray95"),
      panel.spacing = grid::unit(0.6, "lines")
    )
}

plot_transmission_map_support_compare <- function(seed = 4,
                                                  dgp = "spec_linear",
                                                  sex_lab = "M",
                                                  scens = c("up1pc", "freeze", "down1pc", "quit"),
                                                  realistic_raw_dir = OUT_RAW,
                                                  oracle_raw_dir = OUT_RAW_ORACLE,
                                                  realistic_label = "Window-limited SBAPC",
                                                  oracle_label = "Full-support SBAPC") {
  sex_lab <- match.arg(as.character(sex_lab)[1], c("M", "F"))
  metric_levels <- c(
    "Current smoking prevalence",
    "Effective smoking exposure",
    "Annual incident cases",
    "Annual deaths"
  )
  series_levels <- c("Truth", realistic_label, oracle_label)

  g_real <- plot_transmission_map(
    seed = seed, dgp = dgp, sex_lab = sex_lab, scens = scens,
    raw_dir = realistic_raw_dir
  )
  g_oracle <- plot_transmission_map(
    seed = seed, dgp = dgp, sex_lab = sex_lab, scens = scens,
    raw_dir = oracle_raw_dir
  )

  df_real <- tibble::as_tibble(g_real$data)
  df_oracle <- tibble::as_tibble(g_oracle$data)

  truth_df <- df_oracle %>%
    dplyr::filter(as.character(series) == "Truth") %>%
    dplyr::mutate(series = "Truth")
  real_df <- df_real %>%
    dplyr::filter(as.character(series) == "Informed SBAPC") %>%
    dplyr::mutate(series = realistic_label)
  oracle_df <- df_oracle %>%
    dplyr::filter(as.character(series) == "Informed SBAPC") %>%
    dplyr::mutate(series = oracle_label)

  df <- dplyr::bind_rows(truth_df, real_df, oracle_df) %>%
    dplyr::mutate(
      metric = factor(as.character(metric), levels = metric_levels),
      series = factor(as.character(series), levels = series_levels),
      scenario = factor(as.character(scenario), levels = SCEN_LABELS[scens])
    )

  ggplot(df, aes(x = period, y = value, color = series, linetype = series)) +
    facet_grid(metric ~ scenario, scales = "free_y") +
    geom_vline(xintercept = 2022, linetype = "dotted", color = "gray60", linewidth = 0.35) +
    geom_line(linewidth = 0.85, na.rm = TRUE) +
    scale_color_manual(
      values = stats::setNames(c("black", "#D32F2F", "#1565C0"), series_levels),
      breaks = series_levels
    ) +
    scale_linetype_manual(
      values = stats::setNames(c("dashed", "solid", "longdash"), series_levels),
      breaks = series_levels
    ) +
    labs(
      title = sprintf("Smoking-to-Mortality Transmission Map: %s (Seed %d)",
                      ifelse(sex_lab == "M", "Males", "Females"), seed),
      subtitle = "Truth versus full-support and window-limited estimators",
      x = "Year", y = NULL, color = "Series", linetype = "Series"
    ) +
    theme_paper_main(base_size = 10) +
    theme(
      legend.position = "bottom",
      strip.background = element_rect(fill = "gray95"),
      panel.spacing = grid::unit(0.6, "lines")
    )
}

generate_support_transmission_maps <- function(seed = 4,
                                               dgp = "spec_linear",
                                               scens = CANONICAL_SCENS,
                                               realistic_raw_dir = OUT_RAW,
                                               oracle_raw_dir = OUT_RAW_ORACLE,
                                               force_oracle = FALSE) {
  real_files_ok <- all(vapply(scens, function(sc) {
    file.exists(file.path(realistic_raw_dir, sprintf("res_%s_s%d_%s.rds", dgp, seed, sc)))
  }, logical(1)))
  if (!real_files_ok) {
    stop("Missing window-limited RDS files for support transmission map in: ", realistic_raw_dir)
  }

  oracle_files_ok <- all(vapply(scens, function(sc) {
    file.exists(file.path(oracle_raw_dir, sprintf("res_%s_s%d_%s.rds", dgp, seed, sc)))
  }, logical(1)))
  if (!oracle_files_ok || isTRUE(force_oracle)) {
    run_single_seed_replication(
      seed = seed,
      dgp = dgp,
      scens = scens,
      force_rerun = force_oracle,
      information_set = "oracle",
      raw_dir = oracle_raw_dir
    )
  }

  for (sx in c("M", "F")) {
    g <- plot_transmission_map_support_compare(
      seed = seed,
      dgp = dgp,
      sex_lab = sx,
      scens = scens,
      realistic_raw_dir = realistic_raw_dir,
      oracle_raw_dir = oracle_raw_dir,
      realistic_label = "Window-limited SBAPC",
      oracle_label = "Full-support SBAPC"
    )
    save_paper_plot(
      g,
      file.path(OUT_SEC4, sprintf("fig_transmission_map_support_compare_seed%d_%s", seed, sx)),
      width = 13,
      height = 8,
      bg = "white"
    )
  }

  invisible(TRUE)
}

plot_reliability_calibration <- function(data) {
  # data$inc contains the errors per period/sex/seed/dgp
  
  df <- data$inc %>%
    dplyr::filter(period > 2022) %>%
    dplyr::mutate(horizon = period - 2022)
  
  # Calculate endogenous horizons from support_frac thresholds (0.5 and 0.33)
  horizons_info <- df %>%
    group_by(horizon) %>%
    summarise(avg_support = mean(support_frac, na.rm = TRUE), .groups = "drop") %>%
    arrange(horizon)
  
  h_caution <- horizons_info$horizon[which(horizons_info$avg_support < 0.50)[1]]
  h_risky   <- horizons_info$horizon[which(horizons_info$avg_support < 0.33)[1]]
  
  vlines <- c()
  if (!is.na(h_caution)) vlines <- c(vlines, h_caution)
  if (!is.na(h_risky))   vlines <- c(vlines, h_risky)

  # Background sectors
  rects <- tibble::tibble(
    xmin = c(0, h_caution, h_risky),
    xmax = c(h_caution, h_risky, max(df$horizon, na.rm=TRUE)),
    fill = c("#D1E5D1", "#FFF9C4", "#FFCDD2"), # Green, Yellow, Red light
    label = c("Credible", "Caution", "Risky")
  ) %>% filter(is.finite(xmin), is.finite(xmax))

  g <- ggplot() +
    geom_rect(data = rects, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill), alpha = 0.5) +
    stat_summary(data = df, aes(x = horizon, y = abs(rel_error) * 100), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.2, fill = "blue") +
    stat_summary(data = df, aes(x = horizon, y = abs(rel_error) * 100), fun = "mean", geom = "line", linewidth = 1, color = "blue") +
    scale_fill_identity()
    
  if (length(vlines) > 0) {
    g <- g + geom_vline(xintercept = vlines, linetype = "dashed", color = "gray40")
    # Add labels for categories
    for(i in seq_along(rects$label)) {
      g <- g + annotate("text", x = (rects$xmin[i] + rects$xmax[i])/2, y = Inf, label = rects$label[i], vjust = 2, size = 3.5, fontface = "italic", color = "gray30")
    }
    # Add year labels on top of vlines
    for(hv in vlines) {
       g <- g + annotate("label", x = hv, y = Inf, label = sprintf("Year %d", hv + 2022), vjust = 0.8, size = 3)
    }
  }

  g <- g + labs(title = "Reliability Calibration: Error vs. Horizon",
         y = "Mean Absolute Relative Error (%)", x = "Projection Horizon (Years)") +
    theme_paper_main(base_size = 11)
  
  return(g)
}


# =============================================================
# 4. MAIN OUTPUT GENERATOR
# =============================================================

replicate_main_paper <- function() {
  # 1. Scenario Atlas (By Sex)
  for (sx in c("M", "F")) {
    g_atlas <- plot_scenario_atlas_by_sex(seed = 4, sex_lab = sx)
    save_paper_plot(g_atlas, file.path(OUT_SEC4, sprintf("fig_scenario_atlas_seed4_%s", sx)), width = 14, height = 5, bg = "white")
  }
  
  # 2. Waterfall
  g2 <- plot_transmission_waterfall(seed = 4, dgp = "spec_linear", scen = "quit")
  save_paper_plot(g2, file.path(OUT_SEC4, "fig_waterfall_seed4"), width = 8, height = 10, bg = "white")
  
  # 3. Scenario Sensitivity (New)
  g_sens <- plot_scenario_sensitivity_informed(seed = 4, dgp = "spec_linear")
  save_paper_plot(g_sens, file.path(OUT_SEC4, "fig_sensitivity_seed4"), width = 10, height = 6, bg = "white")
  
  # 4. Transmission map
  g_map <- plot_transmission_map(seed = 4, dgp = "spec_linear", sex_lab = "M")
  save_paper_plot(g_map, file.path(OUT_SEC4, "fig_transmission_map_seed4_M"), width = 13, height = 9, bg = "white")
  generate_support_transmission_maps(seed = 4, dgp = "spec_linear", force_oracle = FALSE)

  # 5. Bias Table
  data <- extract_all_metrics()
  export_latex_bias_summary(data$metrics, file.path(OUT_SEC4, "tab_bias_summary.tex"))
  write_csv(data$metrics, file.path(OUT_RAW, "all_metrics.csv"))
  
  # 6. Reliability Plot
  g3 <- plot_reliability_calibration(data)
  save_paper_plot(g3, file.path(OUT_SEC4, "fig_reliability_calibration"), width = 10, height = 6, bg = "white")
  
  # 7. Support Summary
  write_csv(data$support, file.path(OUT_SEC4, "support_summary.csv"))
  
  message("\nMain paper replication files generated in: ", OUT_SEC4)
}

# =============================================================
# 5. APPENDIX C GENERATOR
# =============================================================

generate_appendix_c <- function() {
  data <- extract_all_metrics()
  
  # 1. Distribution of Bias (Boxplots)
  plot_df <- data$metrics %>%
    dplyr::mutate(scenario_f = factor(scenario, levels = c("up1pc", "freeze", "down1pc", "quit"), labels = SCEN_LABELS))
  scen_colors_by_label <- stats::setNames(unname(SCEN_COLORS[names(SCEN_LABELS)]), unname(SCEN_LABELS))
  
  g_bias <- ggplot(plot_df, aes(x = scenario_f, y = proj_bias, fill = scenario_f)) +
    geom_boxplot(alpha = 0.7) +
    facet_wrap(~dgp + sex) +
    scale_fill_manual(values = scen_colors_by_label) +
    labs(title = "Appendix C: Distribution of Projection Bias across Seeds",
         y = "Projection Bias (%)", x = "Scenario") +
    theme_minimal() + theme(legend.position = "none")
  save_paper_plot(g_bias, file.path(OUT_APPENDIX, "fig_bias_distributions"), width = 10, height = 7, bg = "white")
  
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
    save_paper_plot(g_case, file.path(OUT_APPENDIX, sprintf("fig_case_study_%s_s%d", tolower(lbl), s)), width = 10, height = 6, bg = "white")
  }
  
  # 3. Full Detailed Table (CSV)
  write_csv(data$metrics, file.path(OUT_APPENDIX, "full_simulation_matrix.csv"))
  
  message("\nAppendix C replication files generated in: ", OUT_APPENDIX)
}

# =============================================================
# 6. ORCHESTRATOR
# =============================================================

replicate_all_simulations <- function(n_cores = 6, force_rerun = TRUE) {
  message("STARTING FULL REPLICATION WORKFLOW...")
  
  # 1. Run Simulations (FORCE RERUN to overwrite old files)
  # Using 6 cores to avoid memory allocation errors (INLA is memory intensive)
  run_simulation_replication(n_cores = n_cores, force_rerun = force_rerun)
  
  # 2. Section 4
  replicate_main_paper()
  
  # 3. Appendix C
  generate_appendix_c()
  
  message("\nALL REPLICATION TASKS COMPLETED SUCCESSFULLY.")
}
