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
.n_seeds_env <- suppressWarnings(as.integer(Sys.getenv("BAPC_N_SEEDS", "50")))
if (!is.finite(.n_seeds_env) || .n_seeds_env < 1L) .n_seeds_env <- 50L
CANONICAL_SEEDS <- seq_len(.n_seeds_env)
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

SEX_LABELS <- c(
  "M" = "Male",
  "F" = "Female",
  "T" = "Total",
  "Total" = "Total"
)

SEX_COLORS <- c(
  "Male" = "#1F77B4",
  "Female" = "#D32F2F"
)

DGP_LABELS <- c(
  "spec_linear" = "Well-specified design",
  "misspec_tanh" = "Misspecified transmission design"
)

sex_public_label <- function(x) {
  x_chr <- as.character(x)
  out <- unname(SEX_LABELS[x_chr])
  out[is.na(out)] <- x_chr[is.na(out)]
  out
}

dgp_public_label <- function(x) {
  x_chr <- as.character(x)
  out <- unname(DGP_LABELS[x_chr])
  out[is.na(out)] <- x_chr[is.na(out)]
  out
}

MODEL_LABELS <- c(
  truth = "Truth",
  sbapc = "SBAPC",
  sbapc_no_prev = "Incidence-anchored SBAPC",
  bapc = "BAPC benchmark"
)

MODEL_COLORS <- stats::setNames(
  c("black", "#D32F2F", "#EF6C00", "#1976D2"),
  unname(MODEL_LABELS[c("truth", "sbapc", "sbapc_no_prev", "bapc")])
)

MODEL_LINETYPES <- stats::setNames(
  c("dashed", "solid", "longdash", "dotted"),
  unname(MODEL_LABELS[c("truth", "sbapc", "sbapc_no_prev", "bapc")])
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
  if (!nzchar(val)) val <- getOption("BAPC_FIG_FORMAT", BAPC_FIG_FORMAT %||% "both")
  match.arg(as.character(val)[1], c("svg", "pdf", "png", "both"))
}

figure_exts <- function(format = FIG_FORMAT) {
  if (identical(format, "both")) c("svg", "pdf") else format
}

seed_count_label <- function(seeds = CANONICAL_SEEDS) {
  sprintf("%d seed%s", length(unique(seeds)), if (length(unique(seeds)) == 1L) "" else "s")
}

cleanup_inla_temp <- function(tmpdir = BAPC_PATHS$inla_tmp) {
  tmpdir <- normalizePath(tmpdir, winslash = "/", mustWork = FALSE)
  if (!grepl("(^|/)tmp_inla$", tmpdir, ignore.case = TRUE)) {
    stop("Refusing to clean unexpected INLA temp directory: ", tmpdir)
  }
  targets <- list.files(tmpdir, pattern = "^inla\\.model-", full.names = TRUE, recursive = FALSE)
  targets <- targets[dir.exists(targets)]
  if (length(targets)) unlink(targets, recursive = TRUE, force = TRUE)
  invisible(length(targets))
}

save_paper_plot <- function(plot, path_no_ext, width, height, bg = "white", format = FIG_FORMAT, ...) {
  for (ext in figure_exts(format)) {
    args <- list(
      filename = paste0(path_no_ext, ".", ext),
      plot = plot,
      width = width,
      height = height,
      bg = bg,
      ...
    )
    if (identical(ext, "pdf")) args$device <- grDevices::cairo_pdf
    do.call(ggplot2::ggsave, args)
  }
  invisible(path_no_ext)
}

PAPER_EXPORT_PROFILE <- {
  val <- Sys.getenv("BAPC_PAPER_EXPORT_PROFILE", unset = "")
  if (!nzchar(val)) val <- getOption("BAPC_PAPER_EXPORT_PROFILE", "manuscript")
  match.arg(as.character(val)[1], c("manuscript", "legacy"))
}

PAPER_FIG_SPECS <- list(
  manuscript = list(
    scenario_atlas = list(width = 7.2, height = 3.3, base_size = 9.8),
    waterfall = list(width = 7.0, height = 8.7, base_size = 10.2),
    sensitivity = list(width = 6.8, height = 4.0, base_size = 10.0),
    transmission_map = list(width = 7.2, height = 5.6, base_size = 8.8),
    transmission_support = list(width = 7.2, height = 5.3, base_size = 8.8),
    scenario_effect = list(width = 7.2, height = 3.8, base_size = 9.6),
    scenario_effect_bysex = list(width = 7.2, height = 5.6, base_size = 8.9),
    support_window = list(width = 7.2, height = 3.8, base_size = 9.4),
    misspecification = list(width = 7.2, height = 5.2, base_size = 8.8),
    reliability = list(width = 6.8, height = 4.1, base_size = 10.0),
    bias_distribution = list(width = 6.8, height = 4.8, base_size = 9.8),
    case_study = list(width = 6.8, height = 4.2, base_size = 9.8)
  ),
  legacy = list(
    scenario_atlas = list(width = 14, height = 5, base_size = 12),
    waterfall = list(width = 8, height = 10, base_size = 9),
    sensitivity = list(width = 10, height = 6, base_size = 10.5),
    transmission_map = list(width = 13, height = 9, base_size = 10),
    transmission_support = list(width = 13, height = 8, base_size = 10),
    scenario_effect = list(width = 12, height = 5.8, base_size = 10),
    scenario_effect_bysex = list(width = 12, height = 7.5, base_size = 10),
    support_window = list(width = 12, height = 5.8, base_size = 10),
    misspecification = list(width = 12, height = 7.2, base_size = 10),
    reliability = list(width = 10, height = 6, base_size = 11),
    bias_distribution = list(width = 10, height = 7, base_size = 10),
    case_study = list(width = 10, height = 6, base_size = 11)
  )
)

paper_fig_spec <- function(key, profile = PAPER_EXPORT_PROFILE) {
  spec <- PAPER_FIG_SPECS[[profile]][[key]]
  if (is.null(spec)) stop("Unknown paper figure spec: ", key)
  spec
}

paper_fig_base_size <- function(key) paper_fig_spec(key)$base_size

save_profiled_plot <- function(plot, path_no_ext, key, bg = "white", ...) {
  spec <- paper_fig_spec(key)
  save_paper_plot(
    plot = plot,
    path_no_ext = path_no_ext,
    width = spec$width,
    height = spec$height,
    bg = bg,
    ...
  )
}

figure_file_names <- function(stem, format = FIG_FORMAT) {
  paste0(stem, ".", figure_exts(format))
}

write_figure_titles_notes <- function(section = c("section4", "appendixC"), case_seeds = NULL) {
  section <- match.arg(section)
  out_dir <- if (identical(section, "section4")) OUT_SEC4 else OUT_APPENDIX

  entry <- function(stem, title, note) {
    c(
      paste0("## ", stem),
      paste0("Files: ", paste(figure_file_names(stem), collapse = ", ")),
      paste0("Title: ", title),
      paste0("Note: ", note, " Source: Own elaboration."),
      ""
    )
  }

  if (identical(section, "section4")) {
    lines <- c(
      "# Figure Titles and Notes: Section 4",
      "",
      "Use these titles and notes in the LaTeX manuscript. The graphics themselves intentionally omit global titles, subtitles, and notes.",
      "",
      entry(
        "fig_scenario_atlas_seed4_M",
        "Scenario atlas for male mortality projections",
        "Illustrative simulation draw. Panels show smoking-prevalence scenarios. The vertical reference line marks the last historical year, 2022. Lines compare simulated truth, SBAPC, incidence-anchored SBAPC, and the BAPC benchmark."
      ),
      entry(
        "fig_scenario_atlas_seed4_F",
        "Scenario atlas for female mortality projections",
        "Illustrative simulation draw. Panels show smoking-prevalence scenarios. The vertical reference line marks the last historical year, 2022. Lines compare simulated truth, SBAPC, incidence-anchored SBAPC, and the BAPC benchmark."
      ),
      entry(
        "fig_waterfall_seed4",
        "Transmission pathway under the quit scenario",
        "Illustrative simulation draw. Solid lines show the quit scenario and dotted lines show the frozen-prevalence baseline. Panels trace current smoking prevalence, effective smoking exposure, incidence rates, and annual deaths by sex."
      ),
      entry(
        "fig_sensitivity_seed4",
        "Projected mortality sensitivity to smoking-prevalence scenarios",
        "Illustrative simulation draw. Total annual deaths are shown under the four prevalence scenarios. The vertical reference line marks 2022."
      ),
      entry(
        "fig_transmission_map_seed4_M",
        "Smoking-to-mortality transmission map for male projections",
        "Illustrative male simulation draw. Rows trace prevalence, effective exposure, annual incident cases, and annual deaths. Columns show prevalence scenarios. The vertical reference line marks 2022."
      ),
      entry(
        "fig_transmission_map_support_compare_seed4_M",
        "Observed-window and full-support transmission map for male projections",
        "Illustrative male simulation draw. The figure compares simulated truth, Observed-window SBAPC, and Full-support SBAPC across the smoking-to-mortality pathway."
      ),
      entry(
        "fig_scenario_effect_recovery",
        "Mortality scenario-effect recovery",
        "Scenario effects are annual deaths relative to the frozen-prevalence baseline, aggregated across sexes and simulation seeds. Lines show median effects across seeds; ribbons show the 10th-90th percentile range for simulated truth and SBAPC. The BAPC benchmark is scenario-blind and therefore has zero scenario response by construction. Background shading is used only because the horizon-boundary audit found common support-region boundaries for the plotted aggregation."
      )
    )
  } else {
    lines <- c(
      "# Figure Titles and Notes: Appendix C",
      "",
      "Use these titles and notes in the LaTeX supplement. The graphics themselves intentionally omit global titles, subtitles, and notes.",
      "",
      entry(
        "fig_bias_distributions",
        "Distribution of projection bias across simulation seeds",
        "Boxplots summarize projection bias by scenario and sex across simulation seeds."
      ),
      entry(
        "fig_scenario_effect_recovery_bysex",
        "Mortality scenario-effect recovery by sex",
        "Extended diagnostic by sex including the incidence-anchored SBAPC decomposition variant."
      ),
      entry(
        "fig_support_window_comparison",
        "Observed-window and full-support scenario-effect recovery",
        "Appendix diagnostic comparing simulated truth, observed-window SBAPC, and full-support SBAPC for annual mortality scenario effects relative to the frozen-prevalence baseline. Full-support SBAPC is an oracle-style diagnostic, not a feasible empirical estimator."
      ),
      entry(
        "fig_misspecification_scenario_recovery",
        "Scenario-effect recovery under transmission-rule misspecification",
        "Appendix robustness diagnostic comparing the well-specified design with a misspecified monotone transmission design. Panels show annual mortality scenario effects relative to the frozen-prevalence baseline."
      ),
      entry(
        "fig_reliability_calibration",
        "Projection reliability calibration by support horizon",
        "Calibration is treated as an uncertainty-summary diagnostic. The central validation target for Section 4 is recovery of mortality scenario effects, not calibration of predictive summaries."
      )
    )

    if (!is.null(case_seeds) && nrow(case_seeds)) {
      for (i in seq_len(nrow(case_seeds))) {
        lbl <- as.character(case_seeds$label[i])
        seed_i <- as.integer(case_seeds$seed[i])
        stem <- sprintf("fig_case_study_%s_s%d", tolower(lbl), seed_i)
        lines <- c(lines, entry(
          stem,
          sprintf("%s quit-scenario case study", lbl),
          "Illustrative quit-scenario trajectory diagnostic selected by absolute projection bias among male simulations."
        ))
      }
    }
  }

  writeLines(enc2utf8(lines), file.path(out_dir, "figure_titles_notes.md"), useBytes = TRUE)
  invisible(file.path(out_dir, "figure_titles_notes.md"))
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

metrics_cache_file <- function(raw_dir = OUT_RAW,
                               dgps = CANONICAL_DGPS,
                               scens = CANONICAL_SCENS,
                               suffix = NULL) {
  if (!is.null(suffix) && nzchar(as.character(suffix)[1])) {
    return(file.path(raw_dir, paste0("all_extracted_data_", suffix, ".rds")))
  }
  is_default <- identical(normalizePath(raw_dir, winslash = "/", mustWork = FALSE),
                          normalizePath(OUT_RAW, winslash = "/", mustWork = FALSE)) &&
    identical(sort(as.character(dgps)), sort(as.character(CANONICAL_DGPS))) &&
    identical(sort(as.character(scens)), sort(as.character(CANONICAL_SCENS)))
  if (is_default) return(file.path(raw_dir, "all_extracted_data.rds"))
  safe <- paste(c(sort(as.character(dgps)), sort(as.character(scens))), collapse = "_")
  safe <- gsub("[^A-Za-z0-9]+", "_", safe)
  file.path(raw_dir, paste0("all_extracted_data_", safe, ".rds"))
}

extract_all_metrics <- function(seeds = CANONICAL_SEEDS,
                                dgps = CANONICAL_DGPS,
                                scens = CANONICAL_SCENS,
                                force_refresh = FALSE,
                                raw_dir = OUT_RAW,
                                cache_file = NULL,
                                cache_suffix = NULL) {
  if (is.null(cache_file)) cache_file <- metrics_cache_file(raw_dir, dgps, scens, cache_suffix)
  if (file.exists(cache_file) && !isTRUE(force_refresh)) {
    message("Loading cached metrics from: ", cache_file)
    cached <- readRDS(cache_file)
    cached_scens <- sort(unique(as.character(cached$metrics$scenario %||% character(0))))
    requested_scens <- sort(unique(as.character(scens)))
    cached_dgps <- sort(unique(as.character(cached$metrics$dgp %||% character(0))))
    requested_dgps <- sort(unique(as.character(dgps)))
    cached_seeds <- sort(unique(as.integer(cached$metrics$seed %||% integer(0))))
    requested_seeds <- sort(unique(as.integer(seeds)))
    if (identical(cached_scens, requested_scens) &&
        identical(cached_dgps, requested_dgps) &&
        identical(cached_seeds, requested_seeds)) {
      return(cached)
    }
    message("Cached metrics use a different seed/scenario/DGP set; rebuilding extraction cache.")
  }
  
  metrics_list <- list()
  deltas_list  <- list()
  support_list <- list()
  inc_list     <- list()
  mort_list    <- list()
  
  for (seed in seeds) {
    for (dgp in dgps) {
      # First, get freeze mort for deltas
      freeze_rds <- file.path(raw_dir, sprintf("res_%s_s%d_freeze.rds", dgp, seed))
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
        rds_file <- file.path(raw_dir, sprintf("res_%s_s%d_%s.rds", dgp, seed, scen))
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

plot_deconstruction_figure <- function(seed = 4, dgp = "spec_linear", scen = "quit",
                                       base_size = paper_fig_base_size("case_study")) {
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
      `Truth` = deaths_true, 
      `SBAPC` = deaths_hat, 
      `Incidence-anchored SBAPC` = deaths_noP, 
      `BAPC benchmark` = deaths_bapc
    ) %>%
    tidyr::pivot_longer(
      cols = c(`Truth`, `SBAPC`, `Incidence-anchored SBAPC`, `BAPC benchmark`), 
      names_to = "Series", values_to = "Deaths"
    ) %>%
    dplyr::mutate(sex = factor(sex_public_label(sex), levels = c("Male", "Female")))
  
  decomp_levels <- unname(MODEL_LABELS[c("truth", "sbapc", "sbapc_no_prev", "bapc")])
  plot_df$Series <- factor(plot_df$Series, levels = decomp_levels)
  
  last_hist <- rb$meta$last_hist %||% 2022
  
  g <- ggplot(plot_df, aes(x = period, y = Deaths, color = Series, linetype = Series)) +
    geom_line(linewidth = 1) +
    geom_vline(xintercept = last_hist, linetype = "dotted", color = "gray50") +
    facet_wrap(~sex, scales = "free_y") +
    scale_color_manual(values = MODEL_COLORS, breaks = decomp_levels) +
    scale_linetype_manual(values = MODEL_LINETYPES, breaks = decomp_levels) +
    labs(y = "Annual deaths", x = "Year", color = "Series", linetype = "Series") +
    theme_paper_main(base_size = base_size) +
    theme(legend.position = "bottom")
  
  return(g)
}

plot_scenario_atlas_by_sex <- function(seed = 4, sex_lab = "M",
                                       base_size = paper_fig_base_size("scenario_atlas")) {
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
    
    # 2. SBAPC
    df_inf <- data.frame(period = sex_res$annual_anchor$period, deaths = sex_res$annual_anchor$deaths_hat, model = MODEL_LABELS[["sbapc"]])
    
    # 3. Incidence-anchored SBAPC
    df_uninf <- data.frame(period = sex_res$annual_anchor_noP$period, deaths = sex_res$annual_anchor_noP$deaths_hat, model = MODEL_LABELS[["sbapc_no_prev"]])
    
    # 4. BAPC benchmark
    df_bapc <- data.frame(period = sex_res$annual_bapc$period, deaths = sex_res$annual_bapc$deaths_hat, model = MODEL_LABELS[["bapc"]])
    
    combined <- bind_rows(df_truth, df_inf, df_uninf, df_bapc) %>%
      mutate(scenario = SCEN_LABELS[[scen]])
    
    all_data[[scen]] <- combined
  }
  
  df_final <- bind_rows(all_data) %>% filter(period > 1995)
  df_final$scenario <- factor(df_final$scenario, levels = SCEN_LABELS[scens])
  model_levels <- unname(MODEL_LABELS[c("truth", "sbapc", "sbapc_no_prev", "bapc")])
  df_final$model <- factor(df_final$model, levels = model_levels)
  
  # Colors and Types
  pal <- MODEL_COLORS[model_levels]
  types <- MODEL_LINETYPES[model_levels]
  
  g <- ggplot(df_final, aes(x = period, y = deaths, color = model, linetype = model)) +
    facet_wrap(~scenario, nrow = 1, scales = "fixed") +
    geom_vline(xintercept = 2022, linetype = "dotted", color = "gray60") +
    geom_line(linewidth = 0.9) +
    scale_color_manual(values = pal, breaks = model_levels) +
    scale_linetype_manual(values = types, breaks = model_levels) +
    labs(y = "Annual deaths", x = "Year", color = "Model", linetype = "Model") +
    theme_paper_main(base_size = base_size) +
    theme(legend.position = "bottom", strip.background = element_rect(fill = "gray95"))
  
  return(g)
}

# =============================================================

plot_scenario_sensitivity_informed <- function(seed = 4, dgp = "spec_linear",
                                               base_size = paper_fig_base_size("sensitivity")) {
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
    labs(y = "Total annual deaths", x = "Year", color = "Scenario") +
    theme_paper_main(base_size = base_size)
}

plot_transmission_waterfall <- function(seed = 4, dgp = "spec_linear", scen = "quit",
                                        base_size = paper_fig_base_size("waterfall")) {
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
      r_all <- sex_res$inc_fit$rates_all %>% mutate(sex = sex_public_label(sx))
      
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
        mutate(sex = sex_public_label(sx))
      
      res_sx <- sex_res$annual_anchor %>% 
        select(period, deaths_hat) %>%
        mutate(sex = sex_public_label(sx)) %>%
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
  stock_scen$sex <- factor(stock_scen$sex, levels = c("Male", "Female"))
  stock_frz$sex <- factor(stock_frz$sex, levels = c("Male", "Female"))
  
  # Panel A: Current Prevalence Level
  pA <- ggplot(stock_scen, aes(x = period, y = current_prev * 100, color = sex)) +
    geom_line(linewidth = 1) +
    geom_line(data = stock_frz, aes(y = current_prev * 100), linetype = "dotted", alpha = 0.7) +
    scale_color_manual(values = SEX_COLORS, drop = FALSE) +
    labs(title = "Smoking prevalence", y = "Percent", x = NULL, color = "Sex") +
    theme_paper_main(base_size = base_size)
    
  # Panel B: Effective Exposure
  pB <- ggplot(stock_scen, aes(x = period, y = current_q_eff * 100, color = sex)) +
    geom_line(linewidth = 1) +
    geom_line(data = stock_frz, aes(y = current_q_eff * 100), linetype = "dotted", alpha = 0.7) +
    scale_color_manual(values = SEX_COLORS, drop = FALSE) +
    labs(title = "Effective smoking exposure", y = "Percent", x = NULL, color = "Sex") +
    theme_paper_main(base_size = base_size)
  
  # Panel C: Incidence Rate Levels
  pC <- ggplot(stock_scen, aes(x = period, y = inc_rate * 100000, color = sex)) +
    geom_line(linewidth = 1) +
    geom_line(data = stock_frz, aes(y = inc_rate * 100000), linetype = "dotted", alpha = 0.7) +
    scale_color_manual(values = SEX_COLORS, drop = FALSE) +
    labs(title = "Incidence rate", y = "Rate per 100,000", x = NULL, color = "Sex") +
    theme_paper_main(base_size = base_size)
  
  # Panel D: Total Deaths Levels
  pD <- ggplot(stock_scen, aes(x = period, y = deaths_hat, color = sex)) +
    geom_line(linewidth = 1) +
    geom_line(data = stock_frz, aes(y = deaths_hat), linetype = "dotted", alpha = 0.7) +
    scale_color_manual(values = SEX_COLORS, drop = FALSE) +
    labs(title = "Annual deaths", y = "Annual deaths", x = "Year", color = "Sex") +
    theme_paper_main(base_size = base_size)
  
  ((pA | pB) / (pC | pD)) +
    patchwork::plot_layout(guides = "collect") &
    theme_paper_main(base_size = base_size) &
    theme(legend.position = "bottom")
}

plot_transmission_map <- function(seed = 4,
                                  dgp = "spec_linear",
                                  sex_lab = "M",
                                  scens = c("up1pc", "freeze", "down1pc", "quit"),
                                  raw_dir = OUT_RAW,
                                  title_suffix = NULL,
                                  base_size = paper_fig_base_size("transmission_map")) {
  sex_lab <- match.arg(as.character(sex_lab)[1], c("M", "F"))
  series_levels <- unname(MODEL_LABELS[c("truth", "sbapc")])
  metric_levels <- c(
    "Current smoking",
    "Effective exposure",
    "Incident cases",
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
        `Current smoking` = mean(as.numeric(p_curr), na.rm = TRUE) * 100,
        `Effective exposure` = mean(as.numeric(q_eff), na.rm = TRUE) * 100,
        .groups = "drop"
      ) %>%
      tidyr::pivot_longer(-period, names_to = "metric", values_to = "value") %>%
      dplyr::mutate(series = MODEL_LABELS[["truth"]])

    est_stock <- dplyr::bind_rows(
      tibble::as_tibble(sex_res$diag$z_prev_hist %||% tibble::tibble()),
      tibble::as_tibble(sex_res$diag$z_prev_future %||% tibble::tibble())
    ) %>%
      dplyr::filter(as.character(sex) == sex_lab) %>%
      dplyr::group_by(period) %>%
      dplyr::summarise(
        `Current smoking` = mean(as.numeric(p_cur), na.rm = TRUE) * 100,
        `Effective exposure` = mean(as.numeric(q_eff), na.rm = TRUE) * 100,
        .groups = "drop"
      ) %>%
      tidyr::pivot_longer(-period, names_to = "metric", values_to = "value") %>%
      dplyr::mutate(series = MODEL_LABELS[["sbapc"]])

    exposure <- suppressWarnings(as.numeric(rb$meta$args$exposure %||% sim_truth$meta$exposure %||% 100000))[1]
    if (!is.finite(exposure) || exposure <= 0) exposure <- 100000

    truth_inc <- rb$inc_truth_grid %>%
      dplyr::filter(as.character(sex) == sex_lab) %>%
      dplyr::group_by(period) %>%
      dplyr::summarise(value = sum(as.numeric(rateI_scen_true) * exposure, na.rm = TRUE), .groups = "drop") %>%
      dplyr::mutate(metric = "Incident cases", series = MODEL_LABELS[["truth"]])

    est_inc <- sex_res$inc_annual_cond %>%
      dplyr::transmute(period = as.integer(period), value = as.numeric(cases_hat),
                       metric = "Incident cases", series = MODEL_LABELS[["sbapc"]])
    if (!any(est_inc$period <= 2022, na.rm = TRUE) && !identical(scen, "freeze")) {
      if (is.null(freeze_hist_inc_by_seed)) {
        freeze_file <- file.path(raw_dir, sprintf("res_%s_s%d_freeze.rds", dgp, seed))
        rb_freeze <- read_rds_safe(freeze_file)
        if (!inherits(rb_freeze, "try-error")) {
          sex_freeze <- if (identical(sex_lab, "M")) rb_freeze$resM else rb_freeze$resF
          freeze_hist_inc_by_seed <- sex_freeze$inc_annual_cond %>%
            dplyr::filter(as.integer(period) <= 2022) %>%
            dplyr::transmute(period = as.integer(period), value = as.numeric(cases_hat),
                             metric = "Incident cases", series = MODEL_LABELS[["sbapc"]])
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
      dplyr::mutate(metric = "Annual deaths", series = MODEL_LABELS[["truth"]])

    est_mort <- sex_res$annual_anchor %>%
      dplyr::transmute(period = as.integer(period), value = as.numeric(deaths_hat),
                       metric = "Annual deaths", series = MODEL_LABELS[["sbapc"]])

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
    geom_vline(xintercept = 2022, linetype = "dotted", color = "gray45", linewidth = 0.55) +
    geom_line(linewidth = 0.85, na.rm = TRUE) +
    scale_color_manual(values = MODEL_COLORS[series_levels], breaks = series_levels) +
    scale_linetype_manual(values = MODEL_LINETYPES[series_levels], breaks = series_levels) +
    labs(x = "Year", y = NULL, color = "Series", linetype = "Series") +
    theme_paper_main(base_size = base_size) +
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
                                                  realistic_label = "Observed-window SBAPC",
                                                  oracle_label = "Full-support SBAPC",
                                                  base_size = paper_fig_base_size("transmission_support")) {
  sex_lab <- match.arg(as.character(sex_lab)[1], c("M", "F"))
  metric_levels <- c(
    "Current smoking",
    "Effective exposure",
    "Incident cases",
    "Annual deaths"
  )
  series_levels <- c("Truth", realistic_label, oracle_label)

  g_real <- plot_transmission_map(
    seed = seed, dgp = dgp, sex_lab = sex_lab, scens = scens,
    raw_dir = realistic_raw_dir,
    base_size = base_size
  )
  g_oracle <- plot_transmission_map(
    seed = seed, dgp = dgp, sex_lab = sex_lab, scens = scens,
    raw_dir = oracle_raw_dir,
    base_size = base_size
  )

  df_real <- tibble::as_tibble(g_real$data)
  df_oracle <- tibble::as_tibble(g_oracle$data)

  truth_df <- df_oracle %>%
    dplyr::filter(as.character(series) == MODEL_LABELS[["truth"]]) %>%
    dplyr::mutate(series = MODEL_LABELS[["truth"]])
  real_df <- df_real %>%
    dplyr::filter(as.character(series) == MODEL_LABELS[["sbapc"]]) %>%
    dplyr::mutate(series = realistic_label)
  oracle_df <- df_oracle %>%
    dplyr::filter(as.character(series) == MODEL_LABELS[["sbapc"]]) %>%
    dplyr::mutate(series = oracle_label)

  df <- dplyr::bind_rows(truth_df, real_df, oracle_df) %>%
    dplyr::mutate(
      metric = factor(as.character(metric), levels = metric_levels),
      series = factor(as.character(series), levels = series_levels),
      scenario = factor(as.character(scenario), levels = SCEN_LABELS[scens])
    )

  ggplot(df, aes(x = period, y = value, color = series, linetype = series, linewidth = series)) +
    facet_grid(metric ~ scenario, scales = "free_y") +
    geom_vline(xintercept = 2022, linetype = "dotted", color = "gray60", linewidth = 0.35) +
    geom_line(na.rm = TRUE) +
    scale_color_manual(
      values = stats::setNames(c("black", "#D32F2F", "#1565C0"), series_levels),
      breaks = series_levels
    ) +
    scale_linetype_manual(
      values = stats::setNames(c("solid", "solid", "solid"), series_levels),
      breaks = series_levels
    ) +
    scale_linewidth_manual(
      values = stats::setNames(c(0.80, 0.62, 0.62), series_levels),
      breaks = series_levels,
      guide = "none"
    ) +
    labs(x = "Year", y = NULL, color = "Series", linetype = "Series") +
    theme_paper_main(base_size = base_size) +
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
                                               force_oracle = FALSE,
                                               sexes = "M") {
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

  for (sx in sexes) {
    g <- plot_transmission_map_support_compare(
      seed = seed,
      dgp = dgp,
      sex_lab = sx,
      scens = scens,
      realistic_raw_dir = realistic_raw_dir,
      oracle_raw_dir = oracle_raw_dir,
      realistic_label = "Observed-window SBAPC",
      oracle_label = "Full-support SBAPC"
    )
    save_paper_plot(
      g,
      file.path(OUT_SEC4, sprintf("fig_transmission_map_support_compare_seed%d_%s", seed, sx)),
      width = paper_fig_spec("transmission_support")$width,
      height = paper_fig_spec("transmission_support")$height,
      bg = "white"
    )
  }

  invisible(TRUE)
}

build_mortality_scenario_effects <- function(data = NULL,
                                             sex_scope = c("total", "by_sex"),
                                             scens = setdiff(CANONICAL_SCENS, "freeze"),
                                             last_hist = 2022L) {
  sex_scope <- match.arg(sex_scope)
  if (is.null(data)) data <- extract_all_metrics()
  if (is.null(data$mort) || !nrow(data$mort)) {
    stop("No mortality diagnostics available. Run/extract simulation metrics first.")
  }

  mort <- tibble::as_tibble(data$mort) %>%
    dplyr::filter(period > last_hist) %>%
    dplyr::select(seed, dgp, scenario, sex, period,
                  deaths_true, deaths_hat, deaths_noP, deaths_bapc)

  support <- tibble::as_tibble(data$inc %||% tibble::tibble()) %>%
    dplyr::filter(period > last_hist) %>%
    dplyr::select(seed, dgp, scenario, sex, period, support_frac)

  if (identical(sex_scope, "total")) {
    mort <- mort %>%
      dplyr::group_by(seed, dgp, scenario, period) %>%
      dplyr::summarise(
        deaths_true = sum(deaths_true, na.rm = TRUE),
        deaths_hat = sum(deaths_hat, na.rm = TRUE),
        deaths_noP = sum(deaths_noP, na.rm = TRUE),
        deaths_bapc = sum(deaths_bapc, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::mutate(sex = "Total", .before = period)

    if (nrow(support)) {
      support <- support %>%
        dplyr::group_by(seed, dgp, scenario, period) %>%
        dplyr::summarise(support_frac = mean(as.numeric(support_frac), na.rm = TRUE), .groups = "drop") %>%
        dplyr::mutate(sex = "Total", .before = period)
    }
  }

  freeze <- mort %>%
    dplyr::filter(scenario == "freeze") %>%
    dplyr::select(seed, dgp, sex, period,
                  deaths_true_freeze = deaths_true,
                  deaths_hat_freeze = deaths_hat,
                  deaths_noP_freeze = deaths_noP,
                  deaths_bapc_freeze = deaths_bapc)

  wide <- mort %>%
    dplyr::filter(scenario %in% scens) %>%
    dplyr::left_join(freeze, by = c("seed", "dgp", "sex", "period")) %>%
    dplyr::left_join(support, by = c("seed", "dgp", "scenario", "sex", "period")) %>%
    dplyr::mutate(
      delta_truth = deaths_true - deaths_true_freeze,
      delta_sbapc = deaths_hat - deaths_hat_freeze,
      delta_incidence_anchored = deaths_noP - deaths_noP_freeze,
      delta_bapc = deaths_bapc - deaths_bapc_freeze,
      truth_effect_pct = 100 * delta_truth / pmax(abs(deaths_true_freeze), 1e-9),
      horizon = as.integer(period - last_hist),
      horizon_region = dplyr::case_when(
        is.na(support_frac) ~ NA_character_,
        support_frac >= 0.50 ~ "Credible",
        support_frac >= 0.33 ~ "Caution",
        TRUE ~ "Risky"
      ),
      scenario_label = unname(SCEN_LABELS[as.character(scenario)])
    )

  out <- dplyr::bind_rows(
    wide %>%
      dplyr::transmute(seed, dgp, scenario, scenario_label, sex, period, horizon, horizon_region,
                       model = MODEL_LABELS[["sbapc"]],
                       delta_truth, delta_hat = delta_sbapc,
                       freeze_truth = deaths_true_freeze, support_frac),
    wide %>%
      dplyr::transmute(seed, dgp, scenario, scenario_label, sex, period, horizon, horizon_region,
                       model = MODEL_LABELS[["sbapc_no_prev"]],
                       delta_truth, delta_hat = delta_incidence_anchored,
                       freeze_truth = deaths_true_freeze, support_frac),
    wide %>%
      dplyr::transmute(seed, dgp, scenario, scenario_label, sex, period, horizon, horizon_region,
                       model = MODEL_LABELS[["bapc"]],
                       delta_truth, delta_hat = delta_bapc,
                       freeze_truth = deaths_true_freeze, support_frac)
  ) %>%
    dplyr::mutate(
      scenario = factor(as.character(scenario), levels = scens),
      scenario_label = factor(as.character(scenario_label), levels = unname(SCEN_LABELS[scens])),
      model = factor(as.character(model), levels = unname(MODEL_LABELS[c("sbapc", "sbapc_no_prev", "bapc")])),
      horizon_region = factor(horizon_region, levels = c("Credible", "Caution", "Risky")),
      delta_error = delta_hat - delta_truth,
      effect_hat_pct = 100 * delta_hat / pmax(abs(freeze_truth), 1e-9),
      effect_true_pct = 100 * delta_truth / pmax(abs(freeze_truth), 1e-9)
    )

  out
}

summarise_scenario_effect_recovery <- function(effect_df) {
  seed_level <- effect_df %>%
    dplyr::group_by(model, scenario, scenario_label, seed, dgp, sex) %>%
    dplyr::summarise(
      annual_mare_pct = 100 * sum(abs(delta_error), na.rm = TRUE) / pmax(sum(abs(delta_truth), na.rm = TRUE), 1e-9),
      signed_error_pct = 100 * sum(delta_error, na.rm = TRUE) / pmax(sum(abs(delta_truth), na.rm = TRUE), 1e-9),
      cumulative_recovery_pct = 100 * sum(delta_hat, na.rm = TRUE) / pmax(abs(sum(delta_truth, na.rm = TRUE)), 1e-9) *
        sign(sum(delta_truth, na.rm = TRUE)),
      sign_agreement_pct = {
        keep <- is.finite(delta_truth) & abs(delta_truth) > 1e-6 & is.finite(delta_hat)
        if (any(keep)) mean(sign(delta_hat[keep]) == sign(delta_truth[keep])) * 100 else NA_real_
      },
      .groups = "drop"
    )

  seed_level %>%
    dplyr::group_by(model, scenario, scenario_label, sex) %>%
    dplyr::summarise(
      seeds = dplyr::n_distinct(seed),
      annual_mare_pct = mean(annual_mare_pct, na.rm = TRUE),
      signed_error_pct = mean(signed_error_pct, na.rm = TRUE),
      cumulative_recovery_pct = stats::median(cumulative_recovery_pct, na.rm = TRUE),
      cumulative_recovery_p10 = as.numeric(stats::quantile(cumulative_recovery_pct, 0.10, na.rm = TRUE)),
      cumulative_recovery_p90 = as.numeric(stats::quantile(cumulative_recovery_pct, 0.90, na.rm = TRUE)),
      sign_agreement_pct = mean(sign_agreement_pct, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(scenario, model, sex)
}

plot_scenario_effect_recovery <- function(effect_df,
                                          include_models = unname(MODEL_LABELS[c("sbapc", "bapc")]),
                                          title = NULL,
                                          subtitle = NULL,
                                          base_size = paper_fig_base_size("scenario_effect"),
                                          effect_scale = c("counts", "percent"),
                                          show_band = TRUE,
                                          show_horizon_overlay = TRUE) {
  effect_scale <- match.arg(effect_scale)
  value_col <- if (identical(effect_scale, "counts")) "effect_count" else "effect_pct"
  y_lab <- if (identical(effect_scale, "counts")) {
    "Annual mortality effect (deaths relative to freeze)"
  } else {
    "Scenario effect (% of freeze deaths)"
  }

  truth_df <- effect_df %>%
    dplyr::distinct(seed, dgp, scenario, scenario_label, sex, period, delta_truth, effect_true_pct) %>%
    dplyr::mutate(
      model = MODEL_LABELS[["truth"]],
      effect_count = delta_truth,
      effect_pct = effect_true_pct
    )

  model_df <- effect_df %>%
    dplyr::filter(as.character(model) %in% include_models) %>%
    dplyr::mutate(
      effect_count = delta_hat,
      effect_pct = effect_hat_pct
    )

  plot_df <- dplyr::bind_rows(
    truth_df %>% dplyr::select(seed, dgp, scenario, scenario_label, sex, period, model, dplyr::all_of(value_col)),
    model_df %>% dplyr::select(seed, dgp, scenario, scenario_label, sex, period, model, dplyr::all_of(value_col))
  ) %>%
    dplyr::rename(effect_value = dplyr::all_of(value_col)) %>%
    dplyr::mutate(
      model = factor(as.character(model), levels = c(MODEL_LABELS[["truth"]], include_models)),
      sex = factor(sex_public_label(sex), levels = c("Total", "Male", "Female")),
      scenario_label = factor(as.character(scenario_label), levels = unname(SCEN_LABELS[setdiff(CANONICAL_SCENS, "freeze")]))
    )

  sum_df <- plot_df %>%
    dplyr::group_by(scenario_label, sex, period, model) %>%
    dplyr::summarise(
      p10 = as.numeric(stats::quantile(effect_value, 0.10, na.rm = TRUE)),
      med = stats::median(effect_value, na.rm = TRUE),
      p90 = as.numeric(stats::quantile(effect_value, 0.90, na.rm = TRUE)),
      .groups = "drop"
    )

  ribbon_df <- sum_df %>%
    dplyr::filter(as.character(model) %in% c(MODEL_LABELS[["truth"]], MODEL_LABELS[["sbapc"]]))

  boundaries <- effect_df %>%
    dplyr::group_by(seed, dgp, scenario, sex) %>%
    dplyr::summarise(
      caution_start = suppressWarnings(min(period[support_frac < 0.50], na.rm = TRUE)),
      risky_start = suppressWarnings(min(period[support_frac < 0.33], na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      caution_start = dplyr::if_else(is.finite(caution_start), as.integer(caution_start), NA_integer_),
      risky_start = dplyr::if_else(is.finite(risky_start), as.integer(risky_start), NA_integer_)
    )
  identical_boundaries <- nrow(boundaries) > 0 &&
    dplyr::n_distinct(boundaries$caution_start, na.rm = TRUE) <= 1 &&
    dplyr::n_distinct(boundaries$risky_start, na.rm = TRUE) <= 1
  support_lines <- tibble::tibble()
  support_rects <- tibble::tibble()
  if (isTRUE(show_horizon_overlay) && identical_boundaries) {
    x_min <- min(sum_df$period, na.rm = TRUE)
    x_max <- max(sum_df$period, na.rm = TRUE)
    caution_start <- unique(stats::na.omit(boundaries$caution_start))[1]
    risky_start <- unique(stats::na.omit(boundaries$risky_start))[1]
    support_lines <- tibble::tibble(period = c(caution_start, risky_start)) %>%
      dplyr::filter(is.finite(period))
    support_rects <- tibble::tibble(
      xmin = c(x_min, caution_start, risky_start),
      xmax = c(caution_start, risky_start, x_max),
      fill = c("#E8F3EA", "#FFF8D9", "#FBE4E6")
    ) %>%
      dplyr::filter(is.finite(xmin), is.finite(xmax), xmax > xmin)
  }

  color_values <- MODEL_COLORS[levels(plot_df$model)]
  fill_values <- c(
    stats::setNames("#BDBDBD", MODEL_LABELS[["truth"]]),
    stats::setNames("#D32F2F", MODEL_LABELS[["sbapc"]])
  )
  facet_layer <- if (dplyr::n_distinct(plot_df$sex) <= 1) {
    ggplot2::facet_wrap(~scenario_label, nrow = 1, scales = "free_y")
  } else {
    ggplot2::facet_grid(sex ~ scenario_label, scales = "free_y")
  }

  g <- ggplot2::ggplot(sum_df, ggplot2::aes(x = period, y = med, color = model, linetype = model))
  if (nrow(support_rects)) {
    for (i in seq_len(nrow(support_rects))) {
      g <- g + ggplot2::geom_rect(
        data = support_rects[i, ],
        ggplot2::aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
        inherit.aes = FALSE,
        fill = support_rects$fill[i],
        alpha = 0.45,
        color = NA
      )
    }
  }
  if (isTRUE(show_band) && nrow(ribbon_df)) {
    g <- g +
    ggplot2::geom_ribbon(
      data = ribbon_df,
      ggplot2::aes(ymin = p10, ymax = p90, fill = model),
      alpha = 0.12,
      color = NA,
      show.legend = FALSE
    )
  }
  g +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.35, color = "gray60") +
    ggplot2::geom_vline(
      data = support_lines,
      ggplot2::aes(xintercept = period),
      inherit.aes = FALSE,
      linetype = "dotted",
      color = "gray45",
      linewidth = 0.35
    ) +
    ggplot2::geom_line(linewidth = 0.9, na.rm = TRUE) +
    facet_layer +
    ggplot2::scale_color_manual(values = color_values, breaks = levels(plot_df$model)) +
    ggplot2::scale_linetype_manual(values = MODEL_LINETYPES[levels(plot_df$model)], breaks = levels(plot_df$model)) +
    ggplot2::scale_fill_manual(values = fill_values) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "Year",
      y = y_lab,
      color = "Series",
      linetype = "Series"
    ) +
    theme_paper_main(base_size = base_size) +
    ggplot2::theme(
      legend.position = "bottom",
      strip.background = ggplot2::element_rect(fill = "gray95"),
      panel.spacing = grid::unit(0.7, "lines")
    )
}

export_scenario_effect_recovery_table <- function(summary_df,
                                                  file_out,
                                                  models = unname(MODEL_LABELS[c("sbapc", "bapc")]),
                                                  sex = "Total") {
  csv_out <- sub("\\.tex$", ".csv", file_out)
  tab <- summary_df %>%
    dplyr::filter(as.character(model) %in% models, as.character(sex) == !!sex) %>%
    dplyr::mutate(
      scenario_tex = compact_scenario_label(scenario),
      model = as.character(model),
      recovery = cumulative_recovery_pct / 100
    ) %>%
    dplyr::arrange(factor(as.character(scenario), levels = c("up1pc", "down1pc", "quit")),
                   factor(model, levels = models))
  readr::write_csv(tab, csv_out)

  lines <- c(latex_table_open("llrrr"),
             "Scenario & Series & MARE (\\%) & Recovery & Sign (\\%) \\\\",
             "\\midrule")
  last_scenario <- NULL
  for (i in seq_len(nrow(tab))) {
    row <- tab[i, ]
    scen <- if (!identical(last_scenario, as.character(row$scenario))) row$scenario_tex else ""
    if (!is.null(last_scenario) && !identical(last_scenario, as.character(row$scenario))) lines <- c(lines, "\\midrule")
    lines <- c(lines, sprintf(
      "%s & %s & %.1f & %.2f & %.0f \\\\",
      scen, row$model, row$annual_mare_pct, row$recovery, row$sign_agreement_pct
    ))
    last_scenario <- as.character(row$scenario)
  }
  lines <- c(lines, latex_table_close(
    "Annual MARE and cumulative recovery summarize mortality scenario effects relative to the frozen-prevalence baseline. The BAPC benchmark is scenario-blind and has zero scenario response by construction."
  ))
  writeLines(lines, file_out)
  invisible(tab)
}

export_bias_summary_table <- function(metrics_df,
                                      csv_out = file.path(OUT_SEC4, "tab_bias_summary.csv"),
                                      tex_out = file.path(OUT_SEC4, "tab_bias_summary.tex")) {
  tab <- metrics_df %>%
    dplyr::filter(as.character(dgp) == "spec_linear", as.character(sex) %in% c("M", "F")) %>%
    dplyr::mutate(
      scenario_label = compact_scenario_label(scenario),
      sex_label = sex_public_label(sex)
    ) %>%
    dplyr::group_by(scenario, scenario_label, sex_label) %>%
    dplyr::summarise(
      hist_bias_pct = mean(hist_bias, na.rm = TRUE),
      proj_bias_pct = mean(proj_bias, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(factor(as.character(scenario), levels = CANONICAL_SCENS), sex_label)
  readr::write_csv(tab, csv_out)

  wide <- tab %>%
    dplyr::select(scenario, scenario_label, sex_label, hist_bias_pct, proj_bias_pct) %>%
    tidyr::pivot_wider(
      names_from = sex_label,
      values_from = c(hist_bias_pct, proj_bias_pct)
    ) %>%
    dplyr::arrange(factor(as.character(scenario), levels = CANONICAL_SCENS))

  lines <- c(latex_table_open("lrrrr"),
             "Scenario & Male hist. & Male proj. & Female hist. & Female proj. \\\\",
             "\\midrule")
  for (i in seq_len(nrow(wide))) {
    row <- wide[i, ]
    lines <- c(lines, sprintf(
      "%s & %.1f & %.1f & %.1f & %.1f \\\\",
      row$scenario_label,
      row$hist_bias_pct_Male, row$proj_bias_pct_Male,
      row$hist_bias_pct_Female, row$proj_bias_pct_Female
    ))
  }
  lines <- c(lines, latex_table_close(
    "Entries are mean percentage bias across simulation seeds for the well-specified design. Historical and projected periods are summarized separately by sex."
  ))
  writeLines(lines, tex_out, useBytes = TRUE)
  invisible(tab)
}

latex_scenario_labels <- function() {
  c(
    "up1pc" = "$\\uparrow$ 1\\% p.a.",
    "freeze" = "Freeze (2022)",
    "down1pc" = "$\\downarrow$ 1\\% p.a.",
    "quit" = "Quit"
  )
}

fmt_int <- function(x) {
  ifelse(is.finite(x), formatC(round(x), format = "f", digits = 0, big.mark = ","), "")
}

fmt_num <- function(x, digits = 2) {
  ifelse(is.finite(x), formatC(x, format = "f", digits = digits), "")
}

fmt_interval <- function(mid, lo, hi, digits = 2) {
  sprintf("%s [%s, %s]", fmt_num(mid, digits), fmt_num(lo, digits), fmt_num(hi, digits))
}

fmt_int_interval <- function(mid, lo, hi) {
  sprintf("%s [%s, %s]", fmt_int(mid), fmt_int(lo), fmt_int(hi))
}

latex_note_block <- function(note, width = "0.98\\textwidth") {
  c(
    "\\par\\smallskip",
    sprintf("\\begin{minipage}{%s}", width),
    paste0("\\footnotesize\\emph{Note:} ", note, " \\emph{Source:} Own elaboration."),
    "\\end{minipage}",
    "\\endgroup"
  )
}

latex_table_open <- function(colspec) {
  c(
    "\\begingroup",
    "\\renewcommand{\\arraystretch}{1.08}",
    "\\setlength{\\tabcolsep}{4pt}",
    sprintf("\\begin{tabular}{%s}", colspec),
    "\\toprule"
  )
}

latex_table_close <- function(note = NULL) {
  c("\\bottomrule", "\\end{tabular}", "\\endgroup")
}

compact_scenario_label <- function(x) {
  labs <- c(
    "up1pc" = "$\\uparrow$1\\%",
    "freeze" = "Freeze",
    "down1pc" = "$\\downarrow$1\\%",
    "quit" = "Quit"
  )
  out <- unname(labs[as.character(x)])
  out[is.na(out)] <- as.character(x)[is.na(out)]
  out
}

horizon_boundary_audit <- function(data = NULL,
                                   file_out = file.path(OUT_SEC4, "horizon_boundary_audit.csv"),
                                   outcomes = c("Incidence", "Mortality")) {
  if (is.null(data)) data <- extract_all_metrics()
  audit_base <- tibble::as_tibble(data$inc %||% tibble::tibble()) %>%
    dplyr::filter(period > 2022) %>%
    dplyr::group_by(seed, dgp, scenario, sex) %>%
    dplyr::summarise(
      first_caution_year = suppressWarnings(min(period[support_frac < 0.50], na.rm = TRUE)),
      first_risky_year = suppressWarnings(min(period[support_frac < 0.33], na.rm = TRUE)),
      min_support_frac = min(support_frac, na.rm = TRUE),
      mean_support_frac = mean(support_frac, na.rm = TRUE),
      max_projection_year = max(period, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      first_caution_year = dplyr::if_else(is.finite(first_caution_year), as.integer(first_caution_year), NA_integer_),
      first_risky_year = dplyr::if_else(is.finite(first_risky_year), as.integer(first_risky_year), NA_integer_),
      sex = sex_public_label(sex),
      design = dgp_public_label(dgp),
      support_source = "Incidence support fraction"
    )

  audit <- tidyr::crossing(outcome = outcomes, audit_base) %>%
    dplyr::select(design, dgp, seed, scenario, sex, outcome, support_source,
                  first_caution_year, first_risky_year, min_support_frac,
                  mean_support_frac, max_projection_year)
  readr::write_csv(audit, file_out)
  audit
}

add_horizon_region_rows <- function(df,
                                    period_col = "period",
                                    support_col = "support_frac",
                                    include_full = FALSE) {
  out <- df %>%
    dplyr::mutate(
      horizon_region = dplyr::case_when(
        is.na(.data[[support_col]]) ~ NA_character_,
        .data[[support_col]] >= 0.50 ~ "Credible",
        .data[[support_col]] >= 0.33 ~ "Caution",
        TRUE ~ "Risky"
      ),
      horizon_region = factor(horizon_region, levels = c("Credible", "Caution", "Risky"))
    )
  if (isTRUE(include_full)) {
    out <- dplyr::bind_rows(
      out,
      out %>% dplyr::mutate(horizon_region = factor("Full horizon", levels = c("Credible", "Caution", "Risky", "Full horizon")))
    )
  }
  out %>%
    dplyr::mutate(
      horizon_region = factor(as.character(horizon_region),
                              levels = c("Credible", "Caution", "Risky", "Full horizon"))
    )
}

summarise_cumulative_scenario_recovery <- function(effect_df,
                                                   model_label = MODEL_LABELS[["sbapc"]],
                                                   sex = "Total",
                                                   include_full = TRUE) {
  region_df <- effect_df %>%
    dplyr::filter(as.character(model) == model_label, as.character(sex) == !!sex) %>%
    add_horizon_region_rows(include_full = include_full) %>%
    dplyr::filter(!is.na(horizon_region))

  seed_level <- region_df %>%
    dplyr::group_by(scenario, scenario_label, horizon_region, seed, dgp, sex) %>%
    dplyr::summarise(
      true_cumulative = sum(delta_truth, na.rm = TRUE),
      estimated_cumulative = sum(delta_hat, na.rm = TRUE),
      recovery_ratio = estimated_cumulative / dplyr::if_else(abs(true_cumulative) > 1e-9, true_cumulative, NA_real_),
      annual_mare_pct = 100 * sum(abs(delta_hat - delta_truth), na.rm = TRUE) /
        pmax(sum(abs(delta_truth), na.rm = TRUE), 1e-9),
      sign_agreement_pct = {
        keep <- is.finite(delta_truth) & abs(delta_truth) > 1e-6 & is.finite(delta_hat)
        if (any(keep)) mean(sign(delta_hat[keep]) == sign(delta_truth[keep])) * 100 else NA_real_
      },
      .groups = "drop"
    )

  seed_level %>%
    dplyr::group_by(scenario, scenario_label, horizon_region, sex) %>%
    dplyr::summarise(
      seeds = dplyr::n_distinct(seed),
      true_cumulative_median = stats::median(true_cumulative, na.rm = TRUE),
      true_cumulative_q25 = as.numeric(stats::quantile(true_cumulative, 0.25, na.rm = TRUE)),
      true_cumulative_q75 = as.numeric(stats::quantile(true_cumulative, 0.75, na.rm = TRUE)),
      estimated_cumulative_median = stats::median(estimated_cumulative, na.rm = TRUE),
      estimated_cumulative_q25 = as.numeric(stats::quantile(estimated_cumulative, 0.25, na.rm = TRUE)),
      estimated_cumulative_q75 = as.numeric(stats::quantile(estimated_cumulative, 0.75, na.rm = TRUE)),
      recovery_ratio_median = stats::median(recovery_ratio, na.rm = TRUE),
      recovery_ratio_q25 = as.numeric(stats::quantile(recovery_ratio, 0.25, na.rm = TRUE)),
      recovery_ratio_q75 = as.numeric(stats::quantile(recovery_ratio, 0.75, na.rm = TRUE)),
      annual_mare_pct_mean = mean(annual_mare_pct, na.rm = TRUE),
      annual_mare_pct_sd = stats::sd(annual_mare_pct, na.rm = TRUE),
      sign_agreement_pct_mean = mean(sign_agreement_pct, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(factor(as.character(scenario), levels = c("up1pc", "down1pc", "quit")),
                   factor(as.character(horizon_region), levels = c("Credible", "Caution", "Risky", "Full horizon")))
}

export_cumulative_scenario_recovery_table <- function(summary_df,
                                                      csv_out = file.path(OUT_SEC4, "tab_cumulative_scenario_recovery.csv"),
                                                      tex_out = file.path(OUT_SEC4, "tab_cumulative_scenario_recovery.tex")) {
  table_df <- summary_df %>%
    dplyr::filter(as.character(horizon_region) %in% c("Credible", "Caution", "Risky"))
  readr::write_csv(table_df, csv_out)
  tab <- table_df %>%
    dplyr::mutate(
      scenario_tex = compact_scenario_label(scenario),
      horizon = as.character(horizon_region),
      true_display = fmt_int(true_cumulative_median),
      est_display = fmt_int(estimated_cumulative_median),
      recovery_display = fmt_num(recovery_ratio_median, digits = 2)
    )

  lines <- c(latex_table_open("llrrrr"),
             "Scenario & Horizon & Truth & SBAPC & Recovery & MARE (\\%) \\\\",
             "\\midrule")
  last_scenario <- NULL
  for (i in seq_len(nrow(tab))) {
    row <- tab[i, ]
    scen <- if (!identical(last_scenario, as.character(row$scenario))) row$scenario_tex else ""
    if (!is.null(last_scenario) && !identical(last_scenario, as.character(row$scenario))) lines <- c(lines, "\\midrule")
    lines <- c(lines, sprintf(
      "%s & %s & %s & %s & %s & %.1f \\\\",
      scen, row$horizon, row$true_display, row$est_display,
      row$recovery_display, row$annual_mare_pct_mean
    ))
    last_scenario <- as.character(row$scenario)
  }
  lines <- c(lines, latex_table_close(
    "Entries are medians across simulation seeds for cumulative annual mortality effects relative to the frozen-prevalence baseline; MARE is averaged across seeds within each horizon region."
  ))
  writeLines(lines, tex_out, useBytes = TRUE)
  invisible(tab)
}

build_chain_recovery_data <- function(seeds = CANONICAL_SEEDS,
                                      dgp = "spec_linear",
                                      raw_dir = OUT_RAW,
                                      last_hist = 2022L) {
  rows <- list()
  for (seed in seeds) {
    rds_file <- file.path(raw_dir, sprintf("res_%s_s%d_freeze.rds", dgp, seed))
    if (!file.exists(rds_file)) next
    rb <- read_rds_safe(rds_file)
    if (inherits(rb, "try-error") || is.null(rb$inc_truth_grid) || is.null(rb$mort_truth_grid)) next

    sim <- simulate_PIM_data(cause_id = CAUSE_ID, seed = seed, dgp = dgp, scenario_name = "freeze")
    pop <- tibble::as_tibble(sim$pop_all) %>%
      dplyr::mutate(sex = as.character(sex), age = as.integer(age), period = as.integer(period))

    truth_eff <- tibble::as_tibble(sim$z_scen_true) %>%
      dplyr::mutate(sex = as.character(sex), age = as.integer(age), period = as.integer(period)) %>%
      dplyr::left_join(pop, by = c("sex", "age", "period")) %>%
      dplyr::filter(period > last_hist) %>%
      dplyr::group_by(period) %>%
      dplyr::summarise(truth = stats::weighted.mean(as.numeric(q_eff), exposure, na.rm = TRUE) * 100,
                       .groups = "drop")

    est_eff <- dplyr::bind_rows(
      tibble::as_tibble(rb$resM$diag$z_prev_hist %||% tibble::tibble()),
      tibble::as_tibble(rb$resM$diag$z_prev_future %||% tibble::tibble()),
      tibble::as_tibble(rb$resF$diag$z_prev_hist %||% tibble::tibble()),
      tibble::as_tibble(rb$resF$diag$z_prev_future %||% tibble::tibble())
    ) %>%
      dplyr::mutate(sex = as.character(sex), age = as.integer(age), period = as.integer(period)) %>%
      dplyr::left_join(pop, by = c("sex", "age", "period")) %>%
      dplyr::filter(period > last_hist) %>%
      dplyr::group_by(period) %>%
      dplyr::summarise(estimate = stats::weighted.mean(as.numeric(q_eff), exposure, na.rm = TRUE) * 100,
                       support_frac = mean(as.numeric(support_frac), na.rm = TRUE),
                       .groups = "drop")

    effective_exposure <- truth_eff %>%
      dplyr::left_join(est_eff, by = "period") %>%
      dplyr::mutate(object = "Effective exposure")

    incidence <- dplyr::bind_rows(
      rb$resM$inc_annual_cond %>% dplyr::mutate(sex = "M"),
      rb$resF$inc_annual_cond %>% dplyr::mutate(sex = "F")
    ) %>%
      dplyr::filter(period > last_hist) %>%
      dplyr::group_by(period) %>%
      dplyr::summarise(estimate = sum(as.numeric(cases_hat), na.rm = TRUE), .groups = "drop") %>%
      dplyr::left_join(
        tibble::as_tibble(rb$inc_truth_grid) %>%
          dplyr::mutate(sex = as.character(sex), age = as.integer(age), period = as.integer(period)) %>%
          dplyr::left_join(pop %>% dplyr::select(sex, age, period, exposure),
                           by = c("sex", "age", "period")) %>%
          dplyr::filter(period > last_hist) %>%
          dplyr::group_by(period) %>%
          dplyr::summarise(truth = sum(as.numeric(rateI_scen_true) * as.numeric(exposure), na.rm = TRUE),
                           .groups = "drop"),
        by = "period"
      ) %>%
      dplyr::left_join(est_eff %>% dplyr::select(period, support_frac), by = "period") %>%
      dplyr::mutate(object = "Incident cases")

    expected_deaths <- dplyr::bind_rows(
      rb$resM$annual_external_cond,
      rb$resF$annual_external_cond
    ) %>%
      dplyr::filter(period > last_hist) %>%
      dplyr::group_by(period) %>%
      dplyr::summarise(estimate = sum(as.numeric(deaths_ext), na.rm = TRUE), .groups = "drop") %>%
      dplyr::left_join(
        tibble::as_tibble(rb$mort_truth_grid) %>%
          dplyr::filter(period > last_hist) %>%
          dplyr::group_by(period) %>%
          dplyr::summarise(truth = sum(as.numeric(mort_deaths_scen_true), na.rm = TRUE), .groups = "drop"),
        by = "period"
      ) %>%
      dplyr::left_join(est_eff %>% dplyr::select(period, support_frac), by = "period") %>%
      dplyr::mutate(object = "Incidence-linked expected deaths")

    mortality <- dplyr::bind_rows(rb$resM$annual_anchor, rb$resF$annual_anchor) %>%
      dplyr::filter(period > last_hist) %>%
      dplyr::group_by(period) %>%
      dplyr::summarise(estimate = sum(as.numeric(deaths_hat), na.rm = TRUE), .groups = "drop") %>%
      dplyr::left_join(
        tibble::as_tibble(rb$mort_truth_grid) %>%
          dplyr::filter(period > last_hist) %>%
          dplyr::group_by(period) %>%
          dplyr::summarise(truth = sum(as.numeric(mort_deaths_scen_true), na.rm = TRUE), .groups = "drop"),
        by = "period"
      ) %>%
      dplyr::left_join(est_eff %>% dplyr::select(period, support_frac), by = "period") %>%
      dplyr::mutate(object = "Mortality deaths")

    rows[[as.character(seed)]] <- dplyr::bind_rows(effective_exposure, incidence, expected_deaths, mortality) %>%
      dplyr::mutate(seed = seed, dgp = dgp, scenario = "freeze", .before = 1)
  }

  dplyr::bind_rows(rows) %>%
    add_horizon_region_rows(include_full = TRUE) %>%
    dplyr::mutate(
      design = dgp_public_label(dgp),
      object = factor(object, levels = c("Effective exposure", "Incident cases",
                                         "Incidence-linked expected deaths", "Mortality deaths"))
    )
}

summarise_chain_recovery <- function(chain_df) {
  seed_level <- chain_df %>%
    dplyr::filter(!is.na(horizon_region)) %>%
    dplyr::group_by(object, horizon_region, seed, dgp) %>%
    dplyr::summarise(
      true_cumulative = sum(truth, na.rm = TRUE),
      estimated_cumulative = sum(estimate, na.rm = TRUE),
      recovery_ratio = estimated_cumulative / dplyr::if_else(abs(true_cumulative) > 1e-9, true_cumulative, NA_real_),
      annual_mare_pct = 100 * sum(abs(estimate - truth), na.rm = TRUE) / pmax(sum(abs(truth), na.rm = TRUE), 1e-9),
      annual_bias_pct = 100 * sum(estimate - truth, na.rm = TRUE) / pmax(sum(abs(truth), na.rm = TRUE), 1e-9),
      .groups = "drop"
    )

  seed_level %>%
    dplyr::group_by(object, horizon_region) %>%
    dplyr::summarise(
      seeds = dplyr::n_distinct(seed),
      recovery_ratio_median = stats::median(recovery_ratio, na.rm = TRUE),
      recovery_ratio_q25 = as.numeric(stats::quantile(recovery_ratio, 0.25, na.rm = TRUE)),
      recovery_ratio_q75 = as.numeric(stats::quantile(recovery_ratio, 0.75, na.rm = TRUE)),
      annual_mare_pct_mean = mean(annual_mare_pct, na.rm = TRUE),
      annual_bias_pct_mean = mean(annual_bias_pct, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(object, factor(as.character(horizon_region), levels = c("Credible", "Caution", "Risky", "Full horizon")))
}

export_chain_recovery_table <- function(summary_df,
                                        csv_out = file.path(OUT_SEC4, "tab_chain_recovery.csv"),
                                        tex_out = file.path(OUT_SEC4, "tab_chain_recovery.tex")) {
  out_csv <- summary_df %>%
    dplyr::filter(as.character(horizon_region) %in% c("Credible", "Caution", "Risky"))
  readr::write_csv(out_csv, csv_out)
  tab <- summary_df %>%
    dplyr::filter(as.character(horizon_region) %in% c("Credible", "Caution", "Risky")) %>%
    dplyr::mutate(recovery_display = fmt_interval(recovery_ratio_median, recovery_ratio_q25, recovery_ratio_q75, digits = 2))
  lines <- c(latex_table_open("llrrr"),
             "Object & Horizon & Recovery & MARE (\\%) & Bias (\\%) \\\\",
             "\\midrule")
  last_object <- NULL
  for (i in seq_len(nrow(tab))) {
    row <- tab[i, ]
    obj <- if (!identical(last_object, as.character(row$object))) as.character(row$object) else ""
    if (!is.null(last_object) && !identical(last_object, as.character(row$object))) lines <- c(lines, "\\midrule")
    lines <- c(lines, sprintf(
      "%s & %s & %s & %.1f & %.1f \\\\",
      obj, as.character(row$horizon_region), row$recovery_display,
      row$annual_mare_pct_mean, row$annual_bias_pct_mean
    ))
    last_object <- as.character(row$object)
  }
  lines <- c(lines, latex_table_close(
    "Effective exposure is aggregated as an exposure-weighted annual mean across sex-age cells. Incident cases, incidence-linked expected deaths, and mortality deaths are aggregated as annual counts across sex-age cells. Rows report endogenous horizon regions only."
  ))
  writeLines(lines, tex_out, useBytes = TRUE)
  invisible(tab)
}

available_result_seeds <- function(raw_dir, dgp = "spec_linear", scens = CANONICAL_SCENS) {
  seeds <- integer(0)
  files <- list.files(raw_dir, pattern = paste0("^res_", dgp, "_s[0-9]+_.*\\.rds$"), full.names = FALSE)
  if (!length(files)) return(seeds)
  candidates <- sort(unique(as.integer(sub(paste0("^res_", dgp, "_s([0-9]+)_.*$"), "\\1", files))))
  candidates[vapply(candidates, function(seed) {
    all(vapply(scens, function(sc) file.exists(file.path(raw_dir, sprintf("res_%s_s%d_%s.rds", dgp, seed, sc))), logical(1)))
  }, logical(1))]
}

plot_support_window_comparison <- function(realistic_effect,
                                           oracle_effect,
                                           include_band = TRUE,
                                           base_size = paper_fig_base_size("scenario_effect_bysex")) {
  common_keys <- realistic_effect %>%
    dplyr::distinct(seed, dgp, scenario, sex, period)
  oracle_effect <- oracle_effect %>%
    dplyr::semi_join(common_keys, by = c("seed", "dgp", "scenario", "sex", "period"))
  realistic_effect <- realistic_effect %>%
    dplyr::semi_join(oracle_effect %>% dplyr::distinct(seed, dgp, scenario, sex, period),
                     by = c("seed", "dgp", "scenario", "sex", "period"))

  series_levels <- c("Truth", "Observed-window SBAPC", "Full-support SBAPC")
  plot_df <- dplyr::bind_rows(
    realistic_effect %>%
      dplyr::distinct(seed, dgp, scenario, scenario_label, sex, period, delta_truth) %>%
      dplyr::mutate(series = "Truth", value = delta_truth),
    oracle_effect %>%
      dplyr::filter(as.character(model) == MODEL_LABELS[["sbapc"]]) %>%
      dplyr::mutate(series = "Full-support SBAPC", value = delta_hat),
    realistic_effect %>%
      dplyr::filter(as.character(model) == MODEL_LABELS[["sbapc"]]) %>%
      dplyr::mutate(series = "Observed-window SBAPC", value = delta_hat)
  ) %>%
    dplyr::mutate(
      series = factor(series, levels = series_levels),
      scenario_label = factor(as.character(scenario_label), levels = unname(SCEN_LABELS[setdiff(CANONICAL_SCENS, "freeze")]))
    )

  sum_df <- plot_df %>%
    dplyr::group_by(scenario_label, period, series) %>%
    dplyr::summarise(
      p10 = as.numeric(stats::quantile(value, 0.10, na.rm = TRUE)),
      med = stats::median(value, na.rm = TRUE),
      p90 = as.numeric(stats::quantile(value, 0.90, na.rm = TRUE)),
      .groups = "drop"
    )

  bounds <- realistic_effect %>%
    dplyr::group_by(seed, dgp, scenario, sex) %>%
    dplyr::summarise(
      caution_start = suppressWarnings(min(period[support_frac < 0.50], na.rm = TRUE)),
      risky_start = suppressWarnings(min(period[support_frac < 0.33], na.rm = TRUE)),
      .groups = "drop"
    )
  identical_boundaries <- nrow(bounds) > 0 &&
    dplyr::n_distinct(bounds$caution_start, na.rm = TRUE) <= 1 &&
    dplyr::n_distinct(bounds$risky_start, na.rm = TRUE) <= 1
  rects <- tibble::tibble()
  vlines <- tibble::tibble()
  if (identical_boundaries) {
    x_min <- min(sum_df$period, na.rm = TRUE)
    x_max <- max(sum_df$period, na.rm = TRUE)
    caution_start <- unique(stats::na.omit(bounds$caution_start))[1]
    risky_start <- unique(stats::na.omit(bounds$risky_start))[1]
    rects <- tibble::tibble(
      xmin = c(x_min, caution_start, risky_start),
      xmax = c(caution_start, risky_start, x_max),
      fill = c("#E8F3EA", "#FFF8D9", "#FBE4E6")
    ) %>% dplyr::filter(is.finite(xmin), is.finite(xmax), xmax > xmin)
    vlines <- tibble::tibble(period = c(caution_start, risky_start)) %>% dplyr::filter(is.finite(period))
  }

  pal <- c("Truth" = "black", "Observed-window SBAPC" = "#B71C1C", "Full-support SBAPC" = "#6FA8DC")
  ltys <- c("Truth" = "dashed", "Observed-window SBAPC" = "solid", "Full-support SBAPC" = "solid")
  g <- ggplot2::ggplot(sum_df, ggplot2::aes(x = period, y = med, color = series, linetype = series))
  if (nrow(rects)) {
    for (i in seq_len(nrow(rects))) {
      g <- g + ggplot2::geom_rect(
        data = rects[i, ],
        ggplot2::aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
        inherit.aes = FALSE,
        fill = rects$fill[i],
        alpha = 0.45,
        color = NA
      )
    }
  }
  if (isTRUE(include_band)) {
    g <- g + ggplot2::geom_ribbon(
      data = sum_df %>% dplyr::filter(as.character(series) %in% c("Truth", "Observed-window SBAPC", "Full-support SBAPC")),
      ggplot2::aes(ymin = p10, ymax = p90, fill = series),
      alpha = 0.08,
      color = NA,
      show.legend = FALSE
    )
  }
  g +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.35, color = "gray60") +
    ggplot2::geom_vline(
      data = vlines,
      ggplot2::aes(xintercept = period),
      inherit.aes = FALSE,
      linetype = "solid",
      color = "gray65",
      linewidth = 0.25
    ) +
    ggplot2::geom_line(linewidth = 0.72, na.rm = TRUE) +
    ggplot2::facet_wrap(~scenario_label, nrow = 1, scales = "free_y") +
    ggplot2::scale_color_manual(values = pal, breaks = series_levels) +
    ggplot2::scale_linetype_manual(values = ltys, breaks = series_levels) +
    ggplot2::scale_fill_manual(values = pal, breaks = series_levels) +
    ggplot2::labs(x = "Year", y = "Annual mortality effect (deaths relative to freeze)",
                  color = "Series", linetype = "Series") +
    theme_paper_main(base_size = base_size) +
    ggplot2::theme(
      legend.position = "bottom",
      strip.background = ggplot2::element_rect(fill = "gray95"),
      panel.spacing = grid::unit(0.7, "lines")
    )
}

summarise_support_window_comparison <- function(realistic_effect,
                                                oracle_effect,
                                                include_full = TRUE) {
  real <- realistic_effect %>%
    dplyr::filter(as.character(model) == MODEL_LABELS[["sbapc"]]) %>%
    dplyr::select(seed, dgp, scenario, scenario_label, sex, period, delta_truth,
                  delta_observed = delta_hat, support_frac)
  full <- oracle_effect %>%
    dplyr::filter(as.character(model) == MODEL_LABELS[["sbapc"]]) %>%
    dplyr::select(seed, dgp, scenario, sex, period, delta_full = delta_hat)

  joined <- real %>%
    dplyr::inner_join(full, by = c("seed", "dgp", "scenario", "sex", "period")) %>%
    add_horizon_region_rows(include_full = include_full) %>%
    dplyr::filter(!is.na(horizon_region))

  seed_level <- joined %>%
    dplyr::group_by(scenario, scenario_label, horizon_region, seed, dgp, sex) %>%
    dplyr::summarise(
      true_cumulative = sum(delta_truth, na.rm = TRUE),
      full_cumulative = sum(delta_full, na.rm = TRUE),
      observed_cumulative = sum(delta_observed, na.rm = TRUE),
      full_annual_mare_pct = 100 * sum(abs(delta_full - delta_truth), na.rm = TRUE) / pmax(sum(abs(delta_truth), na.rm = TRUE), 1e-9),
      observed_annual_mare_pct = 100 * sum(abs(delta_observed - delta_truth), na.rm = TRUE) / pmax(sum(abs(delta_truth), na.rm = TRUE), 1e-9),
      full_recovery_ratio = full_cumulative / dplyr::if_else(abs(true_cumulative) > 1e-9, true_cumulative, NA_real_),
      observed_recovery_ratio = observed_cumulative / dplyr::if_else(abs(true_cumulative) > 1e-9, true_cumulative, NA_real_),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      mare_difference_pct = observed_annual_mare_pct - full_annual_mare_pct,
      recovery_difference = observed_recovery_ratio - full_recovery_ratio
    )

  seed_level %>%
    dplyr::group_by(scenario, scenario_label, horizon_region, sex) %>%
    dplyr::summarise(
      seeds = dplyr::n_distinct(seed),
      full_annual_mare_pct = mean(full_annual_mare_pct, na.rm = TRUE),
      observed_annual_mare_pct = mean(observed_annual_mare_pct, na.rm = TRUE),
      mare_difference_pct = mean(mare_difference_pct, na.rm = TRUE),
      full_recovery_ratio = stats::median(full_recovery_ratio, na.rm = TRUE),
      observed_recovery_ratio = stats::median(observed_recovery_ratio, na.rm = TRUE),
      recovery_difference = stats::median(recovery_difference, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(factor(as.character(scenario), levels = c("up1pc", "down1pc", "quit")),
                   factor(as.character(horizon_region), levels = c("Credible", "Caution", "Risky", "Full horizon")),
                   sex)
}

export_support_window_table <- function(summary_df,
                                        csv_out = file.path(OUT_APPENDIX, "tab_support_window_comparison.csv"),
                                        tex_out = file.path(OUT_APPENDIX, "tab_support_window_comparison.tex")) {
  table_df <- summary_df %>%
    dplyr::filter(as.character(horizon_region) %in% c("Credible", "Caution", "Risky"))
  readr::write_csv(table_df, csv_out)
  tab <- table_df %>%
    dplyr::filter(as.character(sex) == "Total") %>%
    dplyr::mutate(scenario_tex = compact_scenario_label(scenario))
  lines <- c(latex_table_open("llrrrr"),
             "Scenario & Horizon & Obs. MARE & Full MARE & Obs. rec. & Full rec. \\\\",
             "\\midrule")
  last_scenario <- NULL
  for (i in seq_len(nrow(tab))) {
    row <- tab[i, ]
    scen <- if (!identical(last_scenario, as.character(row$scenario))) row$scenario_tex else ""
    if (!is.null(last_scenario) && !identical(last_scenario, as.character(row$scenario))) lines <- c(lines, "\\midrule")
    lines <- c(lines, sprintf(
      "%s & %s & %.1f & %.1f & %.2f & %.2f \\\\",
      scen, as.character(row$horizon_region),
      row$observed_annual_mare_pct, row$full_annual_mare_pct,
      row$observed_recovery_ratio, row$full_recovery_ratio
    ))
    last_scenario <- as.character(row$scenario)
  }
  lines <- c(lines, latex_table_close(
    "Observed-window SBAPC is the realistic estimator; Full-support SBAPC is an oracle-style diagnostic using the broader latent support available in the simulation truth. MARE is annual absolute relative error for mortality scenario effects."
  ))
  writeLines(lines, tex_out, useBytes = TRUE)
  invisible(tab)
}

filter_extracted_data <- function(data, seeds = NULL, dgps = NULL, scens = NULL) {
  out <- data
  for (nm in names(out)) {
    if (!is.data.frame(out[[nm]])) next
    df <- out[[nm]]
    if (!is.null(seeds) && "seed" %in% names(df)) {
      df <- df %>% dplyr::filter(seed %in% !!seeds)
    }
    if (!is.null(dgps) && "dgp" %in% names(df)) {
      df <- df %>% dplyr::filter(as.character(dgp) %in% !!as.character(dgps))
    }
    if (!is.null(scens) && "scenario" %in% names(df)) {
      df <- df %>% dplyr::filter(as.character(scenario) %in% !!as.character(scens))
    }
    out[[nm]] <- df
  }
  out
}

bind_extracted_data <- function(...) {
  inputs <- list(...)
  keys <- unique(unlist(lapply(inputs, names), use.names = FALSE))
  stats::setNames(lapply(keys, function(key) {
    vals <- lapply(inputs, function(x) x[[key]])
    vals <- vals[vapply(vals, is.data.frame, logical(1))]
    if (!length(vals)) return(tibble::tibble())
    dplyr::bind_rows(vals)
  }), keys)
}

write_support_window_interpretation_notes <- function(summary_df,
                                                      common_seeds,
                                                      file_out = file.path(OUT_APPENDIX, "support_window_interpretation_notes.md")) {
  if (!nrow(summary_df)) {
    lines <- c(
      "# Support-Window Interpretation Notes",
      "",
      "The support-window comparison was not generated because no common complete seed set was available for the observed-window and full-support estimators."
    )
    writeLines(lines, file_out, useBytes = TRUE)
    return(invisible(file_out))
  }

  total <- summary_df %>% dplyr::filter(as.character(sex) == "Total")
  mean_diff <- mean(total$mare_difference_pct, na.rm = TRUE)
  by_region <- total %>%
    dplyr::mutate(
      horizon_region = factor(as.character(horizon_region),
                              levels = c("Credible", "Caution", "Risky", "Full horizon"))
    ) %>%
    dplyr::group_by(horizon_region) %>%
    dplyr::summarise(mare_diff = mean(mare_difference_pct, na.rm = TRUE), .groups = "drop")
  region_txt <- paste(
    sprintf("%s: %.1f percentage points", as.character(by_region$horizon_region), by_region$mare_diff),
    collapse = "; "
  )
  sex_gap <- summary_df %>%
    dplyr::filter(as.character(sex) %in% c("M", "F", "Male", "Female")) %>%
    dplyr::group_by(sex) %>%
    dplyr::summarise(mare_diff = mean(mare_difference_pct, na.rm = TRUE), .groups = "drop")
  sex_txt <- if (nrow(sex_gap)) {
    paste(sprintf("%s: %.1f percentage points", sex_public_label(sex_gap$sex), sex_gap$mare_diff), collapse = "; ")
  } else {
    "Sex-specific support-window summaries were not included in the compact table."
  }

  lines <- c(
    "# Support-Window Interpretation Notes",
    "",
    sprintf("This diagnostic uses %d common complete simulation seed(s) available in both the observed-window and full-support result directories.", length(common_seeds)),
    "",
    "1. Full-support SBAPC is compared with observed-window SBAPC as an oracle-style diagnostic, not as a feasible empirical estimator.",
    sprintf("2. On average, the full-support estimator changes annual MARE by %.1f percentage points relative to observed-window SBAPC, where positive values mean lower error under full support.", mean_diff),
    sprintf("3. By horizon region, the observed-window MARE minus full-support MARE is: %s.", region_txt),
    "4. The comparison should be interpreted as evidence about support truncation only. It does not isolate every remaining source of discrepancy from Truth, including smoothing, APC extrapolation, and the sequential mortality mapping.",
    sprintf("5. Sex-specific differences: %s", if (grepl("[.!?]$", sex_txt)) sex_txt else paste0(sex_txt, "."))
  )
  writeLines(lines, file_out, useBytes = TRUE)
  invisible(file_out)
}

generate_support_window_products <- function(realistic_data = NULL,
                                             oracle_data = NULL,
                                             dgp = "spec_linear") {
  realistic_seeds <- available_result_seeds(OUT_RAW, dgp = dgp)
  oracle_seeds <- available_result_seeds(OUT_RAW_ORACLE, dgp = dgp)
  common_seeds <- sort(intersect(realistic_seeds, oracle_seeds))

  if (!length(common_seeds)) {
    write_support_window_interpretation_notes(tibble::tibble(), integer(0))
    return(invisible(list(status = "missing_common_seeds", seeds = integer(0))))
  }

  if (is.null(realistic_data)) {
    realistic_data <- extract_all_metrics(seeds = common_seeds, dgps = dgp, raw_dir = OUT_RAW,
                                          cache_suffix = paste0("support_realistic_", dgp))
  }
  if (is.null(oracle_data)) {
    oracle_data <- extract_all_metrics(seeds = common_seeds, dgps = dgp, raw_dir = OUT_RAW_ORACLE,
                                       cache_suffix = paste0("support_oracle_", dgp))
  }
  realistic_data <- filter_extracted_data(realistic_data, seeds = common_seeds, dgps = dgp)
  oracle_data <- filter_extracted_data(oracle_data, seeds = common_seeds, dgps = dgp)

  realistic_effect <- build_mortality_scenario_effects(data = realistic_data, sex_scope = "total")
  oracle_effect <- build_mortality_scenario_effects(data = oracle_data, sex_scope = "total")

  g <- plot_support_window_comparison(
    realistic_effect = realistic_effect,
    oracle_effect = oracle_effect,
    include_band = length(common_seeds) >= 5,
    base_size = paper_fig_base_size("support_window")
  )
  save_profiled_plot(g, file.path(OUT_APPENDIX, "fig_support_window_comparison"),
                     key = "support_window", bg = "white")

  support_summary <- summarise_support_window_comparison(realistic_effect, oracle_effect)
  export_support_window_table(support_summary)
  write_support_window_interpretation_notes(support_summary, common_seeds)

  invisible(list(status = "generated", seeds = common_seeds, summary = support_summary))
}

plot_misspecification_scenario_recovery <- function(effect_df,
                                                    show_band = TRUE,
                                                    show_horizon_overlay = TRUE,
                                                    base_size = paper_fig_base_size("misspecification")) {
  truth_df <- effect_df %>%
    dplyr::distinct(seed, dgp, scenario, scenario_label, sex, period, delta_truth, support_frac) %>%
    dplyr::mutate(series = MODEL_LABELS[["truth"]], value = delta_truth)
  sbapc_df <- effect_df %>%
    dplyr::filter(as.character(model) == MODEL_LABELS[["sbapc"]]) %>%
    dplyr::mutate(series = MODEL_LABELS[["sbapc"]], value = delta_hat)
  plot_df <- dplyr::bind_rows(
    truth_df %>% dplyr::select(seed, dgp, scenario, scenario_label, sex, period, support_frac, series, value),
    sbapc_df %>% dplyr::select(seed, dgp, scenario, scenario_label, sex, period, support_frac, series, value)
  ) %>%
    dplyr::mutate(
      design = factor(dgp_public_label(dgp), levels = unname(DGP_LABELS[c("spec_linear", "misspec_tanh")])),
      scenario_label = factor(as.character(scenario_label), levels = unname(SCEN_LABELS[setdiff(CANONICAL_SCENS, "freeze")])),
      series = factor(series, levels = c(MODEL_LABELS[["truth"]], MODEL_LABELS[["sbapc"]]))
    )

  sum_df <- plot_df %>%
    dplyr::group_by(design, scenario_label, period, series) %>%
    dplyr::summarise(
      p10 = as.numeric(stats::quantile(value, 0.10, na.rm = TRUE)),
      med = stats::median(value, na.rm = TRUE),
      p90 = as.numeric(stats::quantile(value, 0.90, na.rm = TRUE)),
      .groups = "drop"
    )

  bounds <- plot_df %>%
    dplyr::group_by(seed, dgp, scenario, sex) %>%
    dplyr::summarise(
      caution_start = suppressWarnings(min(period[support_frac < 0.50], na.rm = TRUE)),
      risky_start = suppressWarnings(min(period[support_frac < 0.33], na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      caution_start = dplyr::if_else(is.finite(caution_start), as.integer(caution_start), NA_integer_),
      risky_start = dplyr::if_else(is.finite(risky_start), as.integer(risky_start), NA_integer_)
    )
  identical_boundaries <- nrow(bounds) > 0 &&
    dplyr::n_distinct(bounds$caution_start, na.rm = TRUE) <= 1 &&
    dplyr::n_distinct(bounds$risky_start, na.rm = TRUE) <= 1
  rects <- tibble::tibble()
  vlines <- tibble::tibble()
  if (isTRUE(show_horizon_overlay) && identical_boundaries) {
    x_min <- min(sum_df$period, na.rm = TRUE)
    x_max <- max(sum_df$period, na.rm = TRUE)
    caution_start <- unique(stats::na.omit(bounds$caution_start))[1]
    risky_start <- unique(stats::na.omit(bounds$risky_start))[1]
    rects <- tibble::tibble(
      xmin = c(x_min, caution_start, risky_start),
      xmax = c(caution_start, risky_start, x_max),
      fill = c("#E8F3EA", "#FFF8D9", "#FBE4E6")
    ) %>% dplyr::filter(is.finite(xmin), is.finite(xmax), xmax > xmin)
    vlines <- tibble::tibble(period = c(caution_start, risky_start)) %>% dplyr::filter(is.finite(period))
  }

  pal <- MODEL_COLORS[c(MODEL_LABELS[["truth"]], MODEL_LABELS[["sbapc"]])]
  ltys <- MODEL_LINETYPES[c(MODEL_LABELS[["truth"]], MODEL_LABELS[["sbapc"]])]
  g <- ggplot2::ggplot(sum_df, ggplot2::aes(x = period, y = med, color = series, linetype = series))
  if (nrow(rects)) {
    for (i in seq_len(nrow(rects))) {
      g <- g + ggplot2::geom_rect(
        data = rects[i, ],
        ggplot2::aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
        inherit.aes = FALSE,
        fill = rects$fill[i],
        alpha = 0.45,
        color = NA
      )
    }
  }
  if (isTRUE(show_band)) {
    g <- g + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = p10, ymax = p90, fill = series),
      alpha = 0.10,
      color = NA,
      show.legend = FALSE
    )
  }
  g +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.35, color = "gray60") +
    ggplot2::geom_vline(
      data = vlines,
      ggplot2::aes(xintercept = period),
      inherit.aes = FALSE,
      linetype = "dotted",
      color = "gray45",
      linewidth = 0.35
    ) +
    ggplot2::geom_line(linewidth = 0.85, na.rm = TRUE) +
    ggplot2::facet_grid(design ~ scenario_label, scales = "free_y") +
    ggplot2::scale_color_manual(values = pal, breaks = names(pal)) +
    ggplot2::scale_linetype_manual(values = ltys, breaks = names(ltys)) +
    ggplot2::scale_fill_manual(values = pal, breaks = names(pal)) +
    ggplot2::labs(x = "Year", y = "Annual mortality effect (deaths relative to freeze)",
                  color = "Series", linetype = "Series") +
    theme_paper_main(base_size = base_size) +
    ggplot2::theme(
      legend.position = "bottom",
      strip.background = ggplot2::element_rect(fill = "gray95"),
      panel.spacing = grid::unit(0.65, "lines")
    )
}

summarise_misspecification_summary <- function(effect_df, include_full = TRUE) {
  region_df <- effect_df %>%
    dplyr::filter(as.character(model) == MODEL_LABELS[["sbapc"]]) %>%
    add_horizon_region_rows(include_full = include_full) %>%
    dplyr::filter(!is.na(horizon_region))

  seed_level <- region_df %>%
    dplyr::group_by(dgp, scenario, scenario_label, horizon_region, seed) %>%
    dplyr::summarise(
      true_cumulative = sum(delta_truth, na.rm = TRUE),
      estimated_cumulative = sum(delta_hat, na.rm = TRUE),
      recovery_ratio = estimated_cumulative / dplyr::if_else(abs(true_cumulative) > 1e-9, true_cumulative, NA_real_),
      annual_mare_pct = 100 * sum(abs(delta_hat - delta_truth), na.rm = TRUE) /
        pmax(sum(abs(delta_truth), na.rm = TRUE), 1e-9),
      sign_agreement_pct = {
        keep <- is.finite(delta_truth) & abs(delta_truth) > 1e-6 & is.finite(delta_hat)
        if (any(keep)) mean(sign(delta_hat[keep]) == sign(delta_truth[keep])) * 100 else NA_real_
      },
      .groups = "drop"
    )

  seed_level %>%
    dplyr::group_by(dgp, scenario, scenario_label, horizon_region) %>%
    dplyr::summarise(
      design = dgp_public_label(dgp[1]),
      seeds = dplyr::n_distinct(seed),
      annual_mare_pct_mean = mean(annual_mare_pct, na.rm = TRUE),
      recovery_ratio_median = stats::median(recovery_ratio, na.rm = TRUE),
      sign_agreement_pct_mean = mean(sign_agreement_pct, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(factor(as.character(dgp), levels = c("spec_linear", "misspec_tanh")),
                   factor(as.character(scenario), levels = c("up1pc", "down1pc", "quit")),
                   factor(as.character(horizon_region), levels = c("Credible", "Caution", "Risky", "Full horizon")))
}

export_misspecification_summary_table <- function(summary_df,
                                                  csv_out = file.path(OUT_APPENDIX, "tab_misspecification_summary.csv"),
                                                  tex_out = file.path(OUT_APPENDIX, "tab_misspecification_summary.tex")) {
  table_df <- summary_df %>%
    dplyr::filter(as.character(horizon_region) %in% c("Credible", "Caution", "Risky"))
  readr::write_csv(table_df, csv_out)
  tab <- table_df %>%
    dplyr::mutate(scenario_tex = compact_scenario_label(scenario))
  lines <- c(latex_table_open("lllrrr"),
             "Design & Scenario & Horizon & MARE (\\%) & Recovery & Sign (\\%) \\\\",
             "\\midrule")
  last_design <- NULL
  last_scenario <- NULL
  for (i in seq_len(nrow(tab))) {
    row <- tab[i, ]
    design <- if (!identical(last_design, as.character(row$design))) row$design else ""
    scenario <- if (!identical(last_design, as.character(row$design)) ||
                    !identical(last_scenario, as.character(row$scenario))) row$scenario_tex else ""
    if (!is.null(last_design) && !identical(last_design, as.character(row$design))) {
      lines <- c(lines, "\\midrule")
    }
    lines <- c(lines, sprintf(
      "%s & %s & %s & %.1f & %.2f & %.1f \\\\",
      design, scenario, as.character(row$horizon_region),
      row$annual_mare_pct_mean, row$recovery_ratio_median, row$sign_agreement_pct_mean
    ))
    last_design <- as.character(row$design)
    last_scenario <- as.character(row$scenario)
  }
  lines <- c(lines, latex_table_close(
    "The table summarizes recovery of mortality scenario effects under the well-specified and misspecified transmission designs. MARE is averaged across seeds; recovery and sign agreement are computed for cumulative scenario effects."
  ))
  writeLines(lines, tex_out, useBytes = TRUE)
  invisible(tab)
}

write_misspecification_interpretation_notes <- function(summary_df = tibble::tibble(),
                                                        file_out = file.path(OUT_APPENDIX, "misspecification_interpretation_notes.md")) {
  if (!nrow(summary_df)) {
    lines <- c(
      "# Misspecification Interpretation Notes",
      "",
      "The misspecification robustness outputs were not generated because the `misspec_tanh` simulation results were not available in the current result directory.",
      "",
      "To generate them, run the misspecified design with the same seed set and then rerun `generate_appendix_c()`."
    )
    writeLines(lines, file_out, useBytes = TRUE)
    return(invisible(file_out))
  }
  full_horizon <- summary_df %>%
    dplyr::filter(as.character(horizon_region) == "Full horizon") %>%
    dplyr::group_by(design) %>%
    dplyr::summarise(
      annual_mare_pct = mean(annual_mare_pct_mean, na.rm = TRUE),
      recovery_ratio = mean(recovery_ratio_median, na.rm = TRUE),
      sign_agreement_pct = mean(sign_agreement_pct_mean, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      design = factor(design, levels = unname(DGP_LABELS[c("spec_linear", "misspec_tanh")]))
    ) %>%
    dplyr::arrange(design)
  design_txt <- paste(
    sprintf("%s: MARE %.1f%%, recovery %.2f, sign agreement %.1f%%",
            full_horizon$design, full_horizon$annual_mare_pct,
            full_horizon$recovery_ratio, full_horizon$sign_agreement_pct),
    collapse = "; "
  )
  lines <- c(
    "# Misspecification Interpretation Notes",
    "",
    "The robustness exercise compares SBAPC recovery of mortality scenario effects under the well-specified design and under a monotone but misspecified transmission rule.",
    sprintf("Full-horizon averages across scenarios are: %s.", design_txt),
    "Interpretation should focus on degradation relative to the well-specified design, preservation of scenario ordering and signs, and whether the largest discrepancies are concentrated in the lower-support projection horizons."
  )
  writeLines(lines, file_out, useBytes = TRUE)
  invisible(file_out)
}

generate_misspecification_products <- function(seeds = CANONICAL_SEEDS) {
  spec_seeds <- available_result_seeds(OUT_RAW, dgp = "spec_linear")
  misspec_seeds <- available_result_seeds(OUT_RAW, dgp = "misspec_tanh")
  common_seeds <- sort(Reduce(intersect, list(seeds, spec_seeds, misspec_seeds)))
  if (!length(common_seeds)) {
    write_misspecification_interpretation_notes()
    return(invisible(list(status = "missing_misspecification_results", seeds = integer(0))))
  }

  spec_data <- extract_all_metrics(
    seeds = common_seeds,
    dgps = "spec_linear",
    raw_dir = OUT_RAW,
    cache_suffix = NULL
  )
  spec_data <- filter_extracted_data(spec_data, seeds = common_seeds, dgps = "spec_linear")
  misspec_data <- extract_all_metrics(
    seeds = common_seeds,
    dgps = "misspec_tanh",
    raw_dir = OUT_RAW,
    cache_suffix = "misspec_tanh"
  )
  data <- bind_extracted_data(spec_data, misspec_data)
  data <- filter_extracted_data(data, seeds = common_seeds, dgps = c("spec_linear", "misspec_tanh"))
  effect_total <- build_mortality_scenario_effects(data = data, sex_scope = "total")

  g <- plot_misspecification_scenario_recovery(
    effect_total,
    show_band = length(common_seeds) >= 5,
    base_size = paper_fig_base_size("misspecification")
  )
  save_profiled_plot(g, file.path(OUT_APPENDIX, "fig_misspecification_scenario_recovery"),
                     key = "misspecification", bg = "white")

  misspec_summary <- summarise_misspecification_summary(effect_total)
  export_misspecification_summary_table(misspec_summary)
  write_misspecification_interpretation_notes(misspec_summary)

  invisible(list(status = "generated", seeds = common_seeds, summary = misspec_summary))
}

write_seed_level_figure_recommendation <- function(file_out = file.path(OUT_APPENDIX, "seed_level_figure_recommendation.md")) {
  lines <- c(
    "# Seed-Level Figure Recommendation",
    "",
    "Recommendation: retain no more than two seed-level figures in Appendix C. Keep them explicitly illustrative and do not use them as evidence for average performance.",
    "",
    "## Retain",
    "",
    "### `fig_case_study_median_s9.svg`",
    "Shows a representative single-seed trajectory diagnostic for the quit scenario. It is useful as a concrete visual complement to aggregate recovery figures. It should remain in Appendix C only.",
    "",
    "## Main Text, Not Duplicated In Appendix C",
    "",
    "- `fig_transmission_map_support_compare_seed4_M.svg`: shows the prevalence-to-effective-exposure-to-incidence-to-mortality chain for one male seed, including Truth, Observed-window SBAPC, and Full-support SBAPC. It is useful enough for the main text support-window discussion and should not be duplicated in Appendix C.",
    "",
    "## Drop Or Keep As Internal Diagnostics",
    "",
    "- `fig_scenario_atlas_seed4_M.svg` and `fig_scenario_atlas_seed4_F.svg`: visually rich but redundant with the aggregate scenario-effect recovery figure.",
    "- `fig_waterfall_seed4.svg`: useful for internal explanation, but the transmission-map figure is a more direct chain diagnostic.",
    "- `fig_sensitivity_seed4.svg`: single-seed scenario sensitivity is redundant once scenario-effect recovery is aggregated across seeds.",
    "- `fig_transmission_map_seed4_M.svg`: superseded by the support-comparison transmission map if that diagnostic is retained.",
    "- `fig_transmission_map_support_compare_seed4_F.svg`: substantively redundant with the male pathway illustration for the current narrative; do not include it in Appendix C unless the text later makes sex-specific pathway differences central.",
    "- `fig_case_study_best_s26.svg` and `fig_case_study_worst_s41.svg`: useful internally for stress-testing, but too anecdotal for the supplement unless the text explicitly discusses heterogeneity across seeds."
  )
  writeLines(lines, file_out, useBytes = TRUE)
  invisible(file_out)
}

write_float_inventories <- function() {
  seed_label <- seed_count_label()
  section4 <- c(
    "# Section 4 Float Inventory",
    "",
    "Recommended main-text set: one central scenario-effect recovery figure and one compact cumulative recovery table. The chain-recovery table is useful if the text explicitly discusses the sequential mechanism; otherwise it can move to Appendix C.",
    "",
    "| Filename | Document | Priority | Seed aggregation | SVG+PDF | Source note | Purpose |",
    "|---|---|---|---|---|---|---|",
    sprintf("| `fig_scenario_effect_recovery` | Main text | Essential | Aggregated across %s | Yes | Yes | Shows recovery of mortality scenario effects relative to freeze for Truth, SBAPC, and the scenario-blind BAPC benchmark. |", seed_label),
    sprintf("| `tab_cumulative_scenario_recovery` | Main text | Essential | Aggregated across %s | Not applicable | Yes | Summarizes cumulative mortality-effect recovery by scenario and endogenous horizon region. |", seed_label),
    "| `fig_transmission_map_support_compare_seed4_M` | Main text | Useful | Illustrative one-seed diagnostic | Yes | Yes | Visualizes the smoking-to-mortality pathway and the observed-window/full-support contrast for Male. |",
    sprintf("| `tab_chain_recovery` | Main text or Appendix C | Useful | Aggregated across %s | Not applicable | Yes | Checks whether the freeze-baseline sequential chain is recovered at intermediate and mortality levels. |", seed_label),
    sprintf("| `tab_bias_summary` | Appendix C or omit | Optional | Aggregated across %s | Not applicable | Yes | Compact historical/projection bias summary by scenario and sex. |", seed_label)
  )
  appendix <- c(
    "# Appendix C Float Inventory",
    "",
    "| Filename | Document | Priority | Seed aggregation | SVG+PDF | Source note | Purpose |",
    "|---|---|---|---|---|---|---|",
    sprintf("| `fig_scenario_effect_recovery_bysex` | Appendix C | Essential | Aggregated across %s | Yes | Yes | Shows by-sex scenario-effect recovery, including the incidence-anchored diagnostic variant. |", seed_label),
    sprintf("| `fig_support_window_comparison` | Appendix C | Useful | Aggregated across %s | Yes | Yes | Compares Truth, Full-support SBAPC, and Observed-window SBAPC for mortality scenario effects. |", seed_label),
    sprintf("| `tab_support_window_comparison` | Appendix C | Useful | Aggregated across %s | Not applicable | Yes | Quantifies the observed-window penalty by horizon region. |", seed_label),
    sprintf("| `fig_misspecification_scenario_recovery` | Appendix C | Useful | Aggregated across %s | Yes | Yes | Assesses degradation under the Misspecified transmission design. |", seed_label),
    sprintf("| `tab_misspecification_summary` | Appendix C | Useful | Aggregated across %s | Not applicable | Yes | Compact numerical summary of misspecification performance. |", seed_label),
    sprintf("| `fig_reliability_calibration` | Appendix C | Useful | Aggregated across %s | Yes | Yes | Calibration diagnostic for predictive summaries; secondary to scenario-effect recovery. |", seed_label),
    "| `fig_case_study_median_s9` | Appendix C | Optional | Illustrative one-seed diagnostic | Yes | Yes | Shows a representative trajectory case study. |",
    sprintf("| `fig_bias_distributions` | Appendix C | Optional | Aggregated across %s | Yes | Yes | Shows bias dispersion across simulations. |", seed_label),
    "| `fig_transmission_map_support_compare_seed4_F` | Not recommended | Optional | Illustrative one-seed diagnostic | Yes if retained | No current Appendix C note | Female counterpart reviewed as substantively redundant for the current supplement narrative. |"
  )
  writeLines(section4, file.path(OUT_SEC4, "section4_float_inventory.md"), useBytes = TRUE)
  writeLines(appendix, file.path(OUT_APPENDIX, "appendixC_float_inventory.md"), useBytes = TRUE)
  invisible(TRUE)
}

generate_scenario_effect_products <- function(data = NULL) {
  if (is.null(data)) data <- extract_all_metrics()

  effect_total <- build_mortality_scenario_effects(data = data, sex_scope = "total")
  effect_bysex <- build_mortality_scenario_effects(data = data, sex_scope = "by_sex")
  summary_total <- summarise_scenario_effect_recovery(effect_total)
  summary_bysex <- summarise_scenario_effect_recovery(effect_bysex)
  horizon_boundary_audit(data)

  readr::write_csv(effect_total, file.path(OUT_SEC4, "scenario_effect_recovery_detail_total.csv"))
  readr::write_csv(summary_total, file.path(OUT_SEC4, "scenario_effect_recovery_summary.csv"))
  readr::write_csv(effect_bysex, file.path(OUT_APPENDIX, "scenario_effect_recovery_detail_bysex.csv"))
  readr::write_csv(summary_bysex, file.path(OUT_APPENDIX, "scenario_effect_recovery_summary_bysex.csv"))

  export_scenario_effect_recovery_table(
    summary_total,
    file.path(OUT_SEC4, "tab_scenario_effect_recovery.tex")
  )
  cumulative_summary <- summarise_cumulative_scenario_recovery(effect_total)
  export_cumulative_scenario_recovery_table(cumulative_summary)

  chain_df <- build_chain_recovery_data()
  chain_summary <- summarise_chain_recovery(chain_df)
  export_chain_recovery_table(chain_summary)

  g_main <- plot_scenario_effect_recovery(
    effect_total,
    include_models = unname(MODEL_LABELS[c("sbapc", "bapc")]),
    base_size = paper_fig_base_size("scenario_effect")
  )
  save_profiled_plot(g_main, file.path(OUT_SEC4, "fig_scenario_effect_recovery"), key = "scenario_effect", bg = "white")

  g_bysex <- plot_scenario_effect_recovery(
    effect_bysex,
    include_models = unname(MODEL_LABELS[c("sbapc", "sbapc_no_prev", "bapc")]),
    base_size = paper_fig_base_size("scenario_effect_bysex")
  )
  save_profiled_plot(g_bysex, file.path(OUT_APPENDIX, "fig_scenario_effect_recovery_bysex"), key = "scenario_effect_bysex", bg = "white")

  invisible(list(effect_total = effect_total, effect_bysex = effect_bysex,
                 summary_total = summary_total, summary_bysex = summary_bysex,
                 cumulative_summary = cumulative_summary, chain_summary = chain_summary))
}

plot_reliability_calibration <- function(data, base_size = paper_fig_base_size("reliability")) {
  # data$inc contains the errors per period/sex/seed/dgp
  
  df <- data$inc %>%
    dplyr::filter(period > 2022) %>%
    dplyr::mutate(
      horizon = period - 2022,
      abs_rel_error_pct = abs(rel_error) * 100
    )
  
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
  max_y <- suppressWarnings(max(df$abs_rel_error_pct, na.rm = TRUE))
  if (!is.finite(max_y)) max_y <- 1
  year_label_y <- max_y * 0.04

  g <- ggplot() +
    geom_rect(data = rects, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill), alpha = 0.5) +
    stat_summary(data = df, aes(x = horizon, y = abs_rel_error_pct), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.2, fill = "blue") +
    stat_summary(data = df, aes(x = horizon, y = abs_rel_error_pct), fun = "mean", geom = "line", linewidth = 1, color = "blue") +
    expand_limits(y = 0) +
    scale_fill_identity()
    
  if (length(vlines) > 0) {
    g <- g + geom_vline(xintercept = vlines, linetype = "dashed", color = "gray40")
    # Add labels for categories
    for(i in seq_along(rects$label)) {
      g <- g + annotate("text", x = (rects$xmin[i] + rects$xmax[i])/2, y = Inf, label = rects$label[i], vjust = 2, size = 3.5, fontface = "italic", color = "gray30")
    }
    for(hv in vlines) {
       g <- g + annotate(
         "text",
         x = hv,
         y = year_label_y,
         label = as.character(hv + 2022),
         angle = 90,
         hjust = 0,
         vjust = -0.35,
         size = 3,
         color = "gray25"
       )
    }
  }

  g <- g + labs(
         y = "Mean absolute relative error (%)",
         x = "Projection horizon (years)"
       ) +
    theme_paper_main(base_size = base_size)
  
  return(g)
}


# =============================================================
# 4. MAIN OUTPUT GENERATOR
# =============================================================

replicate_main_paper <- function() {
  # 1. Scenario Atlas (By Sex)
  for (sx in c("M", "F")) {
    g_atlas <- plot_scenario_atlas_by_sex(seed = 4, sex_lab = sx)
    save_profiled_plot(g_atlas, file.path(OUT_SEC4, sprintf("fig_scenario_atlas_seed4_%s", sx)), key = "scenario_atlas", bg = "white")
  }
  
  # 2. Waterfall
  g2 <- plot_transmission_waterfall(seed = 4, dgp = "spec_linear", scen = "quit")
  save_profiled_plot(g2, file.path(OUT_SEC4, "fig_waterfall_seed4"), key = "waterfall", bg = "white")
  
  # 3. Scenario Sensitivity (New)
  g_sens <- plot_scenario_sensitivity_informed(seed = 4, dgp = "spec_linear")
  save_profiled_plot(g_sens, file.path(OUT_SEC4, "fig_sensitivity_seed4"), key = "sensitivity", bg = "white")
  
  # 4. Transmission map
  g_map <- plot_transmission_map(seed = 4, dgp = "spec_linear", sex_lab = "M")
  save_profiled_plot(g_map, file.path(OUT_SEC4, "fig_transmission_map_seed4_M"), key = "transmission_map", bg = "white")
  generate_support_transmission_maps(seed = 4, dgp = "spec_linear", force_oracle = FALSE)

  # 5. Bias Table
  data <- extract_all_metrics()
  generate_scenario_effect_products(data)
  export_bias_summary_table(data$metrics)
  write_csv(data$metrics, file.path(OUT_RAW, "all_metrics.csv"))
  
  # 6. Support Summary
  write_csv(data$support, file.path(OUT_SEC4, "support_summary.csv"))
  write_figure_titles_notes("section4")
  write_float_inventories()
  
  message("\nMain paper replication files generated in: ", OUT_SEC4)
}

# =============================================================
# 5. APPENDIX C GENERATOR
# =============================================================

generate_appendix_c <- function() {
  data <- extract_all_metrics()
  
  # 1. Distribution of Bias (Boxplots)
  plot_df <- data$metrics %>%
    dplyr::mutate(
      scenario_f = factor(scenario, levels = c("up1pc", "freeze", "down1pc", "quit"), labels = SCEN_LABELS),
      sex_label = factor(sex_public_label(sex), levels = c("Male", "Female")),
      dgp_label = dgp_public_label(dgp)
    )
  scen_colors_by_label <- stats::setNames(unname(SCEN_COLORS[names(SCEN_LABELS)]), unname(SCEN_LABELS))
  facet_vars <- if (dplyr::n_distinct(plot_df$dgp) > 1) {
    ggplot2::vars(dgp_label, sex_label)
  } else {
    ggplot2::vars(sex_label)
  }
  
  g_bias <- ggplot(plot_df, aes(x = scenario_f, y = proj_bias, fill = scenario_f)) +
    geom_boxplot(alpha = 0.7) +
    facet_wrap(facet_vars) +
    scale_fill_manual(values = scen_colors_by_label) +
    labs(y = "Projection bias (%)", x = "Scenario") +
    theme_paper_main(base_size = paper_fig_base_size("bias_distribution")) + theme(legend.position = "none")
  save_profiled_plot(g_bias, file.path(OUT_APPENDIX, "fig_bias_distributions"), key = "bias_distribution", bg = "white")
  
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
    g_case <- plot_deconstruction_figure(seed = s, dgp = "spec_linear", scen = "quit")
    save_profiled_plot(g_case, file.path(OUT_APPENDIX, sprintf("fig_case_study_%s_s%d", tolower(lbl), s)), key = "case_study", bg = "white")
  }
  
  # 3. Reliability calibration is an appendix diagnostic, not a main-text result.
  g_reliability <- plot_reliability_calibration(data)
  save_profiled_plot(g_reliability, file.path(OUT_APPENDIX, "fig_reliability_calibration"),
                     key = "reliability", bg = "white")

  # 4. Support-window and misspecification diagnostics.
  generate_support_window_products(realistic_data = data)
  generate_misspecification_products()

  # 5. Recommendations and full detailed table (CSV)
  write_seed_level_figure_recommendation()
  write_csv(data$metrics, file.path(OUT_APPENDIX, "full_simulation_matrix.csv"))
  write_figure_titles_notes("appendixC", case_seeds = case_seeds)
  write_float_inventories()
  
  message("\nAppendix C replication files generated in: ", OUT_APPENDIX)
}

# =============================================================
# 6. ORCHESTRATOR
# =============================================================

replicate_all_simulations <- function(seeds = CANONICAL_SEEDS, n_cores = 6, force_rerun = TRUE) {
  message("STARTING FULL REPLICATION WORKFLOW...")
  old_seeds <- CANONICAL_SEEDS
  CANONICAL_SEEDS <<- sort(unique(as.integer(seeds)))
  on.exit({ CANONICAL_SEEDS <<- old_seeds }, add = TRUE)
  
  # 1. Run Simulations (FORCE RERUN to overwrite old files)
  # Using 6 cores to avoid memory allocation errors (INLA is memory intensive)
  run_simulation_replication(seeds = CANONICAL_SEEDS, n_cores = n_cores, force_rerun = force_rerun)
  
  # 2. Section 4
  replicate_main_paper()
  
  # 3. Appendix C
  generate_appendix_c()
  
  message("\nALL REPLICATION TASKS COMPLETED SUCCESSFULLY.")
}

replicate_final_simulations <- function(seeds = CANONICAL_SEEDS,
                                        n_cores = 4,
                                        force_rerun = FALSE,
                                        run_oracle = TRUE,
                                        run_misspec = TRUE) {
  message("STARTING FINAL SIMULATION WORKFLOW...")
  old_seeds <- CANONICAL_SEEDS
  CANONICAL_SEEDS <<- sort(unique(as.integer(seeds)))
  on.exit({ CANONICAL_SEEDS <<- old_seeds }, add = TRUE)

  message("Seed set: ", paste(range(CANONICAL_SEEDS), collapse = "-"),
          " (", length(CANONICAL_SEEDS), " seeds)")
  message("Output base: ", OUT_BASE)
  message("Worker count: ", n_cores)

  run_simulation_replication(
    seeds = CANONICAL_SEEDS,
    dgps = "spec_linear",
    scens = CANONICAL_SCENS,
    force_rerun = force_rerun,
    n_cores = n_cores,
    information_set = "realistic",
    raw_dir = OUT_RAW
  )

  if (isTRUE(run_oracle)) {
    run_simulation_replication(
      seeds = CANONICAL_SEEDS,
      dgps = "spec_linear",
      scens = CANONICAL_SCENS,
      force_rerun = force_rerun,
      n_cores = n_cores,
      information_set = "oracle",
      raw_dir = OUT_RAW_ORACLE
    )
  }

  if (isTRUE(run_misspec)) {
    run_simulation_replication(
      seeds = CANONICAL_SEEDS,
      dgps = "misspec_tanh",
      scens = CANONICAL_SCENS,
      force_rerun = force_rerun,
      n_cores = n_cores,
      information_set = "realistic",
      raw_dir = OUT_RAW
    )
  }

  replicate_main_paper()
  generate_appendix_c()

  removed <- cleanup_inla_temp()
  message("Cleaned INLA temporary directories: ", removed)
  message("\nFINAL SIMULATION WORKFLOW COMPLETED SUCCESSFULLY.")
}
