# scratch/final_bias_check.R
source("runs/_source_all.R")
batch_dir <- "results/batches/20260504_1044_FINAL_PRODUCTION"
files <- list.files(batch_dir, pattern = "^res_.*\\.rds$", full.names = TRUE)

if (length(files) == 0) stop("No files found in batch dir")

results_list <- list()
cat("Processing", length(files), "files...\n")

for (f in files) {
  res <- readRDS(f)
  fname <- basename(f)
  parts <- strsplit(gsub("\\.rds$", "", fname), "_")[[1]]
  
  # Format: res_{dgp}_s{seed}_{scn}.rds
  # res_spec_linear_s1_freeze.rds -> parts=["res", "spec", "linear", "s1", "freeze"]
  # res_misspec_tanh_s1_freeze.rds -> parts=["res", "misspec", "tanh", "s1", "freeze"]
  
  if (parts[2] == "spec") {
    dgp_val <- "spec_linear"
    seed_str <- parts[4]
    scn_val <- parts[5]
  } else {
    dgp_val <- "misspec_tanh"
    seed_str <- parts[4]
    scn_val <- parts[5]
  }
  seed_val <- as.integer(gsub("s", "", seed_str))

  # REGENERATE TRUTH
  old_verbose <- BAPC_VERBOSE
  BAPC_VERBOSE <<- FALSE
  sim <- simulate_PIM_data(cause_id = dgp_val, seed = seed_val, dgp = dgp_val, scenario_name = scn_val, beta_mode = "fixed_rr_offset")
  BAPC_VERBOSE <<- old_verbose
  
  exposure_val <- sim$meta$exposure
  
  truth_inc <- sim$inc_truth_grid %>% 
    filter(period > 2022) %>%
    summarise(truth_cases = sum(rateI_scen_true * exposure_val)) %>%
    pull(truth_cases)
    
  # Extract prediction
  pred_m <- if (!is.null(res$resM)) {
    res$resM$inc_fit$rates_all %>%
      filter(period > 2022) %>%
      left_join(sim$pop_all %>% filter(sex == "M"), by = c("sex", "age", "period")) %>%
      summarise(cases = sum(rate_hat * exposure)) %>% pull(cases)
  } else 0
  
  pred_f <- if (!is.null(res$resF)) {
    res$resF$inc_fit$rates_all %>%
      filter(period > 2022) %>%
      left_join(sim$pop_all %>% filter(sex == "F"), by = c("sex", "age", "period")) %>%
      summarise(cases = sum(rate_hat * exposure)) %>% pull(cases)
  } else 0
  
  pred_cases <- pred_m + pred_f
  bias <- (pred_cases / truth_inc - 1) * 100
  
  results_list[[f]] <- data.frame(
    dgp = dgp_val,
    seed = seed_val,
    scenario = scn_val,
    bias_pct = bias
  )
  if (length(results_list) %% 20 == 0) cat(".")
}
cat("\n")

summary_df <- do.call(rbind, results_list) %>%
  group_by(dgp, scenario) %>%
  summarise(
    n = n(),
    mean_bias = mean(bias_pct, na.rm = TRUE),
    sd_bias = sd(bias_pct, na.rm = TRUE),
    .groups = "drop"
  )

print(summary_df)
