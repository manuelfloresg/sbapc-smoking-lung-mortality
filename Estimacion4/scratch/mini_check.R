# scratch/mini_check.R
source("runs/_source_all.R")
batch_dir <- "results/batches/20260504_1044_FINAL_PRODUCTION"
seeds <- 1:5
dgps <- c("spec_linear", "misspec_tanh")
scenarios <- c("freeze", "quit")

results_list <- list()

for (dgp_val in dgps) {
  for (scn_val in scenarios) {
    for (seed_val in seeds) {
      f <- file.path(batch_dir, sprintf("res_%s_s%d_%s.rds", dgp_val, seed_val, scn_val))
      if (!file.exists(f)) next
      
      res <- readRDS(f)
      
      old_verbose <- BAPC_VERBOSE
      BAPC_VERBOSE <<- FALSE
      sim <- simulate_PIM_data(cause_id = dgp_val, seed = seed_val, dgp = dgp_val, scenario_name = scn_val, beta_mode = "fixed_rr_offset")
      BAPC_VERBOSE <<- old_verbose
      
      exposure_val <- sim$meta$exposure
      truth_inc <- sim$inc_truth_grid %>% 
        filter(period > 2022) %>%
        summarise(truth_cases = sum(rateI_scen_true * exposure_val)) %>%
        pull(truth_cases)
        
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
      
      results_list[[length(results_list)+1]] <- data.frame(
        dgp = dgp_val,
        seed = seed_val,
        scenario = scn_val,
        bias_pct = bias
      )
    }
  }
}

summary_df <- do.call(rbind, results_list) %>%
  group_by(dgp, scenario) %>%
  summarise(
    n = n(),
    mean_bias = mean(bias_pct, na.rm = TRUE),
    .groups = "drop"
  )

print(summary_df)
