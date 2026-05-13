# debug_sim.R
source("runs/_source_all.R")
BAPC_VERBOSE <- TRUE
cat(">>> DEBUG: srcref(get_risk_reversal_w):\n")
print(attr(get_risk_reversal_w, "srcref"))
print(environment(get_risk_reversal_w))
print(body(get_risk_reversal_w))

# Test parameters
cause_id <- "spec_linear"
seed <- 1
scenario_name <- "quit"
beta_mode <- "fixed_rr_offset"

# Run simulation
sim <- simulate_PIM_data(
  cause_id = cause_id, 
  seed = seed, 
  dgp = cause_id, 
  scenario_name = scenario_name, 
  beta_mode = beta_mode
)

# Inspect results
inc <- sim$inc_truth_grid
z_base <- sim$zI_base_true_used
z_scen <- sim$zI_scen_true_used

cat("\nMale z_prev for 2022-2025 (Base):\n")
print(z_base %>% filter(sex == "M", age == 60, period %in% 2022:2025) %>% select(period, age, z_prev_used))

cat("\nMale z_prev for 2022-2025 (Scen):\n")
print(z_scen %>% filter(sex == "M", age == 60, period %in% 2022:2025) %>% select(period, age, z_prev_used))

cat("\nMale Incidence Rate (True) for 2022-2025:\n")
print(inc %>% filter(sex == "M", age == 60, period %in% 2022:2025) %>% select(period, age, rateI_scen_true))

cat("\nEffective Prevalence (q_eff) for 2022-2025:\n")
# q_eff is not in inc_truth_grid, let's find it in the function scope or elsewhere
# Actually, I can't easily see it from the outside unless I return it.

# Let's check the rate ratio
r2022 <- (inc %>% filter(sex == "M", age == 60, period == 2022))$rateI_scen_true
r2023 <- (inc %>% filter(sex == "M", age == 60, period == 2023))$rateI_scen_true
cat(sprintf("\nRate 2022: %f\nRate 2023: %f\nRatio: %f\n", r2022, r2023, r2023/r2022))
