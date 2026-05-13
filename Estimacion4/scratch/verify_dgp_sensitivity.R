# Scratch script to verify DGP sensitivity (v2)
source("runs/_source_all.R")

# Simulation settings
seed <- 1
dgp_type <- "spec_linear"
scenario <- "down1pc"
base_year <- 2022

message(">>> Simulating data for seed ", seed, " scenario ", scenario)
sim <- simulate_PIM_data(
  cause_id = "lung", 
  seed = seed, 
  dgp = dgp_type, 
  scenario_name = scenario,
  beta_mode = "fixed_rr_offset"
)

# Inspect the truth grid
truth <- sim$inc_truth_grid

# Extract rates for 2022 and 2050 for a specific age/sex
# Males, age 60
val_2022 <- truth %>% filter(sex == "M", age == 60, period == 2022)
val_2050 <- truth %>% filter(sex == "M", age == 60, period == 2050)
val_2050_base <- truth %>% filter(sex == "M", age == 60, period == 2050) # wait, I need the freeze scenario for comparison

# To compare with freeze, I need to simulate freeze too
sim_freeze <- simulate_PIM_data(
  cause_id = "lung", 
  seed = seed, 
  dgp = dgp_type, 
  scenario_name = "freeze",
  beta_mode = "fixed_rr_offset"
)
val_2050_freeze <- sim_freeze$inc_truth_grid %>% filter(sex == "M", age == 60, period == 2050)

message("\n--- INSIGHTS (Year 2050) ---")
message("  rateI_freeze_true: ", round(val_2050_freeze$rateI_scen_true * 1e5, 2))
message("  rateI_down1_true:   ", round(val_2050$rateI_scen_true * 1e5, 2))
message("  Reduction (%):     ", round(100 * (1 - val_2050$rateI_scen_true / val_2050_freeze$rateI_scen_true), 2), "%")


# Check if it dropped relative to base (unperturbed) in 2023
if (val_2023$rateI_scen_true < val_2023$rateI_base_true) {
  message("SUCCESS: rateI_scen_true is lower than rateI_base_true in 2023!")
} else {
  message("FAILURE: rateI_scen_true is NOT lower than rateI_base_true in 2023!")
}
