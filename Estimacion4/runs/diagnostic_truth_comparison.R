# runs/diagnostic_truth_comparison.R
# =============================================================
# Diagnostic: Compare "Truth" (DGP) trajectories across all scenarios
# =============================================================

source("runs/_runtime_setup.R")
source("runs/_source_all.R")
source("adapters/build_inputs_sim.R")

library(dplyr)
library(ggplot2)
library(tidyr)

# Configuration
SEED <- 4
DGP  <- "spec_linear"
SCENARIOS <- c("freeze", "up1pc", "down1pc", "quit")
CAUSE_ID  <- "lung"

message("Generating Truth comparison for Seed ", SEED, " | DGP: ", DGP)

truth_list <- list()

for (scen in SCENARIOS) {
  message("  Simulating: ", scen)
  sim <- simulate_PIM_data(cause_id = CAUSE_ID, seed = SEED, dgp = DGP, scenario_name = scen, rr_inc = 15)
  
  # Extract mortality truth
  m_true <- sim$mort_truth_grid %>%
    dplyr::filter(age >= 35, age <= 89) %>%
    dplyr::group_by(period, sex) %>%
    dplyr::summarise(deaths_true = sum(mort_deaths_scen_true, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(scenario = scen)
  
  truth_list[[scen]] <- m_true
}

df_plot <- dplyr::bind_rows(truth_list)
df_plot$scenario <- factor(df_plot$scenario, levels = SCENARIOS)

# Create Plot
last_hist <- 2022

g <- ggplot(df_plot, aes(x = period, y = deaths_true, color = scenario, linetype = scenario)) +
  geom_line(linewidth = 1.2) +
  geom_vline(xintercept = last_hist, linetype = "dotted", color = "gray40") +
  facet_wrap(~sex, scales = "free_y") +
  scale_color_manual(values = c(
    "freeze"  = "#d62728", # red
    "up1pc"   = "#9467bd", # purple
    "down1pc" = "#ff7f0e", # orange
    "quit"    = "#1f77b4"  # blue
  )) +
  labs(title = "DGP Ground Truth Comparison: All Scenarios",
       subtitle = sprintf("Seed %d | DGP: %s | Historical + Projection", SEED, DGP),
       y = "Annual Deaths (Truth)", x = "Year",
       color = "Scenario", linetype = "Scenario") +
  theme_minimal(base_family = "sans") +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold", size = 11))

# Save
out_dir <- "results/diagnostics_truth"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(out_dir, sprintf("truth_comparison_s%d_%s.png", SEED, DGP))
ggsave(out_file, g, width = 10, height = 6, bg = "white")

message("Diagnostic plot saved to: ", out_file)
