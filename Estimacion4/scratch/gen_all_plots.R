library(dplyr)
library(ggplot2)
library(readr)

# Load existing helpers
source("R/00_defaults.R")
source("R/31_diagnostics_against_truth.R")
source("runs/replication_diagnostics.R")

# Update Output paths to the new stable folder
OUT_DIR <- "results/20260506_STABLE_V2_PURE"
OUT_SEC4 <- file.path(OUT_DIR, "section4")
OUT_RAW <- file.path(OUT_DIR, "raw_data")

if (!dir.exists(OUT_SEC4)) dir.create(OUT_SEC4, recursive = TRUE)

# 1. New Plotting function for Incidence Deconstruction
plot_incidence_deconstruction <- function(seed = 4, dgp = "spec_linear", scen = "quit") {
  rds_file <- file.path(OUT_RAW, sprintf("res_%s_s%d_%s.rds", dgp, seed, scen))
  if (!file.exists(rds_file)) return(NULL)
  rb <- readRDS(rds_file)
  
  diag_res <- compare_pipeline_to_truth(rb, rb, out_dir = NULL)
  
  df_inc <- diag_res$inc
  if (is.null(df_inc) || nrow(df_inc) == 0) return(NULL)
  
  plot_df <- df_inc %>%
    dplyr::select(period, sex, 
                  Truth = rate_true, 
                  `Informed SBAPC (I | P)` = rate_hat, 
                  `Pure BAPC (I)` = rate_bapc) %>%
    tidyr::pivot_longer(cols = c(Truth, `Informed SBAPC (I | P)`, `Pure BAPC (I)`), 
                        names_to = "Series", values_to = "Rate") %>%
    mutate(Rate = Rate * 100000) # Per 100k
  
  plot_df$Series <- factor(plot_df$Series, levels = c("Truth", "Informed SBAPC (I | P)", "Pure BAPC (I)"))
  
  g <- ggplot(plot_df, aes(x = period, y = Rate, color = Series, linetype = Series)) +
    geom_line(linewidth = 1) +
    geom_vline(xintercept = 2022, linetype = "dotted", color = "gray50") +
    facet_wrap(~sex, scales = "free_y") +
    scale_color_manual(values = c("Truth" = "black", "Informed SBAPC (I | P)" = "#CD5C5C", "Pure BAPC (I)" = "#4682B4")) +
    scale_linetype_manual(values = c("Truth" = "dashed", "Informed SBAPC (I | P)" = "solid", "Pure BAPC (I)" = "dotted")) +
    labs(title = "Incidence Information Gain Deconstruction",
         subtitle = sprintf("Seed %d | DGP: %s | Scenario: %s", seed, dgp, scen),
         y = "Rate per 100,000", x = "Year") +
    theme_minimal() + theme(legend.position = "bottom")
  
  return(g)
}

# 2. Loop through all scenarios for Seed 4
scenarios <- c("quit", "up1pc", "down1pc", "freeze")

for (scen in scenarios) {
  message("Generating plots for: ", scen)
  
  # Mortality
  g_mort <- plot_deconstruction_figure(seed = 4, dgp = "spec_linear", scen = scen)
  if (!is.null(g_mort)) {
    ggsave(file.path(OUT_SEC4, sprintf("fig_deconstruction_seed4_%s.png", scen)), g_mort, width = 10, height = 6, bg = "white")
  }
  
  # Incidence
  g_inc <- plot_incidence_deconstruction(seed = 4, dgp = "spec_linear", scen = scen)
  if (!is.null(g_inc)) {
    ggsave(file.path(OUT_SEC4, sprintf("fig_incidence_deconstruction_seed4_%s.png", scen)), g_inc, width = 10, height = 6, bg = "white")
  }
}

message("DONE")
