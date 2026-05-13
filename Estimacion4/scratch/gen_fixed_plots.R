library(dplyr)
library(ggplot2)
library(readr)

source("R/00_defaults.R")
source("R/31_diagnostics_against_truth.R")
source("runs/replication_diagnostics.R")

OUT_DIR <- "results/20260506_STABLE_V2_PURE"
OUT_SEC4 <- file.path(OUT_DIR, "section4")
OUT_RAW <- file.path(OUT_DIR, "raw_data")

# 1. Extraer benchmarks INVARIANTES (usaremos los del archivo freeze como referencia, 
# pero sabemos que el BAPC puro es el mismo en todos)
seed <- 4
dgp <- "spec_linear"

# Cargamos el freeze solo para sacar los benchmarks (Azul y Naranja)
res_ref <- readRDS(file.path(OUT_RAW, sprintf("res_%s_s%d_freeze.rds", dgp, seed)))
diag_ref <- compare_pipeline_to_truth(res_ref, res_ref, out_dir = NULL)

# Benchmarks de Mortalidad (fijos)
bench_mort <- diag_ref$mort %>%
  dplyr::select(period, sex, 
                `Pure BAPC (M)` = deaths_bapc, 
                `Uninformed SBAPC (M | I)` = deaths_noP)

# Benchmarks de Incidencia (fijos)
bench_inc <- diag_ref$inc %>%
  dplyr::select(period, sex, 
                `Pure BAPC (I)` = rate_bapc)

# 2. Función de ploteo que usa benchmarks EXTERNOS fijos
plot_deconstruction_fixed <- function(scen, bench_m, seed=4, dgp="spec_linear") {
  rds_file <- file.path(OUT_RAW, sprintf("res_%s_s%d_%s.rds", dgp, seed, scen))
  rb <- readRDS(rds_file)
  diag_scen <- compare_pipeline_to_truth(rb, rb, out_dir = NULL)
  
  df_mort <- diag_scen$mort %>%
    dplyr::select(period, sex, Truth = deaths_true, `Informed SBAPC (M | I | P)` = deaths_hat) %>%
    dplyr::left_join(bench_m, by = c("period", "sex"))
    
  plot_df <- df_mort %>%
    tidyr::pivot_longer(cols = -c(period, sex), names_to = "Series", values_to = "Deaths")
  
  plot_df$Series <- factor(plot_df$Series, levels = c("Truth", "Informed SBAPC (M | I | P)", "Uninformed SBAPC (M | I)", "Pure BAPC (M)"))
  
  ggplot(plot_df, aes(x = period, y = Deaths, color = Series, linetype = Series)) +
    geom_line(linewidth = 1) +
    geom_vline(xintercept = 2022, linetype = "dotted", color = "gray50") +
    facet_wrap(~sex, scales = "free_y") +
    scale_color_manual(values = c("Truth"="black", "Informed SBAPC (M | I | P)"="#CD5C5C", "Uninformed SBAPC (M | I)"="#ff7f0e", "Pure BAPC (M)"="#4682B4")) +
    scale_linetype_manual(values = c("Truth"="dashed", "Informed SBAPC (M | I | P)"="solid", "Uninformed SBAPC (M | I)"="dotdash", "Pure BAPC (M)"="dotted")) +
    labs(title = paste("Mortality Deconstruction - Scenario:", scen), subtitle = sprintf("Seed %d | Benchmarks Fixed", seed)) +
    theme_minimal() + theme(legend.position = "bottom")
}

plot_incidence_fixed <- function(scen, bench_i, seed=4, dgp="spec_linear") {
  rds_file <- file.path(OUT_RAW, sprintf("res_%s_s%d_%s.rds", dgp, seed, scen))
  rb <- readRDS(rds_file)
  diag_scen <- compare_pipeline_to_truth(rb, rb, out_dir = NULL)
  
  df_inc <- diag_scen$inc %>%
    dplyr::select(period, sex, Truth = rate_true, `Informed SBAPC (I | P)` = rate_hat) %>%
    dplyr::left_join(bench_i, by = c("period", "sex"))
    
  plot_df <- df_inc %>%
    tidyr::pivot_longer(cols = -c(period, sex), names_to = "Series", values_to = "Rate") %>%
    mutate(Rate = Rate * 100000)
    
  plot_df$Series <- factor(plot_df$Series, levels = c("Truth", "Informed SBAPC (I | P)", "Pure BAPC (I)"))
  
  ggplot(plot_df, aes(x = period, y = Rate, color = Series, linetype = Series)) +
    geom_line(linewidth = 1) +
    geom_vline(xintercept = 2022, linetype = "dotted", color = "gray50") +
    facet_wrap(~sex, scales = "free_y") +
    scale_color_manual(values = c("Truth"="black", "Informed SBAPC (I | P)"="#CD5C5C", "Pure BAPC (I)"="#4682B4")) +
    scale_linetype_manual(values = c("Truth"="dashed", "Informed SBAPC (I | P)"="solid", "Pure BAPC (I)"="dotted")) +
    labs(title = paste("Incidence Deconstruction - Scenario:", scen), y = "Rate per 100k") +
    theme_minimal() + theme(legend.position = "bottom")
}

# Generar 
for (scen in c("quit", "up1pc", "down1pc", "freeze")) {
  message("Plotting: ", scen)
  g_m <- plot_deconstruction_fixed(scen, bench_mort)
  ggsave(file.path(OUT_SEC4, sprintf("fig_decon_mort_fixed_%s.png", scen)), g_m, width=10, height=6, bg="white")
  
  g_i <- plot_incidence_fixed(scen, bench_inc)
  ggsave(file.path(OUT_SEC4, sprintf("fig_decon_inc_fixed_%s.png", scen)), g_i, width=10, height=6, bg="white")
}
message("DONE")
