# scratch/check_z_prev.R
source("runs/_runtime_setup.R")
source("runs/_source_all.R")
source("adapters/build_inputs_sim.R")

library(dplyr)
library(ggplot2)

seed <- 4
dgp <- "spec_linear"
cause_id <- "lung"
sex <- "M"

# Load the simulation input
rds_file <- sprintf("results/20260515_ESTIMATE_V15/raw_data/res_%s_s%s_up1pc.rds", dgp, seed)
if (!file.exists(rds_file)) stop("File not found: ", rds_file)

res <- readRDS(rds_file)
# res$resM$rates_all_full contains the data

df <- res$resM$rates_all_full %>%
  filter(age == 65) %>%
  select(period, z_prev, coef_fc_offset_I, rate_hat)

print(head(df))
print(tail(df))

ggplot(df, aes(x = period, y = z_prev)) +
  geom_line() +
  ggtitle(sprintf("z_prev for Age 65, Seed %d, UP1PC", seed))

ggsave("results/20260515_ESTIMATE_V15/scratch_z_prev_check.png")
