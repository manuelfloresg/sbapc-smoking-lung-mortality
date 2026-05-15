# scratch/inspect_apc_components.R
source("runs/_runtime_setup.R")
source("runs/_source_all.R")

library(dplyr)
library(ggplot2)

seed <- 4
dgp <- "spec_linear"
scen <- "up1pc"

rds_file <- sprintf("results/20260515_ESTIMATE_V15/raw_data/res_%s_s%d_%s.rds", dgp, seed, scen)
if (!file.exists(rds_file)) stop("File not found: ", rds_file)

res <- readRDS(rds_file)
fit <- res$resM$inc_fit$fit_inc # The historical fit object

# Period effects
per_re <- fit$summary.random$period_id
# Cohort effects
coh_re <- fit$summary.random$cohort_id

# Let's see the trends
png("results/20260515_ESTIMATE_V15/scratch_apc_trends.png", width=1000, height=500)
par(mfrow=c(1,2))
plot(per_re$ID, per_re$mean, type="l", main="Period Effect (Informed)")
plot(coh_re$ID, coh_re$mean, type="l", main="Cohort Effect (Informed)")
dev.off()

cat("Period slope (last 5 years):", diff(tail(per_re$mean, 5)), "\n")
cat("Cohort slope (last 5 years):", diff(tail(coh_re$mean, 5)), "\n")
