# generate_final_atlas.R
source("runs/_runtime_setup.R")
source("runs/replication_diagnostics.R")

message(">>> GENERATING FINAL DIAGNOSTIC ATLAS <<<")
message("Base Directory: ", OUT_BASE)

# 1. Section 4 Figures and Table
replicate_main_paper()

# 2. Appendix C Detailed Metrics and Case Studies
generate_appendix_c()

message(">>> ATLAS GENERATION COMPLETED <<<")
