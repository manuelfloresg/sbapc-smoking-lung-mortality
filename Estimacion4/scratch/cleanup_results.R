
dirs <- list.dirs("results", recursive = FALSE)
keep <- "results/20260507_ESTIMATE_V14"
to_delete <- dirs[grepl("^results/202605", dirs) & dirs != keep]
cat("Deleting:\n", paste(to_delete, collapse="\n"), "\n")
unlink(to_delete, recursive = TRUE, force = TRUE)
cat("Cleanup complete.\n")
