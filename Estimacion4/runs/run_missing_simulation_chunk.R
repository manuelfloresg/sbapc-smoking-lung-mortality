# Sequential worker for missing simulation top-ups.
#
# Intended use:
#   Rscript runs/run_missing_simulation_chunk.R mode=both seeds=1-13
#
# The script does not use future/multisession. It is safe to launch several
# independent instances from Windows, each with a disjoint seed range.

args <- commandArgs(trailingOnly = TRUE)
kv <- strsplit(args, "=", fixed = TRUE)
opts <- stats::setNames(vapply(kv, function(x) if (length(x) >= 2L) x[[2L]] else "", character(1)),
                        vapply(kv, function(x) x[[1L]], character(1)))

get_opt <- function(name, default = "") {
  if (name %in% names(opts) && nzchar(opts[[name]])) opts[[name]] else default
}

Sys.setenv(BAPC_OUT_BASE = get_opt("out_base", Sys.getenv("BAPC_OUT_BASE", "results/20260518_PROD_CANDIDATE")))
source("runs/replication_diagnostics.R")

parse_seed_spec <- function(x) {
  if (is.null(x) || !nzchar(x)) return(CANONICAL_SEEDS)
  pieces <- unlist(strsplit(x, ",", fixed = TRUE), use.names = FALSE)
  out <- integer(0)
  for (piece in pieces) {
    piece <- trimws(piece)
    if (grepl("^[0-9]+-[0-9]+$", piece)) {
      ends <- as.integer(strsplit(piece, "-", fixed = TRUE)[[1]])
      out <- c(out, seq.int(ends[[1]], ends[[2]]))
    } else if (grepl("^[0-9]+$", piece)) {
      out <- c(out, as.integer(piece))
    }
  }
  sort(unique(out))
}

mode <- get_opt("mode", "both")
if (!mode %in% c("oracle", "misspec", "both")) {
  stop("Unknown mode: ", mode, ". Use oracle, misspec, or both.")
}
seeds <- parse_seed_spec(get_opt("seeds", ""))

message("Chunk mode: ", mode)
message("Chunk seeds: ", paste(seeds, collapse = ", "))
message("Output base: ", OUT_BASE)
message("INLA temp dir: ", Sys.getenv("INLA_TMPDIR", unset = Sys.getenv("TMPDIR", unset = "")))

run_seed_safe <- function(seed, dgp, information_set, raw_dir) {
  message(sprintf("[%s] seed %d started at %s", dgp, seed, format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  run_single_seed_replication(
    seed = seed,
    dgp = dgp,
    scens = CANONICAL_SCENS,
    force_rerun = FALSE,
    information_set = information_set,
    raw_dir = raw_dir
  )
  message(sprintf("[%s] seed %d finished at %s", dgp, seed, format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  gc()
  invisible(TRUE)
}

if (mode %in% c("oracle", "both")) {
  missing <- intersect(seeds, setdiff(CANONICAL_SEEDS, available_result_seeds(OUT_RAW_ORACLE, dgp = "spec_linear")))
  message("Oracle spec_linear missing in this chunk: ", paste(missing, collapse = ", "))
  for (seed in missing) {
    run_seed_safe(seed, dgp = "spec_linear", information_set = "oracle", raw_dir = OUT_RAW_ORACLE)
  }
}

if (mode %in% c("misspec", "both")) {
  missing <- intersect(seeds, setdiff(CANONICAL_SEEDS, available_result_seeds(OUT_RAW, dgp = "misspec_tanh")))
  message("Realistic misspec_tanh missing in this chunk: ", paste(missing, collapse = ", "))
  for (seed in missing) {
    run_seed_safe(seed, dgp = "misspec_tanh", information_set = "realistic", raw_dir = OUT_RAW)
  }
}

message("Chunk completed at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
