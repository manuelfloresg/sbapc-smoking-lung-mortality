# Runner de Shard para producción SBAPC
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Uso: Rscript run_shard.R <shard_id> <batch_id>")

shard_id <- as.integer(args[1])
batch_id <- args[2]

cat(">>> Cargando fuentes...\n")
source("runs/_source_all.R")
cat(">>> Fuentes cargadas.\n")

# Cargar la grilla original
SEEDS      <- 1:50
DGPS       <- c("spec_linear", "misspec_tanh")
SCENARIOS  <- c("freeze", "quit", "up1pc", "down1pc")
task_grid <- expand.grid(seed = SEEDS, dgp = DGPS, scn = SCENARIOS, stringsAsFactors = FALSE)

# Dividir tareas (6 shards)
num_shards <- 6
total_tasks <- nrow(task_grid)
shard_indices <- split(1:total_tasks, cut(1:total_tasks, num_shards, labels = FALSE))
my_indices <- shard_indices[[shard_id]]

base_out <- file.path(BAPC_PATHS$results, "batches", batch_id)
plot_out <- file.path(base_out, "audit_plots")
dir.create(plot_out, recursive = TRUE, showWarnings = FALSE)

# Configurar INLA aislado para este proceso
worker_tmp <- file.path("C:/tmp_inla", paste0("shard_", shard_id))
dir.create(worker_tmp, recursive = TRUE, showWarnings = FALSE)
options(INLA.tmpdir = worker_tmp)
INLA::inla.setOption(working.directory = worker_tmp)

message(sprintf(">>> SHARD %d iniciado (%d tareas)", shard_id, length(my_indices)))

for (i in my_indices) {
  task <- task_grid[i, ]
  res_path <- file.path(base_out, sprintf("res_%s_s%d_%s.rds", task$dgp, task$seed, task$scn))
  
  if (file.exists(res_path)) {
    message(sprintf("Skipping task %d (already exists)", i))
    next
  }
  
  tryCatch({
    message(sprintf("[%s] Shard %d: Procesando tarea %d (Seed %d, DGP %s, Scn %s)...", 
                    Sys.time(), shard_id, i, task$seed, task$dgp, task$scn))
    
    sim <- simulate_PIM_data(cause_id = task$dgp, seed = task$seed, dgp = task$dgp, scenario_name = task$scn)
    inputs <- build_inputs_sim(sim, cause_id = task$dgp)
    prev_cfg <- make_prev_config(scenario = task$scn)
    
    res <- run_pipeline_both(
      mort_hist_tbl = inputs$mort_hist, pop_all_tbl = inputs$pop_all,
      inc_hist_tbl = inputs$inc_hist, prev_micro_df = inputs$prev_data,
      cause_id_override = task$dgp, beta_mode = "fixed_rr_offset", prev_cfg = prev_cfg,
      emit_prev_diag_write = FALSE
    )
    
    if (!is.null(res)) {
      saveRDS(res, res_path)
      # Gráfico rápido
      prefix_val <- sprintf("%s_s%d_%s", task$dgp, task$seed, task$scn)
      try(compare_pipeline_to_truth(res, sim, out_dir = plot_out, prefix = prefix_val), silent = TRUE)
    }
  }, error = function(e) {
    message(sprintf("!!! SHARD %d ERROR en tarea %d: %s", shard_id, i, e$message))
  })
}

message(sprintf(">>> SHARD %d FINALIZADO", shard_id))
