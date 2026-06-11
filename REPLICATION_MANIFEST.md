# Replication Manifest

## Code Required for Reproduction

Core code:

- `code/R/00_defaults.R`
- `code/R/01_core_helpers.R`
- `code/R/01b_engine_consistency_helpers.R`
- `code/R/01c_prediction_rebuild_helpers.R`
- `code/R/02_stage_models.R`
- `code/R/03_pipeline_sex.R`
- `code/R/04_pipeline_both.R`
- `code/R/04b_rebuilder_helpers.R`
- `code/R/05_postprocess.R`
- `code/R/06_qc.R`
- `code/R/09_figures_maintext.R`
- `code/R/10_diagnostics_methodpaper.R`
- `code/R/31_diagnostics_against_truth.R`

Input adapters:

- `code/adapters/build_inputs_sim.R`
- `code/adapters/build_inputs_real.R`

Main run scripts:

- `code/runs/run_final_simulation_200.R`
- `code/runs/replication_diagnostics.R`
- `code/runs/uruguay_products.R`
- `code/runs/run_real_lung.R`
- `code/runs/run_real_9sites.R`

Public wrappers:

- `scripts/reproduce_simulations.R`
- `scripts/reproduce_uruguay.R`
- `scripts/check_inputs.R`

## Products Included in `output/`

The public `output/` folder should include final PDF/SVG figures, LaTeX table
bodies, CSV table data, figure/table notes, and float inventories for:

- Section 4;
- Appendix C;
- Section 5;
- Appendix D.

Large Monte Carlo detail files and RDS worker outputs are intentionally excluded.

## Do Not Upload

- `code/results/*/raw_data*/`
- `code/runtime/`
- `.Rlib/`
- INLA temporary directories
- `*.rds`, `*.RDS`, `*.rda`, `*.RData`
- private raw source data
- Dropbox-specific source files
- credentials, tokens, API keys, and local environment files

## Public Release Status

- Analysis-ready Uruguay inputs are included under `data/analysis_ready/`.
- Raw institutional source files are not redistributed.
- The code license and repository citation metadata are recorded in `LICENSE`
  and `CITATION.cff`.
- The public archive DOI is recorded in `README.md` and `CITATION.cff`.
