# Replication Manifest

## Code Required for Reproduction

Core code:

- `Estimacion4/R/00_defaults.R`
- `Estimacion4/R/01_core_helpers.R`
- `Estimacion4/R/01b_engine_consistency_helpers.R`
- `Estimacion4/R/01c_prediction_rebuild_helpers.R`
- `Estimacion4/R/02_stage_models.R`
- `Estimacion4/R/03_pipeline_sex.R`
- `Estimacion4/R/04_pipeline_both.R`
- `Estimacion4/R/04b_rebuilder_helpers.R`
- `Estimacion4/R/05_postprocess.R`
- `Estimacion4/R/06_qc.R`
- `Estimacion4/R/09_figures_maintext.R`
- `Estimacion4/R/10_diagnostics_methodpaper.R`
- `Estimacion4/R/31_diagnostics_against_truth.R`

Input adapters:

- `Estimacion4/adapters/build_inputs_sim.R`
- `Estimacion4/adapters/build_inputs_real.R`

Main run scripts:

- `Estimacion4/runs/run_final_simulation_200.R`
- `Estimacion4/runs/replication_diagnostics.R`
- `Estimacion4/runs/uruguay_products.R`
- `Estimacion4/runs/run_real_lung.R`
- `Estimacion4/runs/run_real_9sites.R`

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

- `Estimacion4/results/*/raw_data*/`
- `Estimacion4/runtime/`
- `.Rlib/`
- INLA temporary directories
- `*.rds`, `*.RDS`, `*.rda`, `*.RData`
- private raw source data
- Dropbox-specific source files
- credentials, tokens, API keys, and local environment files

## Manual Decisions Before Public Release

- Confirm whether analysis-ready Uruguay inputs can be redistributed.
- Confirm the final code/data license choice.
- Add the final GitHub repository URL to `CITATION.cff`.
- Add the final Zenodo DOI after creating the release.
