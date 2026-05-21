# SBAPC Simulation Replication Package

This folder contains the replication code for the simulation section of the
SBAPC mortality projection paper. The current workflow is centered on the
simulation DGP, the SBAPC estimator, and diagnostics comparing projected
trajectories against simulated truth.

## Current Method Defaults

The canonical simulation scenarios are:

- `freeze`: prevalence frozen at 2022.
- `up1pc`: prevalence increases at the rate set by `PREV_ANNUAL_RATE_UP`.
- `down1pc`: prevalence decreases at the rate set by `PREV_ANNUAL_RATE_DOWN`.
- `quit`: prevalence moves to the quit/floor scenario.

The current edge-completion default for prevalence outside the observed window is:

```r
PREV_EDGE_COMPLETION_MODE <- "constant_boundary"
```

Alternative sensitivity modes remain available:

- `damped_apc`
- `apc_posterior`
- legacy `carry_states`

The default was changed after multi-seed diagnostics showed that unconstrained
APC posterior extrapolation was too aggressive at the observation-window edges.

Model labels used in paper-facing outputs are:

- `SBAPC`: the full sequential smoking-informed estimator.
- `BAPC benchmark`: an autonomous APC benchmark without the smoking-transmission channel.
- `Incidence-anchored SBAPC`: a decomposition variant used mainly in diagnostics; mortality remains anchored to projected incidence, but the prevalence-channel scenario contribution is disabled.
- `Full-support SBAPC`: an oracle-support diagnostic used to separate support-window effects from other sources of error.

## Directory Structure

```text
Estimacion4/
  adapters/
    build_inputs_real.R
    build_inputs_sim.R
  R/
    00_defaults.R
    01_core_helpers.R
    01b_engine_consistency_helpers.R
    01c_prediction_rebuild_helpers.R
    02_stage_models.R
    03_pipeline_sex.R
    04_pipeline_both.R
    04b_rebuilder_helpers.R
    05_postprocess.R
    06_qc.R
    09_figures_maintext.R
    10_diagnostics_methodpaper.R
    31_diagnostics_against_truth.R
  runs/
    _runtime_setup.R
    _source_all.R
    _source_diagnostics.R
    diagnostic_truth_comparison.R
    replication_diagnostics.R
    run_audit_simulation.R
    run_full_replication_50.R
    run_method_paper.R
    run_real_9sites.R
    run_real_lung.R
```

Generated results are written under `results/` and are intentionally ignored by
Git. Heavy RDS outputs should not be uploaded as journal replication files.

## Runtime Safety

INLA temporary files must not be created inside Dropbox. The runtime setup
redirects R and INLA temporary paths to a local temp directory, currently:

```text
C:/tmp_inla
```

Run scripts from `D:/Git/Bloomberg_2025/Estimacion4`, not from Dropbox.

## R Setup

The working environment currently uses R 4.6 with a local user library:

```text
D:/Git/Bloomberg_2025/.Rlib/4.6
```

The `.Rlib/` folder is local infrastructure and is ignored by Git.

## Main Simulation Entry Point

From `Estimacion4/`:

```r
source("runs/replication_diagnostics.R")
replicate_all_simulations(n_cores = 6, force_rerun = TRUE)
```

The default full run is configured through the `BAPC_N_SEEDS` environment
variable. If it is not set, the simulation hub uses 50 seeds:

```r
Sys.setenv(BAPC_N_SEEDS = "50")
```

For the final 200-seed simulation run used for manuscript and supplement
products:

```r
setwd("D:/Git/Bloomberg_2025/Estimacion4")
Sys.setenv(BAPC_FINAL_N_CORES = "4")  # increase to 6 only if memory is stable
source("runs/run_final_simulation_200.R")
```

This runs the well-specified observed-window design, the full-support diagnostic
design, and the misspecified transmission design, then regenerates all Section 4
and Appendix C products.

To write a production candidate to a fresh folder, set `BAPC_OUT_BASE` before
loading the replication hub:

```r
Sys.setenv(BAPC_OUT_BASE = "results/20260521_FINAL_200SEEDS")
Sys.setenv(BAPC_N_SEEDS = "200")
source("runs/replication_diagnostics.R")
replicate_final_simulations(seeds = 1:200, n_cores = 4, force_rerun = FALSE)
```

`replicate_main_paper()` also generates the seed-4 full-support/oracle run needed
for the support-window transmission map when it is not already present.

## Important Outputs

The main diagnostic products currently used for the simulation section include:

- scenario atlas by sex;
- transmission map comparing truth, window-limited SBAPC, and full-support SBAPC;
- mortality scenario-effect recovery figure;
- mortality scenario-effect recovery table;
- edge-completion sensitivity diagnostics;
- bias and reliability summaries across seeds.

Figures are exported as both SVG and PDF by default through `BAPC_FIG_FORMAT`.

## Replication Packaging Rule

For journal submission, include source code, README, and small deterministic
tables/figures needed for manuscript reproduction. Exclude:

- `results/raw_data*/`
- `runtime/`
- `.Rlib/`
- INLA temp folders
- exploratory scratch scripts and logs
