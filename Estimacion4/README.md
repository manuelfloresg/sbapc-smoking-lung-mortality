# SBAPC Simulation & Estimation Replication Package

This repository contains the complete replication code, simulation scripts, and output generation routines for the **SBAPC (Smoking-adjusted Bayesian Age-Period-Cohort) projection model** manuscript. 

This package is fully configured to reproduce the quantitative results, LaTeX tables, and high-resolution figures presented in **Section 4 (Simulation Validation)** and **Appendix C (Supplemental Sensitivity & Seed Case Studies)**.

---

## 1. Overview of the SBAPC Model

The **SBAPC** model enhances the traditional Bayesian Age-Period-Cohort (BAPC) projection framework by integrating smoking prevalence dynamics. It maps historical and future smoking trends to lung cancer and other tobacco-related cancer incidence, accounting for:
* Stratified Risk Ratios (RRs) by sex and cancer site.
* Smoking cessation risk-reversal schedules (up to a 50-year horizon).
* Stochastic cohort exposure stock calculation (the "effective exposure slide" $q_{eff}$).

By aligning the Data Generating Process (DGP) with the estimator, this codebase proves that SBAPC eliminates the historical level mismatches and systematic biases of uninformed APC models.

---

## 2. Directory Structure

```
Estimacion4/
├── R/                         # Core implementation of the BAPC/SBAPC engine
│   ├── 00_defaults.R          # Configurable defaults, scenario definitions, and colors
│   ├── 01_core_helpers.R      # Risk-reversal calculations and core utility functions
│   ├── 02_stage_models.R      # INLA Model training routines for prevalence and incidence
│   ├── 03_pipeline_sex.R      # Sex-stratified modeling pipeline
│   ├── 04_pipeline_both.R     # Two-sex orchestrator and population aggregator
│   ├── 04b_rebuilder_helpers.R# Scenario reconstruction wrappers
│   ├── 05_postprocess.R       # Projection consolidation and data framing
│   ├── 07_plots_paper.R       # Standard manuscript visualization helpers
│   ├── 09_figures_maintext.R  # Specific Section 4 main text figure builders
│   └── 31_diagnostics_against_truth.R # Multi-seed truth/estimator comparison math
├── adapters/                  # Data standardization and simulation DGP builders
│   ├── build_inputs_real.R    # Real-world data ingestion adapter
│   └── build_inputs_sim.R     # Simulation Data Generating Process (DGP) engine
├── runs/                      # Script entrypoints and parallel batch runners
│   ├── _runtime_setup.R       # PATH configurations and directory auto-creation
│   ├── _source_all.R          # Helper to source all R engine submodules
│   └── replication_diagnostics.R # CENTRAL REPLICATION ENTRYPOINT
└── results/                   # Destination folder for generated tables and figures
    └── 20260515_FINAL_PROD/
        ├── section4/          # Main manuscript Section 4 replication outputs
        └── appendixC/         # Supplemental Appendix C replication outputs
```

---

## 3. Software Dependencies & Requirements

To execute the replication pipeline, you require **R (version >= 4.5.0)** and the following libraries:

### CRAN Packages
* `dplyr`, `tidyr`, `readr` (data manipulation)
* `ggplot2` (visualization)
* `patchwork` (multi-panel figure layout)
* `future`, `future.apply` (robust multi-core parallelization)
* `stringr`, `stringi` (text normalizations)

### Bioconductor / External Packages
* **`INLA` (Integrated Nested Laplace Approximations):** Used for Bayesian APC inference.
  To install INLA, run:
  ```R
  install.packages("INLA", repos = c(gstat = "https://inla.r-inla-download.org/R/stable", CRAN = "https://cloud.r-project.org"), dep = TRUE)
  ```

---

## 4. How to Run the Replication

The entire simulation replication and output generation is automated.

### 4.1 Running the Full Pipeline (50 Seeds)
To run the full 50-seed simulation batch in parallel (highly recommended if you have 6+ cores available) and regenerate all paper outputs, execute:

```bash
cd Estimacion4
Rscript runs/replication_diagnostics.R
```

Inside R, you can trigger the orchestrator manually:
```R
source("runs/replication_diagnostics.R")
replicate_all_simulations()
```

> [!NOTE]
> The simulations utilize 6 parallel workers via `future` to distribute memory footprint safely, as INLA is memory-intensive. Running all 50 seeds takes approximately 20–30 minutes on a modern multi-core workstation.

---

## 5. Description of Generated Replications

Upon completion, all outputs are saved to the `results/20260515_FINAL_PROD/` directory:

The canonical simulation scenarios are `freeze`, `up1pc`, `down1pc`, and `quit`. The moderate increase/decrease rates are controlled in `R/00_defaults.R` via `PREV_ANNUAL_RATE_UP` and `PREV_ANNUAL_RATE_DOWN`.

Figures are exported as SVG by default through `BAPC_FIG_FORMAT <- "svg"` in `R/00_defaults.R`. Set `options(BAPC_FIG_FORMAT = "png")` or environment variable `BAPC_FIG_FORMAT=png` to export PNG instead; use `"both"` to write both formats.

Large intermediate RDS files are not intended for journal upload. Use `compact_replication_rds()` after a production run to write a lighter `raw_data_compact/` cache while preserving the original heavy `raw_data/` directory for local forensic debugging.

### 5.1 Section 4: Main Paper Outputs
* **`tab_bias_summary.tex`**: LaTeX table summarizing historical and projection bias across all scenarios. Proof that historical bias is $\approx 0\%$.
* **`fig_scenario_atlas_seed4_M.png` / `fig_scenario_atlas_seed4_F.png`**: Scenario Atlas showing true vs. estimated deaths for Males and Females across all 4 scenarios on a shared axis.
* **`fig_waterfall_seed4.png`**: The four-stage transmission waterfall showing how prevalence changes propagate to deaths.
* **`fig_reliability_calibration.png`**: Projection reliability errors across a 50-year horizon, divided into *Credible*, *Caution*, and *Risky* zones.
* **`fig_sensitivity_seed4.png`**: Scenario sensitivity chart displaying total aggregated projected deaths.

### 5.2 Appendix C: Supplemental & Distributional Analysis
* **`fig_bias_distributions.png`**: Global boxplot deconstructing the projection bias distributions across all 50 stochastic seeds.
* **`fig_case_study_best_s18.png` / `fig_case_study_best_s29.png`**: Performance atlas for the "Best Case" simulation seeds.
* **`fig_case_study_median_s21.png` / `fig_case_study_median_s44.png`**: Performance atlas for the "Median Case" simulation seeds.
* **`fig_case_study_worst_s22.png` / `fig_case_study_worst_s33.png`**: Performance atlas for the "Worst Case" simulation seeds.
* **`full_simulation_matrix.csv`**: Detailed tabular dataset containing the individual seed metrics, projections, and raw bias calculations.
