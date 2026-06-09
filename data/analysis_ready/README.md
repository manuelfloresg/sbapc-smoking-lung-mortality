# Analysis-Ready Data

Place public, non-identifiable analysis-ready empirical inputs in this folder.

Required filenames for the default Uruguay empirical pipeline:

- `uruguay_mortality_smooth_cancer.csv`
- `uruguay_population_1950_2070.dta`
- `uruguay_smoking_prevalence_aggregated.csv`
- `uruguay_incidence_smooth_1998_2022.csv`

The smoking prevalence file is an aggregated analysis-ready file with columns:

```text
age, period, cohort, sex, inst, y_eff, neff
```

The source harmonized survey microdata are not included in this repository; the
aggregated prevalence cells are sufficient for the SBAPC empirical pipeline.
