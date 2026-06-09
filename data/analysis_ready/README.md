# Analysis-Ready Data

Place public, non-identifiable analysis-ready empirical inputs in this folder.

Required filenames for the default Uruguay empirical pipeline:

- `uruguay_mortality_smooth_cancer.csv`
- `uruguay_population_1950_2070.dta`
- `uruguay_smoking_prevalence_harmonized.dta`
- `uruguay_incidence_smooth_1998_2022.csv`

The smoking prevalence file may also be an aggregated analysis-ready file with
columns:

```text
age, period, cohort, sex, inst, y_eff, neff
```

If the files cannot be redistributed, leave this folder empty except for this
README and set the `BAPC_PATH_*` environment variables to local copies.
