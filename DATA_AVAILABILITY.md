# Data Availability and Reproducibility Notes

## Included Data

The repository is structured to include analysis-ready, non-identifiable
datasets under:

```text
data/analysis_ready/
```

The required empirical inputs are documented in:

```text
data/metadata/uruguay_inputs_required.csv
```

If analysis-ready inputs are not included in the public repository, the
simulation study remains fully reproducible and the empirical Uruguay
application can be reproduced by placing the documented inputs in
`data/analysis_ready/` or by setting the `BAPC_PATH_*` environment variables.

## Non-Redistributed Raw Source Files

Raw source files obtained directly from Uruguayan institutions are not
redistributed unless explicit permission is available. This includes, as
applicable, raw tabulations or files received from:

- the Ministry of Public Health;
- the National Cancer Registry;
- other institutional providers of cancer incidence, mortality, population, or
  survey microdata.

## Fixed External Inputs

The empirical pipeline uses fixed external transmission inputs, including
incidence relative risks, risk-reversal schedules, and post-diagnosis mortality
probabilities. These are documented in the code and in the appendix outputs.
The external-input sensitivity analysis is deterministic and should not be
interpreted as posterior uncertainty.

## Uncertainty Labels

Simulation figures report empirical percentile ranges across simulation
replications. Uruguay mortality figures report approximate 95 percent credible
intervals for expected annual deaths derived from INLA marginal summaries.
External-input sensitivity figures report deterministic sensitivity envelopes.
