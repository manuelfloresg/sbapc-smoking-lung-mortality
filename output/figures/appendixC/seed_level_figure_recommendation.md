# Seed-Level Figure Recommendation

Recommendation: retain no more than two seed-level figures in Appendix C. Keep them explicitly illustrative and do not use them as evidence for average performance.

## Retain

### `fig_case_study_median_s9.svg`
Shows a representative single-seed trajectory diagnostic for the quit scenario. It is useful as a concrete visual complement to aggregate recovery figures. It should remain in Appendix C only.

## Main Text, Not Duplicated In Appendix C

- `fig_transmission_map_support_compare_seed4_M.svg`: shows the prevalence-to-effective-exposure-to-incidence-to-mortality chain for one male seed, including Truth, Observed-window SBAPC, and Full-support SBAPC. It is useful enough for the main text support-window discussion and should not be duplicated in Appendix C.

## Drop Or Keep As Internal Diagnostics

- `fig_scenario_atlas_seed4_M.svg` and `fig_scenario_atlas_seed4_F.svg`: visually rich but redundant with the aggregate scenario-effect recovery figure.
- `fig_waterfall_seed4.svg`: useful for internal explanation, but the transmission-map figure is a more direct chain diagnostic.
- `fig_sensitivity_seed4.svg`: single-seed scenario sensitivity is redundant once scenario-effect recovery is aggregated across seeds.
- `fig_transmission_map_seed4_M.svg`: superseded by the support-comparison transmission map if that diagnostic is retained.
- `fig_transmission_map_support_compare_seed4_F.svg`: substantively redundant with the male pathway illustration for the current narrative; do not include it in Appendix C unless the text later makes sex-specific pathway differences central.
- `fig_case_study_best_s26.svg` and `fig_case_study_worst_s41.svg`: useful internally for stress-testing, but too anecdotal for the supplement unless the text explicitly discusses heterogeneity across seeds.
