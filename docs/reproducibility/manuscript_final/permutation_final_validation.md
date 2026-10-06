# Permutation final validation

- Active source: `scripts/portable_full_pipeline/bridge_to_layer2/publication_extra/01_label_permutation_null.R`.
- Input: `data/derived/figure_inputs/step09_split_sample_level.tsv`, restricted by the final role/evaluation state in `step09_split_eval.tsv` and `manuscript_dataset_role_table.tsv`.
- Final universe: 68 contrasts (40 locked anchor, 14 primary acute validation, 11 extended validation, 3 secondary support); no human HepG2 contrasts.
- Permutations: 1,000 within each split contrast; frozen IMRSz values retained and delivery/control labels permuted by random subsets of the observed delivery size.
- Interval: empirical 2.5th and 97.5th percentiles using R quantile type 7; outside means observed < q025 or observed > q975 (strict inequalities).
- Two-sided empirical p-value: `(1 + sum(abs(null_delta) >= abs(observed_delta))) / (B + 1)`.
- Multiple testing: Benjamini-Hochberg adjustment across exactly the 68 final two-sided empirical p-values.
- Validated result: 49/68 outside the 95% null interval; 43/68 BH-adjusted two-sided p < 0.05.
- Exact stored-boundary ties: 5; strict interval classification leaves those ties inside the interval.
- Discrepancy: the stale 48/42 state was a filtered 70-contrast run. The two invalid HepG2 contrasts had already consumed 2,000 random draws, so post hoc row removal did not reproduce a final-68 run; BH was also calculated in the stale 70-test universe.
- Classification changes are recorded in `permutation_discrepancy_resolution.tsv`: one GSE264344 interval-status change and three GSE279372 BH-status changes (net +1 BH-significant contrast).
- Manuscript correction: not needed; the validated final values are 49/68 and 43/68.
