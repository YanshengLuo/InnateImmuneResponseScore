# Reviewer Reproduction Status

Generated for the IMRS repository-script package in `<repository_root>`.

## 1. Reproduction status summary

The reviewer-facing reproduction layer is functional for the staged manuscript-output package.

The controlled reviewer-facing runner successfully regenerated:

- Supplementary Figure S2 enrichment outputs
- Supplementary Tables S1-S5
- NAR G&B readiness materials

The final reproducibility equivalence check passed:

- `full equivalence_check_success`: yes
- `unexpected differences`: 0
- Remaining differences are expected metadata, path, workbook, or rendering differences, not result-value differences.

## 2. What the runner regenerates

`run_all_manuscript_outputs_v6.R` regenerates the active reviewer-facing outputs from staged derived inputs.

Executed steps:

- `supplementary_figure_s2_enrichment`
- `supplementary_tables_s1_s5`
- `nar_gb_readiness`

Generated output folders:

- `results_release_templates\Priority3_gene_program_enrichment`
- `results_release_templates\supplementary_tables`
- `results_release_templates\NAR_GB_readiness`
- `results_release_templates\logs`

## 3. What the runner intentionally does not run

The reviewer-facing runner does not run:

- raw-data retrieval
- alignment/counting
- full framework reconstruction
- frozen-gene reconstruction
- HiPerGator/SLURM scripts
- archive/legacy scripts
- `full_pipeline` scripts

These steps are documented separately and may require public data retrieval, additional configuration, or HiPerGator/HPC.

## 4. Optional skipped steps

The following steps remain optional and skipped in the current reviewer-facing run:

- `main_figure_regeneration`
- `targeted_v6_figure_regeneration`

Reason: these depend on legacy figure-generation inputs or script dependencies that are documented but are not part of the current reviewer-facing execution path.

The submitted/main figure outputs should be provided as release artifacts, while the active runner currently focuses on regenerating Supplementary Figure S2, Supplementary Tables S1-S5, and readiness materials.

## 5. Inputs required for reviewer-facing reproduction

Released derived inputs are staged under:

`data_release_templates\derived`

Key staged inputs include:

- `frozen_gene_weights.tsv`
- `gene_power.tsv`
- `gene_symbol_mapping.tsv`
- `supplement_dataset_split_provenance_v7.tsv`
- `weak_dataset_paper_context_audit.tsv`
- `label_permutation_null_summary.tsv`
- `leave_one_gene_out_summary.tsv`
- `gene_dominance_summary.tsv`
- `threshold_sensitivity_summary.tsv`
- `leave_one_anchor_out_summary.tsv`
- `baseline_signature_summary_by_group.tsv`
- `coefficient_sensitivity_summary.tsv`

Staged source notes:

- `frozen_gene_weights.tsv` was staged from `<local_project_root>\05_score\anchors\gene_weights.tsv`.
- `gene_power.tsv` was staged from `<local_project_root>\05_score\anchors\gene_power.tsv`.
- `gene_symbol_mapping.tsv` was staged from `<local_project_root>\05_score\anchors\gene_symbol_mapping.tsv`.
- `baseline_signature_summary_by_group.tsv` and `coefficient_sensitivity_summary.tsv` were staged from `<local_project_root>\05_score\publication_extra_analyses\results`.

## 6. Frozen-gene provenance

The reviewer-facing manuscript reproduction path starts from frozen derived inputs to regenerate figures and supplementary tables. The full framework reconstruction path documents how `frozen_gene_weights.tsv` can be regenerated from upstream clean-gene and locked-anchor inputs. These are intentionally separated to keep the reviewer reproduction path stable while preserving provenance for the frozen IMRS coefficient table.

The frozen-gene reconstruction script is not called by default.

## 7. Equivalence validation

The final equivalence check confirmed:

- Supplementary Figure S2 matched original v6 enrichment outputs, including term IDs, database labels, gene counts, and adjusted p-values.
- Supplementary Tables S1, S2, S3, S4, and S5 matched original v6 outputs after staging comparator benchmarking and coefficient sensitivity sources.
- S4 robustness summary now has 7 rows and includes `Comparator benchmarking` and `Coefficient sensitivity`.
- NAR G&B readiness files are present and nonzero; differences from original v6 are expected because repository-local paths, timestamps, session info, and manifest scope differ.

The final comparison covered 38 files:

- Matching or nonzero structural matches: 29
- Expected differences: 9
- Unexpected differences: 0

## 8. Safety validation

During controlled execution:

- no `hpc_hypergator` scripts were run
- no `archive_legacy` scripts were run
- no executable `full_pipeline` scripts were run
- frozen-gene reconstruction was not run
- raw-data/full-framework reconstruction was not run
- no original `revised_plots_v6` files were modified
- no IMRS scoring, validation, robustness, enrichment statistics, figure-data values, or supplementary-table source values were intentionally changed

## 9. Remaining manual actions before public release

- Replace GitHub repository URL placeholder.
- Replace reviewer-access link placeholder.
- Replace Zenodo/archive DOI placeholder.
- Choose a real open-source license instead of `LICENSE_PLACEHOLDER.txt`.
- Add `renv.lock` or `environment.yml`.
- Decide raw-data redistribution versus accession-only retrieval.
- Review optional missing `manuscript_dataset_role_table.tsv` warning if needed.
- Document optional main-figure regeneration status clearly in README/release notes.

## 10. Recommended reviewer command

From `<repository_root>`, after confirming `config\config.yml` points to the staged repository folders:

```sh
Rscript run_all_manuscript_outputs_v6.R
```

Expected required outputs are written under `results_release_templates`.

## 11. Final status statement

Overall, the repository package now provides a functional reviewer-facing reproduction layer for the active manuscript-output package, with full raw-data/framework reconstruction documented separately.

