# Codex Porting Audit Report

## Audit Scope

This report records both the discarded initial implementation and the repaired
current pipeline. The files under `scripts/portable_full_pipeline/original_ported/`
were compared with the source scripts under:

```text
D:/IMRS_Project/Hypergator_scripts/InnateImmuneResponseScore/R/
D:/IMRS_Project/Hypergator_scripts/InnateImmuneResponseScore/R/modeling/
```

The branch checkout contains the historical repository copies under
`scripts/full_pipeline/provenance_legacy/`. Their local default paths are
sanitized in this update, while the executable portable route uses the
copied files under `scripts/portable_full_pipeline/original_ported/`.

## Bottom Line

The first attempt recreated core computational stages from scratch. Those
files were deleted and are disclosed below as
`REIMPLEMENTED_FROM_SCRATCH_REQUIRES_REVIEW`.

The current executable core now calls copied original scripts. The repairs
restored the original minimum-anchor behavior in Steps 06A and 06B and
removed dataset-selection policy from Step08. Dataset staging and invalid
single-anchor test-mode rejection now occur in the runner layer.

All nine current core stages are classification `A`: direct copies with
path/config/call-interface edits only. Step 06A additionally persists its
already-computed `support_by_dataset` object for packaging; this adds an
output file but does not change selection logic or any calculation.

The subsequently located Step09-to-Layer2 handoff is distributed rather than
a single core stage. Path-parameterized copies of its audit and
publication-extra scripts are now recorded in
`docs/original_script_porting_inventory.tsv` and orchestrated by
`scripts/portable_full_pipeline/run_step09_to_layer2_inputs.R`. These
components retain their original computations and write into a generated
Layer 2 staging package; inputs for which no directly ported generator exists
are explicitly labelled released references.

## 1. Files Created From Scratch

Current non-core files created from scratch:

| File | Role |
| --- | --- |
| `scripts/portable_full_pipeline/lib/portable_pipeline_utils.R` | Config/path/preflight utilities. |
| `scripts/portable_full_pipeline/00_preflight_inputs.R` | Preflight validation and reports. |
| `scripts/portable_full_pipeline/compare_generated_frozen_weights.R` | Read-only regenerated-versus-released QA comparison. |
| `scripts/portable_full_pipeline/package_generated_outputs_for_v6_figures.R` | Output staging/packaging adapter. |
| `scripts/portable_full_pipeline/run_count_to_v6_outputs.R` | Runner, logging, manifests, input staging, and stop policy. |
| `config/full_pipeline_config.yml` | Public-safe configuration. |
| `config/local_d_drive_config_TEMPLATE.yml` | Local example configuration. |
| `docs/count_level_full_pipeline.md` | Workflow documentation. |
| `docs/derived_input_lineage.tsv` | Derived-input lineage. |
| `docs/original_script_porting_inventory.tsv` | Original-to-port inventory. |
| `docs/codex_porting_audit_report.md` | This audit. |
| `docs/codex_generated_vs_ported_inventory.tsv` | File-level classification. |

Earlier core scripts created from scratch and then removed:

```text
scripts/portable_full_pipeline/01_prepare_scoring_and_split_designs.R
scripts/portable_full_pipeline/02_run_deseq2_contrasts.R
scripts/portable_full_pipeline/03_reconstruct_frozen_gene_weights.R
scripts/portable_full_pipeline/04_score_samples_with_frozen_weights.R
scripts/portable_full_pipeline/05_run_step09_evaluation.R
```

Each removed core file is classified as:

```text
REIMPLEMENTED_FROM_SCRATCH_REQUIRES_REVIEW
```

## 2. Files Copied From Original Or Provenance Scripts

| Current portable file | Original source |
| --- | --- |
| `original_ported/01_designs/02_adding_control_to_design.R` | `R/02_adding_control_to_design.R` |
| `original_ported/01_designs/03_Final_design_file_autogroup.R` | `R/03_Final_design_file_autogroup.R` |
| `original_ported/02_deseq2_contrasts/05_DEseq_contrast_delivery_control.R` | `R/modeling/05_DEseq_contrast_delivery_control.R` |
| `original_ported/03_modeling/06A_core_gene_set.R` | `R/modeling/06A_core_gene_set.R` |
| `original_ported/03_modeling/06B_gene_heterogeneity.R` | `R/modeling/06B_gene_heterogeneity.R` |
| `original_ported/03_modeling/06C_Power_analysis.R` | `R/modeling/06C_Power_analysis.R` |
| `original_ported/03_modeling/07_weight_estimation.R` | `R/modeling/07_weight_estimation.R` |
| `original_ported/04_scoring/08_score_samples.R` | `R/modeling/08_score_samples.R` |
| `original_ported/05_evaluation/09_calibration_evaluation.R` | `R/modeling/09_calibration_evaluation.R` |

The portable runner calls these copied-original files, not the existing
in-place repository provenance copies.

## 3. Files Copied And Path-Parameterized

All nine `original_ported/` scripts have portable input/output root handling.
Their original calculation blocks, filters, thresholds, formulas, and output
schemas remain authoritative.

Step 06A has one output-only porting note:

```r
# PORTING NOTE: Persisting already-computed support_by_dataset object; no selection logic changed.
```

This table was computed by the original algorithm before the write was added.

## 4. Files With Newly Written Algorithmic Logic

No current `original_ported/` core script contains newly written
modeling/scoring/evaluation logic.

`compare_generated_frozen_weights.R` is newly written and calculates overlap,
correlation, and absolute differences for QA reporting. It neither estimates
weights nor changes generated results.

The removed scratch stage files contained newly written core logic and remain
disclosed as `REIMPLEMENTED_FROM_SCRATCH_REQUIRES_REVIEW`; they are not in
the execution path.

## 5. Wrapper, Parser, Preflight, Logger, Or Packaging Files

| File | Classification |
| --- | --- |
| `run_count_to_v6_outputs.R` | Wrapper/runner; stages configured Step08 inputs and enforces test-mode stop behavior. |
| `lib/portable_pipeline_utils.R` | Config/path/preflight utility. |
| `00_preflight_inputs.R` | Preflight/check/report component. |
| `package_generated_outputs_for_v6_figures.R` | Packaging adapter only. |
| `compare_generated_frozen_weights.R` | QA comparison adapter; calculation does not feed the model. |

## 6. Was Original Logic Preserved?

Yes for the current core route, with the output-persistence qualification for
Step 06A described above. The port does not alter DESeq2 contrasts, retained
gene selection, heterogeneity, power, frozen weight estimation, Step08
scoring, or Step09 evaluation logic.

Test mode no longer changes Steps 06A or 06B. If it has fewer than the locked
anchor set and attempts weight reconstruction, the runner stops before
Step06A rather than making the original algorithm operate outside its
supported inputs.

## 7. Core Stage Classification

Classification `A` means direct port with path/config/call-interface edits
only and no algorithmic logic modification.

| Core stage | Original script | Current classification | Finding |
| --- | --- | --- | --- |
| Design/scoring metadata | `02_adding_control_to_design.R` | A | Paths parameterized. |
| Split designs | `03_Final_design_file_autogroup.R` | A | Paths parameterized. |
| DESeq2 contrasts | `05_DEseq_contrast_delivery_control.R` | A | Paths parameterized; original DESeq2 calculations retained. |
| Core/support construction | `06A_core_gene_set.R` | A | Original guard and selection restored; computed support table additionally persisted. |
| Gene heterogeneity | `06B_gene_heterogeneity.R` | A | Original minimum-anchor stop restored. |
| Gene power | `06C_Power_analysis.R` | A | Paths parameterized. |
| Frozen weight estimation | `07_weight_estimation.R` | A | Paths parameterized. |
| Sample scoring | `08_score_samples.R` | A | Paths parameterized; dataset scoping moved to runner staging. |
| Step09 evaluation | `09_calibration_evaluation.R` | A | Paths parameterized; original delta-IMRSz/AUC logic retained. |

No current core stage is `B`, `C`, or `D`.

## 8. Logic Modification Review

| Component | Algorithmic logic modified in current original port? | Evidence/handling |
| --- | --- | --- |
| DESeq2 contrast generation | NO | Path roots only. |
| 06A core/support gene construction | NO | Original `K_present < 2` skip restored; output persistence only. |
| 06B gene heterogeneity | NO | Original `K_present < 3` stop restored. |
| 06C gene power | NO | Path root only. |
| 07 frozen weight estimation | NO | Path root only. |
| 08 sample scoring | NO | Internal configured-dataset filter removed; runner stages input roots. |
| 09 calibration/evaluation | NO | Path roots only. |

## 9. Recreated Stages And Why

The initial reimplementation was an error in approach, not a response to a
missing original dependency. Original source scripts were available. The
scratch scripts were removed and replaced with copied originals before the
current repair was finalized.

## 10. Missing Dependencies Or Assumptions

| Item | Finding |
| --- | --- |
| Original helper files | No separately sourced helper required by the nine core scripts was found missing. |
| R packages | The original route requires its declared packages, including `DESeq2`, tidyverse components, `ggplot2`, and `pROC`. |
| Metadata scope | The count-level pipeline assumes verified metadata exists; upstream raw metadata creation is outside this workflow. |
| Robustness tables | Steps 02-09 do not generate these tables directly; the ported distributed Layer 2 bridge now regenerates available audit/publication-extra tables after `all_scored`. Active five-anchor leave-one-anchor-out figure inputs without a matching ported generator remain labelled released references. |
| Step06A locked candidate list | The original 06A source includes `GSE262515`; this is preserved because changing it would change the original script. Runs should use a clean generated output root to avoid stale additional DE files. |
| Single-anchor test mode | The original workflow cannot reconstruct frozen weights from one anchor. The runner now refuses that reconstruction. |

## Verification After Repair

Verification completed on May 25, 2026:

| Check | Result |
| --- | --- |
| Parse runner and all nine `original_ported` R scripts | PASS |
| Local-template dry run | PASS; all configured count and metadata inputs found; only preflight executed. |
| Full `all_scored` count-level execution | PASS; weights used the five locked anchors only, all available configured dataset arms were scored/evaluated, and canonical frozen weights matched within tolerance. |
| Ported Step09-to-Layer2 bridge execution | PASS; generated input package and manifest written below configured generated output root. |
| Layer 2 generated-package dry run | PASS; required manuscript figure, enrichment, and supplementary table steps reported ready. |
| Forbidden hard-coded path scan of portable scripts/public default config | PASS; no `D:/IMRS_Project`, `D:\IMRS_Project`, `setwd(`, `/blue/`, or `C:\` occurrences. Local-template and porting-inventory source-location documentation are intentionally excluded. |
| Step06A/06B bypass search | PASS; no test-mode continuation branches remain. |
| Step08 scope-filter search | PASS; no `IMRS_PORTED_DATASETS` filter remains inside the original-port scoring script. |

## Required Repair Plan

The required repairs have been applied:

| Previously problematic stage | Required repair | Status |
| --- | --- | --- |
| Step06A test-mode continuation | Restore original skip behavior and move policy outside original code. | Completed. |
| Step06B test-mode continuation | Restore original minimum-anchor stop. | Completed. |
| Step08 configured-dataset filter | Remove from original port; stage selected inputs in wrapper. | Completed. |
| Single-anchor test reconstruction | Stop in runner before modeling. | Completed. |

The removed scratch implementations remain documented historical artifacts
only. No replacement core computation remains in the current executable path.
