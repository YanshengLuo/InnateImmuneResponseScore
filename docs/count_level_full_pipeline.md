# Prepared Count-Level Reconstruction Workflow

## Scope

This workflow ports the original IMRS count-level scripts into
`scripts/portable_full_pipeline/original_ported/` and invokes them in their
historical order through a config-driven runner. The original design,
DESeq2, retained-gene, heterogeneity, power, frozen-weight, scoring, and
Step09 evaluation calculations remain in those copied scripts.

The workflow starts from prepared featureCounts-style gene-count matrices and
verified metadata. It does not retrieve FASTQ or BAM files, align reads, or
run feature counting.

The quick reviewer-facing runner, `Rscript run_all_manuscript_outputs_v6.R`,
starts from released derived inputs. This count-level runner starts earlier,
reconstructs derived inputs into a staging area, and may then invoke the
existing v6 manuscript runner against those regenerated inputs.

## Validated Canonical Scope

On May 25, 2026, the canonical Layer 1 count-level workflow completed for the
five locked anchors configured in the local-example profile. It reran the
ported original design/split preparation, DESeq2 anchor contrasts, Steps
06A-07 weight reconstruction, frozen-weight comparison, Step08 scoring,
Step09 anchor-level evaluation, and figure-input packaging. The comparison
reported:

```text
PASS: regenerated frozen weights match the released canonical table within tolerance.
```

Canonical mode validates locked-anchor frozen-model construction. It does not
generate new count-derived score tables for non-anchor manuscript validation
datasets. Those full manuscript scoring/evaluation inputs remain released
Layer 2 inputs for figure and table reproduction.

`all_scored` mode has a different application scope: it reconstructs and
checks the frozen model from the same five locked anchors only, then applies
the fixed weights to every configured manuscript dataset with available
count/design inputs. It is the prerequisite for regenerating the distributed
Step09-to-Layer2 bridge tables.

## Ported Original Order

| Stage | Ported original script | Purpose |
| --- | --- | --- |
| 01 | `original_ported/01_designs/02_adding_control_to_design.R` | Build scoring designs with `CONTROL`/`DELIVERY`. |
| 01 | `original_ported/01_designs/03_Final_design_file_autogroup.R` | Build delivery-versus-control split designs. |
| 02 | `original_ported/02_deseq2_contrasts/05_DEseq_contrast_delivery_control.R` | Produce raw-count DESeq2 anchor/calibration contrasts. |
| 03 | `original_ported/03_modeling/06A_core_gene_set.R` | Construct anchor retained-gene support. |
| 03 | `original_ported/03_modeling/06B_gene_heterogeneity.R` | Estimate cross-anchor heterogeneity. |
| 03 | `original_ported/03_modeling/06C_Power_analysis.R` | Estimate gene power/information. |
| 03 | `original_ported/03_modeling/07_weight_estimation.R` | Estimate frozen anchor weights. |
| 04 | `original_ported/04_scoring/08_score_samples.R` | Apply frozen coefficients to sample counts. |
| 05 | `original_ported/05_evaluation/09_calibration_evaluation.R` | Evaluate delivery-minus-control IMRSz by split. |

`package_generated_outputs_for_v6_figures.R` is a thin filename/staging
adapter, not an analytical step. `compare_generated_frozen_weights.R` is a
read-only canonical comparison adapter. The full original-to-port mapping is
recorded in `docs/original_script_porting_inventory.tsv`.

## Required Inputs

The public-safe default layout is:

```text
data/counts/<DATASET>/featurecounts/validation/gene_counts_clean.tsv
data/curated_metadata/verified_metadata/<DATASET>_design.tsv
```

Each count matrix must have a gene identifier column followed by sample
columns. Counts must be finite, nonnegative integers representing prepared raw
feature counts; TPM, CPM, or transformed values are rejected by preflight and
by the original DESeq2/scoring scripts.

Verified metadata must identify samples and provide enough curated treatment
information for the original scripts to produce exact `CONTROL` and
`DELIVERY` labels. Sample identifiers must overlap count-matrix column names.

## Locked Anchors And Modes

Canonical reconstruction requires all five locked anchors:

```text
GSE167521
GSE264344
GSE279372
GSE279744
GSE39129
```

These five acute mouse datasets are the production retained-gene and frozen-
coefficient construction set. `GSE262515` is validation/secondary support only
and does not contribute to Steps 06A-07.

The strict-3 set (`GSE39129`, `GSE167521`, and `GSE264344`) is used only for
supporting sensitivity or ablation analyses. In particular, threshold-
sensitivity analyses that use strict-3 inputs are robustness checks; they do
not define, replace, or refit the production five-anchor frozen model.

In `canonical` mode, missing anchor counts or metadata cause preflight to
fail. Generated `gene_weights.tsv` is compared with
`data/derived/frozen_gene_weights.tsv`, when present, and differences are
reported under the configured output root. The released table is never
overwritten.

In `all_scored` mode, the five locked anchors remain mandatory and are the
only inputs to retained-gene construction and weight estimation. Step08 and
Step09 additionally run for available configured validation datasets using
the already frozen weights. The wrapper stages `GSE262515_cell_line` and
`GSE262515_tissue` from the curated arms of the prepared `GSE262515` input;
the copied original scoring/evaluation code is unchanged.

In `test` mode, a subset such as `GSE39129` may be inspected through
preflight, design generation, and optionally DESeq2. Outputs are placed below
`<output_root>/test_only/` and marked `TEST_ONLY_NON_CANONICAL`; no
canonical-equivalence claim is made. Because the historical 06A and 06B
scripts require multiple anchors, the runner stops before model reconstruction
when the locked anchors are incomplete. A scoring-only test may be run from
Stage 04 only after an existing frozen `gene_weights.tsv` has been explicitly
placed in the test output anchor directory and the deliberate reuse is
acknowledged with `--force --start-stage 04`; Step08 never refits weights.

## Scoring And Evaluation

The ported `08_score_samples.R` remains authoritative for scoring. The runner
stages configured count and scoring-design inputs before it is invoked, so
dataset-selection policy is not embedded in the original algorithm file. It loads
the frozen anchor coefficient table once, uses raw integer counts, estimates
DESeq2 size factors, applies `log2(normalized_count + 1)`, forms control-only
gene z-scores with its fixed SD floor, computes weighted sample scores, and
control-standardizes them to IMRSz. It checks control count and gene coverage
and writes sample score, QC, and top-contributor files. It does not refit
weights or modify retained genes using validation labels.

When the coefficient table contains `beta_meta`, Step08 applies `beta_meta`.
It copies the penalized `weight` field into `beta_meta` only when `beta_meta`
is absent. The released canonical table contains `beta_meta`, so released
scores use `beta_meta`; `weight` remains an audit/fallback field.

The ported `09_calibration_evaluation.R` is evaluation only. It computes
delivery-minus-control IMRSz differences, group counts and means,
directionality, effect summaries, and secondary AUC where its original
criteria permit. It does not alter the scoring model.

In the validated canonical run, Step09 evaluated configured anchor score
files. Non-anchor validation datasets without newly generated Step08 score
files were expectedly skipped because they are outside canonical
locked-anchor reconstruction scope.

## All-Scored Mode And Generated Layer2 Handoff

There is no single historical bridge script. The portable bridge runner
invokes path-parameterized copies of the original scripts from
`audit/scripts/` and `05_score/publication_extra_analyses/scripts/` through
the original distributed sequence. It must be run after `all_scored` Step09
so manuscript roles, context audits, robustness summaries, and comparator
tables see generated validation evaluations.

Bridge outputs are written only beneath:

```text
<output_root>/layer2_generated_inputs/data/derived/
<output_root>/layer2_generated_inputs/data/derived/figure_inputs/
<output_root>/layer2_generated_inputs/layer2_bridge_manifest.tsv
```

The bridge manifest marks any active Layer 2 input copied from released
reference data because no directly ported generator was found, including the
active five-anchor leave-one-anchor-out figure inputs. It does not present
released references as regenerated bridge results.

In short, `canonical` validates frozen-weight reconstruction, while
`all_scored` applies those frozen weights to all available configured
manuscript datasets. After `all_scored`, run
`run_step09_to_layer2_inputs.R` to construct the generated Layer 2 input
package, then point `run_all_manuscript_outputs_v6.R` to that generated
configuration when a count-derived manuscript-output run is desired.

## Outputs

All regenerated products are written below the configured `paths.output_root`,
including:

```text
00_preflight/
01_designs/scoring/
01_designs/splited/
04_de/
05_score/anchors/gene_weights.tsv
05_score/transfer/scores/
05_score/transfer/qc/
05_score/transfer/eval/
06_derived_for_figures/data/derived/figure_inputs/
layer2_generated_inputs/data/derived/figure_inputs/
logs/
manifests/
```

Packaging creates staged v6 figure inputs below
`<output_root>/06_derived_for_figures/data/derived/`. It does not write into
released `data/derived/` or `data/derived/figure_inputs/`. When manuscript
execution is enabled in a local configuration, the staged config routes the
existing v6 plotting/table runner to this generated derived-input root.

## Commands

Run from the repository root:

```sh
Rscript scripts/portable_full_pipeline/run_count_to_v6_outputs.R --config config/full_pipeline_config.yml --dry-run
Rscript scripts/portable_full_pipeline/run_count_to_v6_outputs.R --config config/full_pipeline_config.yml --mode canonical
Rscript scripts/portable_full_pipeline/run_count_to_v6_outputs.R --config config/local_d_drive_config_TEMPLATE.yml --mode canonical
Rscript scripts/portable_full_pipeline/run_count_to_v6_outputs.R --config config/local_d_drive_config_TEMPLATE.yml --mode all_scored
Rscript scripts/portable_full_pipeline/run_step09_to_layer2_inputs.R --config config/local_d_drive_config_TEMPLATE.yml --mode all_scored
Rscript scripts/portable_full_pipeline/run_count_to_v6_outputs.R --config config/local_d_drive_config_TEMPLATE.yml --mode test
```

After an all-scored bridge run, Layer 2 can be checked against the generated
input package without changing its default released-input configuration:

```sh
Rscript run_all_manuscript_outputs_v6.R --config <output_root>/layer2_generated_inputs/layer2_generated_inputs_config.yml
```

Optional controls are `--start-stage 04`, `--stop-stage 05`, `--force`, and
`--skip-manuscript`. Dry run performs preflight, lists configured and
discovered inputs, prints the ported original execution order and planned
outputs, and does not run DESeq2, modeling, scoring, Step09, or packaging.

## Troubleshooting

| Problem | Resolution |
| --- | --- |
| Sample mismatch | Inspect `00_preflight/sample_match_report.tsv`; count column names must match verified metadata sample IDs. |
| Missing metadata | Supply the verified `<DATASET>_design.tsv` below `paths.verified_metadata_dir` or correct the configured directory. |
| Missing anchor | Canonical mode needs all five anchors. Test mode with an incomplete anchor set may run preparatory stages, but cannot reconstruct frozen weights with the original workflow. |
| Single-anchor test modeling stop | Supply all locked anchors for reconstruction, or supply an existing frozen `gene_weights.tsv` and run scoring-only test use with `--force --start-stage 04`. |
| Non-integer counts | Provide prepared raw feature-count matrices, not TPM, CPM, or log-normalized data. |
| DESeq2 missing | Install Bioconductor `DESeq2` in the R library used for the run and rerun after preflight passes. |
| Canonical weight mismatch | Review `manifests/frozen_gene_weight_comparison.tsv`, `frozen_gene_weight_differences.tsv`, and upstream DE/design logs; mismatches are surfaced rather than overwritten. |
