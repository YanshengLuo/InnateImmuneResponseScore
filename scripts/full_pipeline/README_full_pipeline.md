# Full Computational Provenance Pipeline

## Purpose

This folder documents Layer 1 of IMRS reproducibility: the path from public gene-count matrices and verified sample metadata to frozen gene weights, sample-level IMRS scores, and Step09 split-level evaluation tables.

Layer 1 provides provenance for the frozen, transfer-oriented IMRS scoring framework. It is not the default reviewer-facing run. The default quick reproduction path is Layer 2, implemented by `run_all_manuscript_outputs_v6.R`, which starts from released derived inputs and regenerates manuscript figures, Supplementary Tables S1-S5, and Priority3 gene-program enrichment outputs.

No script in this folder is called by the default Layer 2 runner.

## Scope

The intended Layer 1 computation is:

```text
public gene count matrices + verified metadata
  -> scoring designs and delivery-versus-control split definitions
  -> DESeq2 delivery-versus-control contrasts
  -> locked-anchor support, heterogeneity, and information-content filtering
  -> frozen gene weights
  -> sample-level frozen IMRS scores
  -> Step09 split-level delivery-minus-control Delta IMRSz evaluation tables
  -> released derived inputs consumed by Layer 2
```

Upstream retrieval, alignment, or feature counting from public FASTQ/BAM resources belongs to raw-data/HPC provenance. Account-specific HPC submission scripts are not distributed in this public release candidate; future public HPC materials should be supplied as sanitized templates. This upstream work is separate from both this count-level reconstruction description and the default reviewer run.

## Principal Inputs

Layer 1 expects the following input classes after any public-data retrieval/counting work has been completed:

| Input class | Expected project-stage location in original workflow | Release-facing representation |
| --- | --- | --- |
| Public count matrices | `03_counts/<DATASET>/featurecounts/validation/gene_counts_clean.tsv` | Obtain from public data processing or provide in an archival artifact where appropriate |
| Verified sample metadata | `00_metadata/verified_metadata/` | `data/curated_metadata/` |
| Treatment/control and split definitions | `00_metadata/verified_metadata/scoring/` and `00_metadata/verified_metadata/splited/` | `data/split_designs/` and curation documentation |
| Locked-anchor dataset definition | Embedded in Step06 scripts and documented curation tables | Documentation and released derived/audit tables |

Some scripts preserve historical local project paths as provenance copies. Their current packaging status is documented in `full_pipeline_script_inventory.tsv`; path cleanup and an explicit environment configuration are required before a portable Layer 1 execution attempt.

## Pipeline Sequence

### 1. Metadata preparation and design definition

`METADATA_BUILD.R` reads per-dataset SRA metadata tables and derives dataset-level sample/design tables with `condition_simple` assignments. `01_run_all_metadata_build.R` is a wrapper for applying that build across datasets.

`02_adding_control_to_design.R` constructs dataset-level scoring designs, including `CONTROL` and `DELIVERY` labels used by frozen scoring. `03_Final_design_file_autogroup.R` constructs DESeq2-ready delivery-versus-control split designs stratified by tissue, timepoint, batch, and delivery group, including its documented time-zero control fallback behavior.

These steps encode metadata and contrast definitions. Manual curation documentation is released separately because public datasets do not express every group, tissue/timepoint, manuscript role, or boundary-context label in a uniform machine-readable form.

### 2. Delivery-versus-control differential expression

`05_DEseq_contrast_delivery_control.R` starts from raw integer gene counts and verified/split design files. It performs phase-aware DESeq2 normalization and delivery-versus-control contrasts, writing normalized count artifacts and DE result tables for anchor and calibration contexts.

This stage supplies differential-expression tables to the locked-anchor coefficient construction steps. It is a count-level analysis step and is not rerun during reviewer-facing manuscript-output reproduction.

### 3. Locked-anchor retained-gene construction

The core retained-gene reconstruction sequence is:

| Script | Function | Principal output |
| --- | --- | --- |
| `06A_core_gene_set.R` | Builds locked-anchor support/core gene sets from Step05 DE outputs | `05_score/anchors/core_gene_set.tsv`, support and contrast-count tables |
| `06B_gene_heterogeneity.R` | Summarizes cross-anchor effect estimates and heterogeneity | `05_score/anchors/gene_heterogeneity.tsv` |
| `06C_Power_analysis.R` | Calculates information/power annotation from anchor estimates | `05_score/anchors/gene_power.tsv` |
| `07_weight_estimation.R` | Applies the existing frozen coefficient estimation rule | `05_score/anchors/gene_weights.tsv` |

The release-facing canonical coefficient table is `data/derived/frozen_gene_weights.tsv`, which corresponds to the frozen `gene_weights.tsv` used for manuscript scoring and interpretation. It is not edited manually in the workflow.

The `frozen_gene_reconstruction/` subfolder retains a focused copy of Steps 06A-07 plus an optional comparison template. `check_frozen_gene_weights_reproducibility_TEMPLATE.R` is intended to compare a future regenerated coefficient table against the released canonical file. It is never called by the default reviewer runner.

### 4. Frozen sample-level scoring

`08_score_samples.R` applies frozen anchor weights to datasets without refitting coefficients. It reads count matrices, scoring designs, and `gene_weights.tsv`; performs dataset-internal normalization and control-referenced standardization; and writes sample-level IMRS score and QC/contributor tables under the original `05_score/transfer/` workflow location.

The relevant framework concept is frozen scoring: validation or transfer-evaluation datasets do not update the retained gene list or frozen coefficients.

### 5. Step09 split-level evaluation

`09_calibration_evaluation.R` combines Step08 sample-level scores with verified split designs and writes:

- `step09_split_eval.tsv`
- `step09_split_summary.tsv`
- `step09_split_sample_level.tsv`

This stage evaluates delivery-minus-control Delta IMRSz, directionality, and secondary AUC by split contrast. It is evaluation only and does not refit IMRS or alter frozen weights.

Released Step09-derived and supporting tables are the bridge from Layer 1 into Layer 2 manuscript-output reproduction.

## Canonical Mode Scope

The portable canonical runner has been validated using the local-example
configuration profile:

```sh
Rscript scripts/portable_full_pipeline/run_count_to_v6_outputs.R \
  --config config/local_d_drive_config_TEMPLATE.yml \
  --mode canonical
```

That configuration selects the five locked anchor datasets:

```text
GSE167521
GSE264344
GSE279372
GSE279744
GSE39129
```

The successful canonical run executed input checks, original metadata/design
and split generation, DESeq2 anchor contrasts, Steps 06A-07 frozen-weight
construction, released-weight comparison, Step08 scoring, Step09 anchor-level
evaluation, and v6 figure-input packaging. It reported:

```text
PASS: regenerated frozen weights match the released canonical table within tolerance.
```

Because canonical mode is scoped to the locked-anchor reconstruction,
non-anchor validation datasets such as `GSE119119`, `GSE139529`, `GSE166655`,
`GSE178313`, `GSE279743`, `GSE314070`, and
`GSE262515_cell_line`/`GSE262515_tissue` do not receive newly generated
Step08 score files in that run and are consequently skipped by Step09. This
is expected behavior, not an error.

Generated products are written below the configured `output_root` (for the
local-example profile, `portable_full_pipeline_v6`). Released `data/derived/`
tables remain read-only and are not overwritten. This reconstructs the
locked-anchor frozen model and packages derived inputs; full manuscript
output reproduction is performed by Layer 2 from released derived
scoring/evaluation tables.

## All-Scored Mode And Step09-To-Layer2 Bridge

All-scored mode first reconstructs weights from exactly the same five locked
anchors, including the released-weight comparison, and then runs the original
Step08 and Step09 application/evaluation logic for every configured dataset
whose prepared counts and designs are available. It does not add validation
datasets to weight estimation. The wrapper stages `GSE262515` as the curated
`GSE262515_cell_line` and `GSE262515_tissue` manuscript arms used by the
original downstream scripts.

```sh
Rscript scripts/portable_full_pipeline/run_count_to_v6_outputs.R \
  --config config/local_d_drive_config_TEMPLATE.yml \
  --mode all_scored

Rscript scripts/portable_full_pipeline/run_step09_to_layer2_inputs.R \
  --config config/local_d_drive_config_TEMPLATE.yml \
  --mode all_scored
```

The second command ports the historically distributed handoff across the
audit and publication-extra scripts, producing manuscript-facing inputs
under:

```text
<output_root>/layer2_generated_inputs/data/derived/
<output_root>/layer2_generated_inputs/data/derived/figure_inputs/
```

A Layer 2 preflight can be pointed to the generated package with:

```sh
Rscript run_all_manuscript_outputs_v6.R \
  --config <output_root>/layer2_generated_inputs/layer2_generated_inputs_config.yml
```

The generated bridge manifest distinguishes regenerated tables from active
Layer 2 inputs that still require released reference tables. Released
`data/derived/` files are never overwritten by either mode.

### 6. Optional diagnostic/context scripts

The following scripts preserve analyses or checks useful for provenance, but they are not required to derive the canonical frozen coefficient table or Step09 tables:

- `Post-09_Check.R`: post-evaluation summaries/pass-check outputs.
- `audit_gse264344.R`: dataset-specific stratified audit.
- `compare_imrs_to_literature.R`: historical literature-context comparison.
- `post_09A_flag_low_imrs_datasets.R`, `post_09B_diagnose_weak_datasets.R`, `post_09C_Noise_Rescore.R`: historical diagnostic/sensitivity materials.

Several of these scripts use historical wording that is not preferred for manuscript-facing text. They are retained as provenance material and should be terminology-reviewed before public interpretation or execution.

## Outputs Handed To Layer 2

By default, Layer 2 consumes released, read-only derived inputs rather than
rerunning Layer 1. After an explicit all-scored and bridge run, it can instead
be configured to read the separate generated input package. Important handoff
files include:

- `data/derived/frozen_gene_weights.tsv`
- `data/derived/gene_power.tsv`
- `data/derived/gene_symbol_mapping.tsv`
- `data/derived/supplement_dataset_split_provenance_v7.tsv`
- `data/derived/weak_dataset_paper_context_audit.tsv`
- `data/derived/label_permutation_null_summary.tsv`
- `data/derived/leave_one_gene_out_summary.tsv`
- `data/derived/gene_dominance_summary.tsv`
- `data/derived/threshold_sensitivity_summary.tsv`
- `data/derived/leave_one_anchor_out_summary.tsv`
- `data/derived/baseline_signature_summary_by_group.tsv`
- `data/derived/coefficient_sensitivity_summary.tsv`
- `data/derived/figure_inputs/`, including Step09 and other figure-source tables

These released tables support manuscript-output reproduction without changing the IMRS score, validation results, robustness results, enrichment statistics, or manuscript claims.

## Execution Status and Portability

The portable Layer 1 canonical count-level path has been executed successfully
for the five locked anchors and validated against the released canonical
frozen-weight table. Historical scripts retained directly in this folder
remain provenance materials; the supported portable entry point is
`scripts/portable_full_pipeline/run_count_to_v6_outputs.R`.

The historical provenance scripts now use a relative/current-directory project-root
default or accept an explicit project-root argument. For an executable,
configuration-driven reconstruction from prepared counts, use
`scripts/portable_full_pipeline/run_count_to_v6_outputs.R`; it executes
path-parameterized copies of the original scripts under
`scripts/portable_full_pipeline/original_ported/`, adds canonical anchor checks
and staged outputs, and treats released derived tables as read-only.

`Post-09_Check.R` remains an optional historical audit script and its dependency
behavior should be reviewed before independent use. `.Rhistory` is local state
and is not pipeline code.

## Safety and Interpretation Scope

Layer 1 reconstructs a frozen, transfer-oriented IMRS scoring framework for acute delivery-associated innate transcriptional responses. It does not establish causal pathways, cell-type sources, clinical reactogenicity prediction, or delivery-platform safety ranking.
