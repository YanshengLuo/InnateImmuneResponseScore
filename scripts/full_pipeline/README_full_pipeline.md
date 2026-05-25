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

### 6. Optional diagnostic/context scripts

The following scripts preserve analyses or checks useful for provenance, but they are not required to derive the canonical frozen coefficient table or Step09 tables:

- `Post-09_Check.R`: post-evaluation summaries/pass-check outputs.
- `audit_gse264344.R`: dataset-specific stratified audit.
- `compare_imrs_to_literature.R`: historical literature-context comparison.
- `post_09A_flag_low_imrs_datasets.R`, `post_09B_diagnose_weak_datasets.R`, `post_09C_Noise_Rescore.R`: historical diagnostic/sensitivity materials.

Several of these scripts use historical wording that is not preferred for manuscript-facing text. They are retained as provenance material and should be terminology-reviewed before public interpretation or execution.

## Outputs Handed To Layer 2

Layer 2 consumes released, read-only derived inputs rather than rerunning Layer 1. Important handoff files include:

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

This folder is documentation-first. The scripts are included to preserve the computational path, but full execution has not been performed as part of release packaging.

Known portability issues include:

- Several scripts retain the historical default `D:/IMRS_Project` path or other original workflow directory assumptions.
- `03_Final_design_file_autogroup.R` has a directly assigned absolute metadata input path.
- `01_run_all_metadata_build.R` sources metadata code from the historical script-tree layout.
- `Post-09_Check.R` attempts package installation during execution.
- The frozen-weight comparison template still contains historical staged-path defaults that should be aligned with `data/derived/` before future use.
- `.Rhistory` is a local-state artifact and should not be committed as pipeline code.

Before a full Layer 1 rerun, these scripts should be updated to use a configuration file or command-line project root, tested against released metadata/split definitions and count sources, and provided with an environment lockfile.

## Safety and Interpretation Scope

Layer 1 reconstructs a frozen, transfer-oriented IMRS scoring framework for acute delivery-associated innate transcriptional responses. It does not establish causal pathways, cell-type sources, clinical reactogenicity prediction, or delivery-platform safety ranking.
