# Portable Dataset Reliability Framework with IMRS Plugin

## Purpose
Build a **portable dataset-reliability core** for transcriptomic analyses that does **not** depend on IMRS outputs, while allowing **IMRS-specific evidence** to be added as an optional plugin.

This keeps the method portable across external pipelines and preserves IMRS-specific strengths as extensions rather than requirements.

---

## Design Principles

1. **Portable core first**
   - Required inputs should be available from most RNA-seq workflows after quantification.
   - The core must not require a downstream score.

2. **Plugins are optional**
   - IMRS-specific features such as sample scores, outlier-rescore behavior, and frozen gene weights belong in a plugin.
   - Other analyses can later add their own plugins.

3. **Simple v1 structure**
   - Keep file contracts explicit.
   - Use numbered scripts and numbered input files.
   - Avoid overengineering.

4. **Dataset-level output, multi-level evidence**
   - Final output is one confidence/reliability summary per dataset.
   - Evidence may come from sample-level expression behavior, dataset metadata, and optional analysis-specific outputs.

---

## Top-Level Architecture

### Layer 1: Portable Core
Required or broadly available inputs:
- dataset manifest
- sample metadata
- expression matrix
- optional QC summary
- optional gene effects

The portable core computes dataset reliability from:
- replicate/sample coherence
- variance structure
- outlier burden in expression space
- metadata/domain fit
- optional biological concordance from gene-level effects

### Layer 2: Optional Plugins
Plugins add analysis-specific evidence.

Initial plugin:
- **IMRS plugin**

IMRS plugin can use:
- IMRS sample scores
- IMRS outlier-rescore summary
- IMRS frozen gene weights
- IMRS-specific concordance metrics

---

## Folder Structure

```text
scripts/
├── confidence_core/
│   ├── 01_build_dataset_manifest.R
│   ├── 02_build_sample_metadata.R
│   ├── 03_build_expression_matrix.R
│   ├── 04_build_dataset_qc_summary.R
│   ├── 05_build_gene_effects.R
│   ├── 06_compute_core_confidence.R
│   └── 07_plot_core_confidence.R
│
└── confidence_plugins/
    └── imrs/
        ├── 01_build_imrs_sample_scores.R
        ├── 02_build_imrs_outlier_summary.R
        ├── 03_build_imrs_gene_weights.R
        └── 04_compute_imrs_extension.R
```

```text
05_score/
└── confidence_framework/
    ├── core/
    │   ├── inputs/
    │   ├── intermediate/
    │   ├── outputs/
    │   ├── plots/
    │   └── logs/
    │
    └── plugins/
        └── imrs/
            ├── inputs/
            ├── intermediate/
            ├── outputs/
            ├── plots/
            └── logs/
```

---

## Portable Core Input Contracts

### 01_dataset_manifest.tsv
One row per dataset.

Suggested columns:
- `dataset_id`
- `dataset_class`
- `organism`
- `tissue`
- `assay_type`
- `time_hr`
- `treatment_type`
- `grouping_variable`
- `include_in_analysis`
- `notes`

Purpose:
- define the dataset universe
- carry basic biological and experimental context
- support metadata/domain-fit logic

### 02_sample_metadata.tsv
One row per sample.

Suggested columns:
- `dataset_id`
- `sample_id`
- `group`
- `replicate`
- `batch` (optional)
- `time_hr` (optional if sample-specific)
- `notes` (optional)

Purpose:
- define sample grouping
- support within-group variance and replicate coherence analyses

### 03_expression_matrix.tsv
Gene-by-sample quantified expression matrix.

Format:
- first column: `gene_id`
- remaining columns: sample IDs

Purpose:
- compute sample correlation structure
- compute within-group spread
- detect expression-space outliers
- support PCA/centroid-based coherence metrics

Preferred source:
- post-quantification, cleaned gene-level matrix
- can be normalized expression if consistently defined
- should not be raw FASTQ

### 04_dataset_qc_summary.tsv
Optional but recommended.
One row per sample.

Suggested columns:
- `dataset_id`
- `sample_id`
- `read_depth`
- `mapping_rate`
- `assigned_rate`
- `detected_genes`
- `duplication_rate` (optional)

Purpose:
- support technical QC-aware reliability summaries
- should remain optional because not all external pipelines expose the same QC fields

### 05_gene_effects.tsv
Optional.
One row per gene per dataset.

Suggested columns:
- `dataset_id`
- `gene_id`
- `effect`
- `padj` (optional)
- `stat` (optional)
- `effect_type`

Examples of `effect`:
- log2FC
- standardized contrast effect

Purpose:
- support biological coherence/concordance metrics in the core
- should not be required for the portable core to run in minimal mode

---

## Core Outputs

### dataset_core_confidence.tsv
One row per dataset.

Suggested columns:
- `dataset_id`
- `confidence_core`
- `confidence_metadata`
- `confidence_variance`
- `confidence_outlier`
- `confidence_biological_core`
- `confidence_class`
- `primary_failure_mode`

### dataset_core_flags.tsv
Suggested columns:
- `dataset_id`
- `flag_type`
- `flag_value`
- `flag_message`

### dataset_core_for_plotting.tsv
Merged plotting-ready summary.

---

## IMRS Plugin Input Contracts

### 01_imrs_sample_scores.tsv
One row per sample.

Suggested columns:
- `dataset_id`
- `sample_id`
- `condition`
- `imrs_raw`
- `imrs_z`

Purpose:
- evaluate score-specific stability and separation

### 02_imrs_outlier_summary.tsv
One row per dataset.

Suggested columns:
- `dataset_id`
- `before_delta_imrs`
- `after_delta_imrs`
- `delta_imrs_change_abs`
- `remove_n`
- `removed_samples`
- `rescore_status`

Purpose:
- evaluate IMRS-specific outlier sensitivity

### 03_imrs_gene_weights.tsv
One row per gene.

Suggested columns:
- `gene_id`
- `weight`
- `variance` (optional)
- `source` (optional)

Purpose:
- support weight-aware IMRS diagnostics such as dominance and concordance

---

## IMRS Plugin Outputs

### dataset_imrs_extension.tsv
One row per dataset.

Suggested columns:
- `dataset_id`
- `confidence_imrs_extension`
- `imrs_score_stability`
- `imrs_outlier_sensitivity`
- `imrs_weight_concordance`
- `imrs_dominance_penalty`

### dataset_imrs_flags.tsv
Flag-style plugin diagnostics.

---

## Combined Reporting Strategy

The framework may report:
- `confidence_core`
- `confidence_imrs_extension` (if plugin available)
- `confidence_combined` (optional later)

For v1, the safest reporting is:
- always output the core
- output plugin extension separately when available
- do not force a combined score unless the behavior is clearly justified

---

## Current Project Decision

For this project:
- the **portable core** is the primary architecture
- anything using IMRS score outputs belongs in the **IMRS plugin**
- we will now begin by constructing the **portable core inputs**

Construction order:
1. `01_dataset_manifest.tsv`
2. `02_sample_metadata.tsv`
3. `03_expression_matrix.tsv`
4. `04_dataset_qc_summary.tsv` (if available)
5. `05_gene_effects.tsv` (optional)

After that, we can define the IMRS plugin input builders.

---

## Key Constraint to Preserve

The framework must remain portable even if an external analysis has:
- no downstream score
- no gene weights
- no rescore workflow

Therefore:
- scores and weights must never be required by the portable core
- they are plugin-only inputs

