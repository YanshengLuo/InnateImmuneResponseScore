# Reviewer Reproduction Status

## Current Status

The reviewer-facing reproduction layer is functional for the staged manuscript-output package. The public release candidate uses a two-layer structure:

- Layer 1 documents full computational provenance from public gene-count and verified-metadata inputs to frozen IMRS scoring and Step09-derived result tables.
- Layer 2 is the default reviewer-facing quick run from released derived inputs to manuscript figures, Supplementary Figure S2 Priority3 gene-program enrichment outputs, and Supplementary Tables S1-S5.

Layer 2 contains its active implementation within this repository, including the figure helper code formerly sourced from an external `v5_helpers.R` file.

## Default Reviewer Run

Run the repository-root script:

```sh
Rscript run_all_manuscript_outputs_v6.R
```

The default configuration keeps `execute_active_scripts: false`, which writes a preflight checklist without generating scientific outputs. After reviewing that checklist, setting `execute_active_scripts: true` runs the required Layer 2 steps:

- `manuscript_figures`
- `priority3_gene_program_enrichment`
- `supplementary_tables_s1_s5`

Generated outputs are written below `results_release_templates/`.

## Optional/Internal Material

`internal_readiness` is optional and disabled by default through `run_internal_readiness: false`. NAR G&B readiness packaging is retained only as optional internal documentation and is not required for reviewer reproduction.

The default Layer 2 runner does not execute raw-data retrieval, alignment/counting, HiPerGator/SLURM work, Layer 1 full-pipeline reconstruction, archive/legacy scripts, or frozen-gene reconstruction.

## Validation

Layer 2 validation recorded in `docs/reproducibility/layer2_runner_validation_summary.txt` confirms:

- checklist mode passed with all required Layer 2 steps ready and internal readiness disabled;
- controlled Layer 2 execution regenerated manuscript figures, Priority3 gene-program enrichment outputs, and Supplementary Tables S1-S5 under repository-local output folders;
- no full-pipeline, HPC, archive/legacy, raw-data, or frozen-gene reconstruction steps were executed.

A final pre-push checklist run on the Git worktree also passed in non-executing mode with the three required Layer 2 steps ready and internal readiness disabled.

## Scope

The repository provides a reproducible manuscript-output layer from released derived inputs while preserving documentation for upstream reconstruction. Manual curation supplies metadata, split, role, and context definitions; it does not manually alter delivery-minus-control Delta IMRSz values.

No IMRS scoring, validation, robustness, enrichment statistics, figure-data values, or supplementary-table source values were intentionally changed as part of release packaging.
