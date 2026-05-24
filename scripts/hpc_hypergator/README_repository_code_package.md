# IMRS Repository Code Reproducibility Package

Generated: 2026-05-23T17:22:52.128Z

## Purpose

This package inventories code used or available to reproduce the IMRS workflow and manuscript outputs. It documents the frozen, transfer-oriented IMRS scoring framework from public-data processing through locked-anchor coefficient construction, frozen scoring, transfer evaluation, robustness checks, gene-program enrichment, figure generation, supplementary table packaging, and reviewer-facing reproducibility documentation. This is an inventory and documentation package only; it does not rerun scoring, validation, robustness, enrichment, figure, or table-generation analyses.

## How to Set project_root

Before public repository release, replace local absolute paths such as `<local_project_root>` and HiPerGator paths such as `<HiPerGator_project_path> with a portable project root. Recommended patterns are a `config/config_template.yml`, a command-line `--project_root` argument, or in R scripts `here::here()` anchored at the repository root. Versioned output roots such as `revised_plots_v6` should be command-line or config values.

## What Scripts Reproduce

The inventory spans scripts from the main HiPerGator source tree and v6 manuscript-output folders. The detected workflow covers metadata processing, split design generation, anchor differential expression, retained gene selection, frozen coefficient estimation, frozen IMRS scoring, split-level transfer evaluation, manuscript role cleanup, boundary-context audits for late or context-shifted settings, robustness analyses supporting score stability and non-degeneracy, comparator signatures, retained-gene enrichment, figure generation, supplementary table packaging, and NAR Genomics and Bioinformatics readiness documentation.

## Suggested Run Order

Use `run_order_candidate.tsv` as a static run-order draft. It orders scripts by workflow stage but does not certify that every script can be executed end to end from a clean checkout. The recommended public run path should be refined to a smaller tested set of active scripts, with legacy/prototype/diagnostic scripts moved to an archive folder.

## Public Raw Data

Early workflow stages require retrieval of public raw sequencing data and associated metadata, primarily from public accession records. Raw FASTQ/count retrieval and alignment/counting scripts are represented in the inventory, but public release should provide accession tables and retrieval instructions rather than redistributing public raw data unless permissions allow it.

## Curated and Manual Inputs

Several stages require curated metadata and split-design inputs: delivery/control labels, tissue/timepoint annotations, manuscript analysis groups, and boundary-context annotations. These inputs should be released as versioned tables under `data/curated_metadata/` and `data/split_designs/`, with provenance and curation rationale documented separately.

## Regenerating Manuscript Figures and Tables

The v6 figure-generation scripts can regenerate manuscript-facing figures from derived analysis inputs. The Priority3 enrichment script produces Supplementary Figure S2 and Supplementary Table S5 enrichment outputs. The supplementary-table builder packages manuscript-ready Supplementary Tables S1-S5 from existing derived outputs. Captions assembled outside R should be documented and updated in manuscript or slide/word-processing sources, not forced into scoring scripts.

## HiPerGator Requirements

Scripts with `.slurm` extensions, `#SBATCH` directives, module-loading commands, or `<HiPerGator_project_path> paths are marked as requiring HiPerGator or an equivalent SLURM/HPC environment. Public repository documentation should either retain these as HPC-specific examples or provide local/containerized alternatives for downstream derived-input workflows.

## Known Limitations

- Runnable status is inferred by static text inspection and has not been execution-tested in this packaging task.
- Some scripts contain local absolute paths or version-specific output roots that should be moved to configuration before public release.
- Prototype, diagnostic, GUI, and older support scripts are inventoried but should not be presented as the primary manuscript run path without review.
- Some non-script plot/document outputs were found within script source trees and should be placed under `results/figures` or `results/supplementary_figures` in a clean public repository.
- The package documents code provenance and reproducibility; it does not modify analysis values or manuscript claims.

## Inventory Summary

- Scripts and operational notes inventoried: 73
- Hardcoded path issues: 208
- Scripts requiring HiPerGator/SLURM: 18
- Plot/document outputs found under script source trees: 4

## Stage Counts

```text
00_setup_environment: 16
01_metadata_processing: 2
02_design_and_split_generation: 2
03_anchor_differential_expression: 2
04_locked_anchor_gene_selection: 4
05_weight_estimation: 1
06_frozen_IMRS_scoring: 2
07_split_contrast_transfer_evaluation: 4
08_dataset_role_cleanup: 0
09_boundary_context_audit: 3
10_robustness_permutation: 0
11_robustness_leave_one_gene_out: 0
12_gene_dominance: 1
13_threshold_sensitivity: 0
14_leave_one_anchor_out: 0
15_comparator_signatures: 0
16_gene_program_enrichment: 1
17_figure_generation: 7
18_supplementary_table_generation: 1
19_NAR_GB_readiness: 3
20_manual_or_unknown: 24
```

