# IMRS Reproducibility Release Candidate

IMRS is a frozen, transfer-oriented transcriptomic scoring framework for acute delivery-associated innate transcriptional responses.

IMRS is not a mechanistic pathway model, clinical reactogenicity predictor, or delivery-platform safety ranking tool.

## Two-layer reproducibility structure

### Layer 1: Full computational provenance

This layer documents and scripts the path from public RNA-seq count matrices and verified metadata to frozen gene weights, sample-level IMRS scores, and Step09 split-level evaluation tables. It is not the default reviewer run because it may require large public data files, DESeq2 processing, and longer-running reconstruction steps.

Layer 1 scripts are organized under `scripts/full_pipeline/`. Raw-data retrieval and alignment/counting on HiPerGator/SLURM are upstream provenance operations within this full reconstruction context and are not invoked by the default runner. Account-specific HPC submission scripts are intentionally excluded from this public release candidate and can be distributed later only as sanitized templates.

### Layer 2: Reviewer-facing manuscript-output reproduction

This is the default quick run. It regenerates manuscript figures, supplementary tables, and Priority3 gene-program enrichment outputs from released frozen derived inputs and scoring/evaluation tables. It does not rerun raw-data retrieval, HiPerGator/HPC jobs, DESeq2 anchor reconstruction, or full frozen-gene reconstruction.

The current active manuscript figure runner was refactored so that the formerly external `v5_helpers.R` implementation is now included inside the clean release repository under `scripts/active_manuscript/lib/figure_helpers_v6.R`. Its repo-contained panel-builder source and workflow builder are in the same library folder. The data read by this implementation are released under `data/derived/figure_inputs/`.

## Layer 2 Runner

Use `run_all_manuscript_outputs_v6.R` for the reviewer-facing quick run. With active execution enabled it runs:

- `scripts/active_manuscript/00_generate_manuscript_figures_v6.R`
- `scripts/active_manuscript/02_run_priority3_gene_program_enrichment_v6.R`
- `scripts/active_manuscript/01_build_supplementary_tables_v6.R`

Generated outputs are written beneath `results_release_templates/` by default:

- `results_release_templates/figures/`: manuscript figure assemblies and figure manifests
- `results_release_templates/priority3_gene_program_enrichment/`: Supplementary Figure S2 and enrichment tables
- `results_release_templates/supplementary_tables/`: Supplementary Tables S1-S5
- `results_release_templates/logs/` and `results_release_templates/manifests/`: run records and preflight checklist

NAR G&B readiness packaging is retained as optional internal documentation under `scripts/optional_internal/`. It is not required by the default reviewer-facing run and runs only when `run_internal_readiness: true` is set explicitly in a local config.

## Configuration

Before running, copy:

```sh
cp config/config_template.yml config/config.yml
```

Edit `config/config.yml` only if repository-relative defaults need adjustment. Keep `execute_active_scripts: false` for preflight/checklist mode:

```sh
Rscript run_all_manuscript_outputs_v6.R
```

After the checklist confirms all required released inputs are present, set `execute_active_scripts: true` for controlled Layer 2 regeneration. Keep `run_internal_readiness: false` for the default reviewer path.

## Released Inputs

Frozen scoring and enrichment inputs are provided under `data/derived/`, including `frozen_gene_weights.tsv`, `gene_power.tsv`, and the existing robustness/provenance source tables used to package supplementary results.

Figure-generation tables are provided under `data/derived/figure_inputs/`. This explicit bundle includes the released anchor-support, Step09 scoring/evaluation, robustness, comparator, and role/context tables that the repo-contained figure implementation reads. It avoids dependencies on local `revised_plots_v2`, `revised_plots_v3`, or `revised_plots_v4` folders.

Curated metadata and split-design records are included for transparency under `data/curated_metadata/`, `data/split_designs/`, and `docs/manual_curation/`. They define treatment/control group context, tissue/timepoint labels, manuscript roles, and boundary-setting annotations; they do not manually alter delivery-minus-control &Delta;IMRSz values.

## Priority3 Gene-Program Enrichment

Priority3 gene-program enrichment is a required Layer 2 step. It uses the retained frozen IMRS gene set and released background inputs from `data/derived/` and writes Supplementary Figure S2 plus enrichment tables to the configured repository output folder. If an input or R package is unavailable, the step stops or records an explicit diagnostic; it is not silently skipped.

## Scope and Release Notes

The full framework reconstruction path is separate from the reviewer-facing Layer 2 run and may require public-data retrieval, additional storage, and HiPerGator/HPC configuration. Legacy and internal diagnostic material is not part of the active execution path.

Repository code is released under the MIT License; see `LICENSE` and `DATA_LICENSE.md`. Citation metadata identify the public repository without asserting an archival DOI that has not yet been issued. The included `renv.lock` captures the R environment for reproduction. Raw public sequencing data should generally be referenced by accession rather than redistributed unless redistribution is allowed and practical.
