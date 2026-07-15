# IMRS Reproducibility Release Candidate

IMRS is a frozen, transfer-oriented transcriptomic scoring framework for acute delivery-associated innate transcriptional responses.

IMRS is not a mechanistic pathway model, clinical reactogenicity predictor, or delivery-platform safety ranking tool.

## Two-layer reproducibility structure

### Layer 1: Full computational provenance

This layer documents and scripts the path from public RNA-seq count matrices and verified metadata to frozen gene weights, sample-level IMRS scores, and Step09 split-level evaluation tables. It is not the default reviewer run because it may require large public data files, DESeq2 processing, and longer-running reconstruction steps.

Layer 1 scripts are organized under `scripts/full_pipeline/`. Raw-data retrieval and alignment/counting on HiPerGator/SLURM are upstream provenance operations within this full reconstruction context and are not invoked by the default runner. Account-specific HPC submission scripts are intentionally excluded from this public release candidate and can be distributed later only as sanitized templates.

## Three Practical Run Routes

### A. Reviewer quick run from released derived inputs

```sh
Rscript run_all_manuscript_outputs_v6.R
```

This default route checks or regenerates manuscript outputs from the released,
read-only `data/derived/` bundle according to the Layer 2 config.

### B. Canonical count-level reconstruction

```sh
Rscript scripts/portable_full_pipeline/run_count_to_v6_outputs.R --config config/full_pipeline_config.yml --mode canonical
```

This route reconstructs the frozen IMRS model from the five locked anchors
and validates the released canonical coefficient table. Its Step08/Step09
scope is the locked-anchor set.

### C. Full generated-input route

```sh
Rscript scripts/portable_full_pipeline/run_count_to_v6_outputs.R --config config/full_pipeline_config.yml --mode all_scored
Rscript scripts/portable_full_pipeline/run_step09_to_layer2_inputs.R --config config/full_pipeline_config.yml --mode all_scored
Rscript run_all_manuscript_outputs_v6.R --config <output_root>/layer2_generated_inputs/layer2_generated_inputs_config.yml
```

This route reconstructs weights from the locked anchors only, scores
available manuscript validation datasets with those frozen coefficients, and
regenerates available manuscript-facing Layer 2 inputs through the ported
distributed bridge. Validation datasets never enter Steps 06A-07 and never
refit the model.

### Reviewer count-level validation from cleaned matrices

The reviewer package includes cleaned raw integer count matrices, not
normalized expression values. These are featureCounts-derived count matrices
cleaned for a portable rerun under
`data/counts/{dataset}/featurecounts/validation/gene_counts_clean.tsv`;
reviewers do not need to rerun FASTQ download, alignment, or featureCounts.

Use `--mode all_scored` for reviewer validation: frozen weights are
reconstructed only from the five locked anchors, and all configured scored
validation datasets are then scored with those frozen weights without
refitting coefficients.

First check the count-level run plan:

```bat
Rscript scripts\portable_full_pipeline\run_count_to_v6_outputs.R --config config\full_pipeline_config.yml --mode all_scored --dry-run
```

If the dry-run passes, execute the reconstruction:

```bat
Rscript scripts\portable_full_pipeline\run_count_to_v6_outputs.R --config config\full_pipeline_config.yml --mode all_scored --force
```

This count-level route is separate from the Layer 2 reviewer quick run:

```bat
Rscript run_all_manuscript_outputs_v6.R --config config\config.yml
```

## Layer 1 canonical count-level reconstruction

The default reviewer runner regenerates manuscript figures and tables from
released derived inputs. The count-level runner at
`scripts/portable_full_pipeline/run_count_to_v6_outputs.R` starts from
prepared gene-count matrices and verified metadata. The portable Layer 1
canonical mode reconstructs the locked-anchor frozen IMRS model from
configured public count matrices and verified metadata. It ports and connects
the original design generation, DESeq2 anchor contrasts, Step06/Step07
frozen-weight reconstruction, Step08 scoring, and Step09 anchor-level split
evaluation rather than reimplementing those calculations.

A canonical Layer 1 run for the five locked anchors completed successfully and
reported: "PASS: regenerated frozen weights match the released canonical table
within tolerance." Canonical mode is intended to validate the frozen-weight
construction and regenerate anchor-level scoring/evaluation outputs and
derived figure-input packages. It does not recompute all primary, extended,
or secondary manuscript validation score tables from counts. Full manuscript
figure/table reproduction remains handled by Layer 2 from released derived
scoring/evaluation tables.

This workflow does not retrieve FASTQ/BAM files, align reads, or perform
feature counting. It never overwrites canonical released tables under
`data/derived/`; regenerated products are written beneath
the configured `output_root` (by default `results/full_pipeline_v6/`).
Canonical reconstruction requires the five locked anchor datasets. A subset
run under `--mode test` is labelled non-canonical/test-only and does not
reconstruct frozen weights unless the complete original anchor requirement is
met; single-anchor modeling stops before Step06A. See
`docs/count_level_full_pipeline.md` and
`docs/layer1_canonical_validation_summary.md` for the input contract,
validated scope, configuration, and exact run commands.

### Layer 1 all-scored application and generated Layer 2 inputs

`--mode all_scored` reconstructs the same frozen model from the five locked
anchors and then applies those frozen coefficients to every configured
manuscript dataset with available prepared counts and scoring designs.
Validation datasets do not participate in Steps 06A-07 and do not refit or
alter the frozen model. For `GSE262515`, the portable wrapper stages its
curated cell-line and tissue manuscript arms for the unmodified Step08/Step09
scripts.

### Production, strict-3 sensitivity, and scoring coefficient scope

The production frozen model uses exactly five locked acute mouse anchors:
`GSE39129`, `GSE167521`, `GSE264344`, `GSE279372`, and `GSE279744`.
Production core-gene selection follows this five-anchor workflow;
`GSE262515` is validation/secondary support only and does not contribute to
Steps 06A-07.

The strict-3 set is `GSE39129`, `GSE167521`, and `GSE264344`. Strict-3 and
strict-3 threshold-sensitivity analyses are supporting sensitivity/ablation
checks only; they do not define, replace, or refit the production five-anchor
frozen model.

Step08 applies `beta_meta` whenever that column is present and copies the
penalized `weight` field into `beta_meta` only when `beta_meta` is absent. The
released canonical table contains `beta_meta`, so released scores use
`beta_meta`; `weight` remains an audit/fallback field. Validation datasets do
not refit or update these coefficients.

After an `all_scored` run, the distributed historical Step09-to-Layer2
handoff can be regenerated with the portable bridge runner:

```sh
Rscript scripts/portable_full_pipeline/run_count_to_v6_outputs.R --config config/full_pipeline_config.yml --mode all_scored
Rscript scripts/portable_full_pipeline/run_step09_to_layer2_inputs.R --config config/full_pipeline_config.yml --mode all_scored
Rscript run_all_manuscript_outputs_v6.R --config <output_root>/layer2_generated_inputs/layer2_generated_inputs_config.yml
```

The bridge invokes path-parameterized copies of the original audit and
publication-extra analysis scripts and writes a separate generated package
below `<output_root>/layer2_generated_inputs/`. Inputs required by the active
Layer 2 figures that do not have a directly ported bridge generator are
identified as released references in its manifest and input contract. The
default Layer 2 reviewer run continues to use released `data/derived/`
tables unless it is explicitly configured to use this generated package.

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
