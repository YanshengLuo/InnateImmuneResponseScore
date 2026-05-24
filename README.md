# IMRS Repository Release Candidate

IMRS is a frozen, transfer-oriented transcriptomic scoring framework for acute delivery-associated innate transcriptional responses.

IMRS is not a mechanistic pathway model, clinical reactogenicity predictor, or delivery-platform safety ranking tool.

## Reproducibility Levels

This release candidate separates reproducibility into three levels.

Level 1 is the default reviewer-facing manuscript-output reproduction path. It starts from staged derived inputs in `data/derived`, not from raw metadata, FASTQ, BAM, or full count reconstruction. The Level 1 runner regenerates Supplementary Figure S2, Supplementary Tables S1-S5, and NAR G&B readiness materials.

Level 2 is analytic reconstruction from gene count matrices, verified metadata, split definitions, and locked-anchor inputs. It documents how count-level inputs feed into locked-anchor differential expression, retained gene selection, frozen IMRS gene coefficients, sample-level scoring, delivery-minus-control &Delta;IMRSz, validation, and robustness outputs.

Level 3 is raw-data/HPC reconstruction from public GEO/SRA accessions or raw public data through FASTQ/BAM retrieval, alignment/counting, and gene count matrix creation. This path may require HiPerGator/HPC and is documented separately from the default reviewer runner.

## Reviewer-Facing Runner

Use `run_all_manuscript_outputs_v6.R` for Level 1 reviewer-facing reproduction. The active runner regenerates:

- Supplementary Figure S2 retained-gene enrichment outputs
- Supplementary Tables S1-S5
- NAR G&B readiness/reproducibility materials

Generated outputs are written to `results_release_templates/` by default and are ignored by Git.

## Configuration

Before running, copy:

```sh
cp config/config_template.yml config/config.yml
```

Edit `config/config.yml` if needed. Keep `execute_active_scripts: false` for checklist/dry-run mode. Set `execute_active_scripts: true` only for controlled regeneration.

The clean-release layout expects staged Level 1 derived inputs under `data/derived`.

## Recommended Reviewer Command

```sh
Rscript run_all_manuscript_outputs_v6.R
```

## Curated Metadata and Split Definitions

Curated metadata and split-design files are included for transparency under `data/curated_metadata`, `data/split_designs`, and `docs/manual_curation`. They document treatment/control group definitions, tissue/timepoint labels, manuscript role assignments, boundary categories, and publication-context mapping. The Level 1 runner does not directly rebuild all metadata/split definitions from raw public metadata; it consumes staged derived manuscript inputs that already encode these decisions. Level 2 reconstruction documents how metadata and raw/public data feed into scoring and validation.

Manual curation does not manually alter delivery-minus-control &Delta;IMRSz values.

## Full Framework and Raw-Data Reconstruction

Full count-to-score reconstruction and raw-data/HPC reconstruction are documented separately. The reviewer-facing runner does not run raw-data retrieval, alignment/counting, frozen-gene reconstruction, HiPerGator/SLURM scripts, archive/legacy scripts, or executable `full_pipeline` scripts.

## Main Figure Outputs

Main figure regeneration is optional and legacy-dependent in this release. Submitted main figure outputs should be provided as release artifacts. The active runner currently focuses on regenerating Supplementary Figure S2, Supplementary Tables S1-S5, and readiness materials.

## Release Notes

Before public release, replace placeholder repository/Zenodo links, choose a real license, and add `renv.lock` or `environment.yml` for software environment capture. Raw public sequencing data should generally be referenced by accession rather than redistributed unless redistribution is allowed and practical.
