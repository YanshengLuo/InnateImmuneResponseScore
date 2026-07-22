# Innate Immune Response Score (IMRS) reproducibility package

This package contains the cleaned count-level inputs, curated metadata, executable R scripts, configuration templates, release documentation, and final publication-facing supplementary TSV files supporting the IMRS manuscript.

## Scope

The executable workflow begins from cleaned integer gene-count matrices and curated sample metadata. Raw sequencing retrieval, alignment, and gene-level quantification are documented for provenance but are not rerun by the default count-to-output route.

The frozen production model uses five acute mouse anchor datasets. Validation and transfer datasets are scored without coefficient refitting, tuning, or recalibration. The publication-facing analysis contains 70 unique matched delivery-versus-control contrasts; two additional rows in Supplementary Table S2 are metadata-only and are not scored.

## Starting points

- `run_all_manuscript_outputs_v6.R` regenerates staged manuscript outputs when the required inputs and R environment are available.
- `run_full_pipeline_from_public_data_TEMPLATE.R` is the portable count-level pipeline template.
- `config/config_template.yml` and `config/full_pipeline_config.yml` provide portable configuration examples.
- `scripts/` contains the implementation and supporting utilities.
- `data/` contains included count-level inputs, curated metadata, and derived source tables.
- `publication_outputs/` contains the synchronized final machine-readable Supplementary Tables S1-S5 and notes.

## Environment

Package versions are recorded in `renv.lock`. Restore dependencies with `renv::restore()` in a suitable R 4.3 environment. The library cache itself is intentionally excluded from this archive.

## Interpretation boundary

IMRS is a frozen bulk RNA-seq framework for measuring an acute delivery-associated innate transcriptional response axis. It is not a mechanistic pathway model, clinical reactogenicity predictor, adverse-event predictor, or universal delivery-platform safety metric.

## Licenses

See `LICENSE` and `DATA_LICENSE.md`. Public source datasets remain subject to their originating repositories and publications.
