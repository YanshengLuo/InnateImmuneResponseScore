# Innate Immune Response Score (IMRS) reproducibility package

This package contains the cleaned count-level inputs, curated metadata, executable R scripts, configuration templates, release documentation, and final publication-facing supplementary TSV files supporting the IMRS manuscript.

## Scope

The executable workflow begins from cleaned integer gene-count matrices and curated sample metadata. Raw sequencing retrieval, alignment, and gene-level quantification are documented for provenance but are not rerun by the default count-to-output route.

The frozen production model uses five acute mouse anchor datasets. Validation and transfer datasets are scored without coefficient refitting, tuning, or recalibration. The final manuscript figures use 68 scored contrasts (40 locked-anchor, 14 primary acute-validation, 11 extended-validation, and 3 secondary-support contrasts). Released upstream tables retain historical rows; some supplementary summaries still describe the earlier 70-contrast state. See [release notes](docs/reproducibility/REPRODUCE_PACKAGE_RELEASE.md) for that known consistency limitation.

## Starting points

- `RUN_FIGURES_v6.R` is the simplest one-click entry point for manuscript figure regeneration. It discovers the repository location automatically in RStudio or Rscript; no manual `setwd()` is required.
- `run_all_manuscript_outputs_v6.R` regenerates staged manuscript outputs when the required inputs and R environment are available.
- `run_full_pipeline_from_public_data_TEMPLATE.R` is the portable count-level pipeline template.
- `config/config_template.yml` and `config/full_pipeline_config.yml` provide portable configuration examples.
- `scripts/` contains the implementation and supporting utilities.
- `data/` contains included count-level inputs, curated metadata, and derived source tables.
- `publication_outputs/` contains supplied Supplementary Tables S1-S5, notes, and the final R-assembled PowerPoint.
- `publication_outputs/figures/` contains Figures 1-6 and Supplementary Figures S1-S2 in SVG, PNG, and PDF.
- `data/derived/gene_program_enrichment/` contains the retained enrichment results and mapping/background audits.

## Reproduce-package branch

This branch preserves the structure and history of `v6-reproducibility-release` and incorporates the latest available local `MANUSCRIPT_FINAL` scripts, corrected permutation results, and final figure snapshots. Temporary authoring files, duplicate ZIP packages, and draft manuscripts are excluded.

```sh
git clone --branch reproduce-package --single-branch https://github.com/YanshengLuo/InnateImmuneResponseScore.git
cd InnateImmuneResponseScore
Rscript run_all_manuscript_outputs_v6.R
```

The command above performs the default preflight. For figure regeneration after restoring dependencies, run `Rscript RUN_FIGURES_v6.R`. It regenerates Figures 1-6 and S1 and restores the released S2 snapshot if no generated S2 is present. To recompute enrichment and supplementary tables, enable `execute_active_scripts: true` in a local `config/config.yml` copied from the template, then run the full manuscript-output runner. Generated files are written under the ignored `results_release_templates/`; committed snapshots remain under `publication_outputs/`.

See [the release inventory](docs/reproducibility/reproduce_package_inventory.tsv) for file sizes and SHA-256 checksums of included scripts, inputs, and outputs.

## Environment

Package versions are recorded in `renv.lock`. Restore dependencies with `renv::restore()` in a suitable R 4.3 environment. The library cache itself is intentionally excluded from this archive.

## Interpretation boundary

IMRS is a frozen bulk RNA-seq framework for measuring an acute delivery-associated innate transcriptional response axis. It is not a mechanistic pathway model, clinical reactogenicity predictor, adverse-event predictor, or universal delivery-platform safety metric.

## Licenses

See `LICENSE` and `DATA_LICENSE.md`. Public source datasets remain subject to their originating repositories and publications.
