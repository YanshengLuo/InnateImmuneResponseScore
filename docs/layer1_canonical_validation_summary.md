# Layer 1 Canonical Validation Summary

## Validation Date

May 25, 2026

## Command

The validation used the local-example configuration profile, which points to
locally available prepared count matrices and verified metadata:

```sh
Rscript scripts/portable_full_pipeline/run_count_to_v6_outputs.R \
  --config config/local_d_drive_config_TEMPLATE.yml \
  --mode canonical
```

The default public configuration remains repository-relative; the local
template is not the reviewer-facing default.

## Scope

Canonical mode configured the five locked anchor datasets:

```text
GSE167521
GSE264344
GSE279372
GSE279744
GSE39129
```

The completed run executed:

- Preflight input checks.
- Original metadata/design preparation and split-design generation.
- DESeq2 contrasts for the five locked anchors.
- Step06A core gene set reconstruction.
- Step06B gene heterogeneity.
- Step06C power analysis.
- Step07 frozen weight estimation.
- Comparison with the released canonical frozen-weight table.
- Step08 frozen sample scoring for the configured locked anchors.
- Step09 split evaluation for those anchor datasets.
- Packaging of regenerated figure inputs for v6 figures.

## Validation Result

```text
PASS: regenerated frozen weights match the released canonical table within tolerance.
```

The run validates reconstruction of the canonical locked-anchor frozen IMRS
model. Generated products are written under the configured output root and do
not overwrite released `data/derived/` inputs.

## Expected Limitation

Canonical mode is not a count-level regeneration of every manuscript
validation score table. Because only the five locked anchors were configured
for scoring, Step09 expectedly skipped non-anchor validation datasets lacking
newly generated Step08 scores, including `GSE119119`, `GSE139529`,
`GSE166655`, `GSE178313`, `GSE279743`, `GSE314070`, and the
`GSE262515` cell-line/tissue validation splits.

Full manuscript figure and table reproduction remains the responsibility of
Layer 2, which consumes the released derived scoring/evaluation tables.
