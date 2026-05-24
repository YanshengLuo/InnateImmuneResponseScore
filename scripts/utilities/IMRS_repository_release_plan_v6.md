# IMRS Repository Release Plan v6

## Recommended Repository Structure

```text
IMRS/
  README.md
  LICENSE
  CITATION.cff
  renv.lock or environment.yml
  data/
    external_accessions/
    metadata/
    derived/
  scripts/
    01_metadata/
    02_design/
    03_anchor_model/
    04_scoring/
    05_transfer_evaluation/
    06_robustness/
    07_figures/
    08_supplementary_tables/
  results/
    figures/
    supplementary_figures/
    supplementary_tables/
    logs/
  docs/
    reproducibility/
    data_dictionary/
```

## Public GitHub Contents

Include analysis scripts, metadata-processing scripts, split-design builders, locked-anchor coefficient construction code, frozen scoring code, transfer-evaluation code, robustness scripts, comparator-signature benchmarking, retained-gene enrichment, figure generation, supplementary table packaging, manifests, data dictionaries, and documentation.

## Zenodo or Equivalent Archive

Archive a versioned release matching the submitted manuscript. Include derived tables, frozen gene weights, verified metadata and split-design definitions where redistribution is permitted, figure outputs, supplementary tables, reproducibility manifests, logs, and session information.

## Data Not to Redistribute Directly

Do not redistribute raw public GEO data if the preferred journal/repository practice is to retrieve them from their original public accessions. Instead, provide accession numbers, retrieval scripts, metadata curation, split definitions, and checksums or provenance notes sufficient to rebuild the analysis inputs.

## Review Access

If the repository remains private before acceptance, provide a reviewer-access link or archive. The reviewer package should include the v6 outputs, scripts, manifests, data/code availability statements, and this NAR G&B readiness folder.

## Scope Statement

The release documents the frozen, transfer-oriented IMRS scoring framework and delivery-minus-control ΔIMRSz analyses. It does not frame IMRS as a mechanistic pathway model, clinical reactogenicity predictor, or delivery-platform safety ranking tool.
