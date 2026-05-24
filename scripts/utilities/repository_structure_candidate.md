# Candidate Public Repository Structure

The following layout is recommended for a public IMRS repository release. It separates executable code, curated/manual inputs, derived results, logs, and reproducibility documentation.

```text
IMRS/
  README.md
  LICENSE
  CITATION.cff
  renv.lock or environment.yml
  config/
    config_template.yml
  data/
    accessions/
    curated_metadata/
    split_designs/
    derived/
  scripts/
    00_setup_environment/
    01_metadata_processing/
    02_design_and_split_generation/
    03_anchor_model/
    04_frozen_scoring/
    05_transfer_evaluation/
    06_robustness/
    07_comparator_signatures/
    08_gene_program_enrichment/
    09_figures/
    10_supplementary_tables/
  results/
    figures/
    supplementary_figures/
    supplementary_tables/
    manifests/
    logs/
  docs/
    reproducibility/
    manual_curation/
    data_dictionary/
```

## Notes

- Keep raw public-data retrieval instructions and accession tables under `data/accessions/`.
- Keep curated delivery/control labels, tissue/timepoint annotations, manuscript roles, and split definitions under `data/curated_metadata/` and `data/split_designs/`.
- Keep frozen scoring outputs and derived evaluation tables under `data/derived/` or `results/manifests/`, depending on release preference.
- Move final plot outputs currently located inside script-oriented folders into `results/figures/` or `results/supplementary_figures/`.
- Use `config/config_template.yml` to replace local absolute paths and versioned manuscript output roots.
