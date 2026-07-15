# Released Figure Inputs

This folder contains the released derived score, support, robustness, comparator, and audit tables read by `scripts/active_manuscript/00_generate_manuscript_figures_v6.R`.

These tables are inputs for Layer 2 manuscript-figure regeneration. They are not raw RNA-seq counts and are not used to refit or redefine the frozen IMRS score. Layer 1/full-pipeline scripts document upstream reconstruction from gene counts and verified metadata.

For public release hygiene, historical local source-path prefixes in provenance fields were replaced with `<upstream_project_root>`. This replacement affects provenance text only; no numeric/statistical values, term identifiers, gene lists, row counts, or column counts were changed.
