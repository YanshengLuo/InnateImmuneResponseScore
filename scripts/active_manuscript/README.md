# Active Manuscript Layer

This directory implements Layer 2 reviewer-facing manuscript-output reproduction from released derived inputs.

## Entry Points

All active entry scripts now locate the repository from their own file path (including RStudio Run/Source) and do not require the working directory to be changed manually. For figure-only regeneration, `RUN_FIGURES_v6.R` at the repository root is the simplest entry point.

- `00_generate_manuscript_figures_v6.R`: regenerates manuscript figure panels and assembled figures from `data/derived/figure_inputs/`.
- `02_run_priority3_gene_program_enrichment_v6.R`: regenerates Supplementary Figure S2 and retained frozen-gene enrichment tables from released frozen inputs.
- `01_build_supplementary_tables_v6.R`: packages Supplementary Tables S1-S5 after Priority3 enrichment outputs are present.

## Local Implementation Library

- `lib/active_config.R`: resolves repository-relative configuration.
- `lib/figure_helpers_v6.R`: the formerly external `v5_helpers.R` implementation, now included in the active release layer.
- `lib/panel_builders_v6.R`: the panel-construction functions loaded by the helper, refactored to read released figure inputs only.
- `lib/merged_workflow_v6.R`: Figure 1A workflow renderer.

The default runner at repository root calls only the three active entry points above. Optional internal readiness packaging and prior launchers are outside this directory and are not required for reviewer execution.
