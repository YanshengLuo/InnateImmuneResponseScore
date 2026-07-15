# Layer 2 Figure Implementation Audit

## Resolved Implementation Source

The prior active manuscript figure script was a launcher: it sourced `v5_helpers.R` and delegated generation, assembly, manifest construction, wording audit, validation, and logging. The current v6 source copy of `v5_helpers.R` and the previously staged `scripts/utilities/v5_helpers.R` copy had identical SHA-256 hash:

`332FF38E18BCFD6FB756F2CE618539B73457CE90A1E94D6B1B17E3668D8AE0E6`

The implementation was therefore retained as the current source logic and moved into the active release layer as `scripts/active_manuscript/lib/figure_helpers_v6.R`.

## Direct Functions Used By The Active Entry Point

The active figure entry point uses:

- `stop_if_missing_packages_v5()`
- `newest_file_snapshot_v5()`
- `log_msg_v5()`
- `render_v5_panels()`
- `assemble_v5_figures()`
- `manifest_long_v5()`
- `wording_audit_rows_v5()`
- `write_tsv_v5()`
- `validate_v5_outputs()`
- `write_generation_log_v5()`

## Internal Figure Helper Responsibilities

`render_v5_panels()` uses the repo-contained panel builder and helper routines for:

- input-table loading and plot construction through the `make_Figure*()` functions in `lib/panel_builders_v6.R`;
- split support-panel construction through `render_figure2_support_split_v5()`;
- corrected landscape display through `render_corrected_figure1b_landscape_v5()`;
- workflow rendering through `render_figure1a_merged_workflow_v5()` and `lib/merged_workflow_v6.R`;
- PNG/PDF/SVG output through `save_plot_all_formats_v5()` and `save_grid_all_formats_v5()`;
- panel and figure plans through `make_v5_panel_plan()` and `make_v5_figure_plan()`.

The assembly and validation functions generate figure manifests, wording audits, output checks, and generation logs without refitting IMRS or changing released numeric source tables.

## Refactor Decision

The local v6 script folder did not contain a different newer helper implementation; it contained the same helper logic together with the launcher and workflow renderer. The active release path now includes that logic and a refactored panel builder under `scripts/active_manuscript/lib/`, configured through repository-relative paths. Released figure source tables are staged explicitly at `data/derived/figure_inputs/`, so active execution does not require local `<local_project_root>` paths or `revised_plots_v2`, `revised_plots_v3`, or `revised_plots_v4` directories.
