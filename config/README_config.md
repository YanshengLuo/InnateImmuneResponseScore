# Configuration

Copy `config_template.yml` to `config.yml` before running repository scripts. The template uses repository-relative paths; `project_root: "."` is appropriate when running from the repository root.

The Layer 2 manuscript-output runner expects released inputs under `data/derived/`, including the explicit `data/derived/figure_inputs/` bundle used by the repo-contained figure implementation. It writes generated artifacts beneath `results_release_templates/`.

Keep `execute_active_scripts: false` to run the preflight/checklist only. Set it to `true` only for controlled regeneration of manuscript figures, Priority3 enrichment outputs, and supplementary tables.

`run_internal_readiness: false` excludes optional NAR G&B readiness packaging from the default reviewer-facing run. Full computational provenance and HiPerGator/HPC operations remain separate from Layer 2.

For reconstruction from prepared raw gene-count matrices and verified
metadata, use `full_pipeline_config.yml` with
`scripts/portable_full_pipeline/run_count_to_v6_outputs.R`.
`local_d_drive_config_TEMPLATE.yml` is an example local-source profile and is
not the default configuration. The count-level runner orchestrates the
path-parameterized copies of the original Steps 02 through 09 scripts under
`scripts/portable_full_pipeline/original_ported/`; it is not an independent
statistical reimplementation.

Use `--mode canonical` to validate the five-anchor frozen-weight
reconstruction and anchor scoring only. Use `--mode all_scored` to
reconstruct the same locked-anchor model and apply its frozen weights to all
available configured manuscript datasets. After an all-scored run,
`scripts/portable_full_pipeline/run_step09_to_layer2_inputs.R` writes a
separate generated Layer 2 package below `<output_root>/layer2_generated_inputs/`.
The default reviewer configuration continues to consume released
`data/derived/` inputs unless explicitly changed to the generated package.
