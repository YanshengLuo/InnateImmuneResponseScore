# IMRS v1.1 portability fixes

This copy was audited specifically for the "extract the repository and run it" workflow.

## Recommended figure workflow

1. Extract the ZIP without moving files out of the repository structure.
2. Open `RUN_FIGURES_v6.R` in RStudio.
3. Use **Source** or **Run**. No manual `setwd()` is required.
4. Generated figure outputs are written under `results_release_templates/figures/` according to `config/config_template.yml` unless a repository `config/config.yml` or `IMRS_ACTIVE_CONFIG` override is supplied.

The original active entry point remains supported:

` scripts/active_manuscript/00_generate_manuscript_figures_v6.R `

It now also discovers the repository root automatically when run from RStudio or Rscript.

## Portability fixes applied

- Removed the duplicated `scripts/active_manuscript/scripts/active_manuscript/...` failure mode.
- Active manuscript entry scripts now detect their own location under Rscript, source/sys.source, and RStudio Run/Source.
- Repository-root discovery now walks upward until repository markers are found instead of assuming a fixed working directory depth.
- Root `run_all_manuscript_outputs_v6.R` uses the same working-directory-independent discovery logic.
- Supported portable-pipeline and bridge entry scripts no longer depend exclusively on the Rscript `--file` argument.
- `scripts/reviewer/check_clean_count_inputs.R` can now also be invoked from RStudio.
- Fixed the optional NAR G&B readiness output path to use the configured `internal_readiness_dir`.
- Added an RStudio project file for convenient root-based project opening.

## Intentionally unchanged

- `config/local_d_drive_config_TEMPLATE.yml` contains example `D:/IMRS_Project/...` paths by design. It is a local configuration template and is not used by the default figure regeneration route.
- Legacy/provenance launchers remain archival and are not used by `RUN_FIGURES_v6.R` or the active manuscript figure route.
- Package installation is not performed automatically. Versions are recorded in `renv.lock`; use `renv::restore()` when the required R packages are not already installed.
- The root `run_all_manuscript_outputs_v6.R` retains its controlled-regeneration setting (`execute_active_scripts: false` in the template). Use `RUN_FIGURES_v6.R` for direct figure-only regeneration.
