# Final GitHub Release-Readiness Inspection After Fixes

Inspection date: 2026-05-25

## Executive Conclusion

**Status: Ready for public GitHub publication as the IMRS release candidate.**

The public candidate now contains a finalized MIT code license, an explicit data/documentation release note, completed repository citation metadata without an unissued DOI, public-facing code/data availability text using the actual GitHub repository URL, and a repository-contained Layer 2 reviewer runner. Account-specific HPC submission scripts have been excluded from the public candidate. Local-source path text in released figure-input provenance fields has been portabilized without changing statistical or scientific values.

## Fixes Applied

- Added `LICENSE` using the MIT License for repository code.
- Added `DATA_LICENSE.md` to distinguish repository code/derived documentation from original public RNA-seq source terms.
- Updated `CITATION.cff` with the repository URL and MIT license; no DOI or unknown release date was asserted.
- Updated code/data availability drafts to use `https://github.com/YanshengLuo/InnateImmuneResponseScore` and to state that no archival DOI has been assigned in this release candidate.
- Updated the optional internal readiness generator so it will not recreate obsolete availability placeholders if it is invoked later.
- Removed `scripts/hpc_hypergator/` from this public candidate and added it to `.gitignore`; original staging/project copies were not modified.
- Updated public documentation to state that account-specific HPC scripts are excluded and may be supplied later only as sanitized templates.
- Replaced historical local-root prefixes in nine `data/derived/figure_inputs/*.tsv` provenance fields with `<upstream_project_root>`.
- Updated `data/derived/figure_inputs/figure_input_manifest.tsv` and `data/derived/figure_inputs/README.md` with the public-release portability note.
- Refreshed `docs/RELEASE_CONTENTS.tsv` to omit ignored/generated/local files and excluded HPC material.
- Retained `config/config.yml` only as ignored local state with active execution disabled.

## Value-Preservation Validation

The path portabilization affected provenance text only.

- Path-bearing figure-input tables sanitized: 9.
- Path-prefix occurrences replaced: 3,625.
- Files with remaining original local-root prefix in `data/derived/figure_inputs/`: 0.
- Row and column counts before/after: unchanged for all nine files.
- Replacement validation: each edited table matched the expected content obtained by replacing only the local path prefix.

The validation record is `docs/reproducibility/figure_input_path_portabilization_audit.tsv`.

No IMRS scoring values, validation values, robustness statistics, enrichment statistics, gene lists, term identifiers, figure-data numeric fields, or supplementary-table source numeric fields were changed.

## Remaining Non-Blocking Findings

- Historical Layer 1 executables are isolated under `scripts/full_pipeline/provenance_legacy/` because they retain original local-root defaults. The Layer 1 README and inventory state that they are documentation/provenance and are not intended to run as-is from a fresh clone.
- Two non-default utility scripts with historical local-root assumptions (`scripts/utilities/project_tree.R` and `scripts/utilities/temporal_code.R`) were excluded from the public branch because they are not invoked by the Layer 2 runner.
- The prior pre-fix inspection report intentionally describes findings that have since been addressed; it is marked as superseded.
- Local `results_release_templates/` and `local_test_logs/` still contain generated execution traces from validation. They are ignored and are not recommended for the public source commit.
- `config/config.yml` is a local ignored configuration file and must not be committed.

## Privacy And Security Rescan

Public-content scanning after fixes found:

- No personal institutional or Gmail addresses in public candidate content outside ignored generated folders.
- No public distributed HPC submission folder or account-specific cluster script paths.
- No local filesystem prefix remains in released Layer 2 figure-input data.
- No repository URL, DOI, license, or release-date placeholder token remains in citation or public availability documents.
- No password, API-key, client-secret, bearer-token, or private-key credential was detected.

Matches for `token` in data and code refer to split/group token fields or variable names. Historical path defaults in Layer 1/provenance utilities are documented and not required by the default runner.

## Layer 2 Dry-Run Validation

Dry run was executed with:

- `execute_active_scripts: false`
- `run_internal_readiness: false`

Checklist result:

| Step | Classification | Status |
| --- | --- | --- |
| `manuscript_figures` | required | ready |
| `priority3_gene_program_enrichment` | required | ready |
| `supplementary_tables_s1_s5` | required | ready |
| `internal_readiness` | optional | disabled_by_config |

The checklist was written to `results_release_templates/manifests/active_manuscript_step_checklist.tsv`. No active manuscript script executed during this post-fix dry run, and checked existing manuscript-output files were unchanged.

## Files Removed Or Ignored

Removed from this public candidate:

- `scripts/hpc_hypergator/`, because it contained account-specific submission information and is not required by Layer 2.
- Local R session debris identified during inspection.
- Archived citation placeholder text.

Ignored and not recommended for commit:

- `config/config.yml`
- `results_release_templates/`
- `local_test_logs/`
- any local logs, R session files, or reintroduced `scripts/hpc_hypergator/` content

## Files Recommended For Commit

- `README.md`, `.gitignore`, `LICENSE`, `DATA_LICENSE.md`, `CITATION.cff`, and `renv.lock`.
- `config/config_template.yml` and `config/README_config.md`.
- `run_all_manuscript_outputs_v6.R` and `run_full_pipeline_from_public_data_TEMPLATE.R`.
- `scripts/active_manuscript/` and its repository-contained helpers.
- `scripts/full_pipeline/`, with historical executables grouped under `provenance_legacy/`, as documented Layer 1 provenance.
- Reviewed `scripts/utilities/` and `scripts/optional_internal/` with their stated non-default roles.
- `data/accessions/`, `data/curated_metadata/`, `data/split_designs/`, and the sanitized `data/derived/` Layer 2 release inputs.
- Public-facing documentation beneath `docs/`.

## Recommended Commit Message

```text
Add IMRS two-layer reproducibility release candidate
```
