# Final GitHub Release-Readiness Inspection

Inspection date: 2026-05-25

**Status note:** This document records the pre-fix inspection state. See `FINAL_GITHUB_RELEASE_READINESS_INSPECTION_AFTER_FIXES.md` for the final public-release assessment after hygiene fixes were applied.

## A. Executive Conclusion

**Status: Ready with required minor release fixes before public GitHub publication.**

The Layer 2 reviewer-facing reproduction path is functional and repository-contained. A dry run completed with all three required steps ready, followed by a controlled execution that regenerated manuscript figures, Priority3 gene-program enrichment outputs, and Supplementary Tables S1-S5 beneath `results_release_templates/`. Optional internal readiness packaging remained disabled.

The repository should not be published unchanged because license/citation placeholders remain, the HPC provenance scripts contain account-specific cluster paths and a personal notification address, and nine released figure-input tables retain historical local-source path strings. These findings do not affect scientific result values or Layer 2 execution, but they should be addressed or deliberately documented before a public push.

## Repository Inventory

| Folder | Intended purpose | Release assessment |
| --- | --- | --- |
| `config/` | Repository-relative Layer 2 configuration template and local test config | Commit template/docs; exclude ignored `config/config.yml`. |
| `data/` | Public accessions, curated metadata/split definitions, and released derived Layer 2 inputs | Appropriate in principle; sanitize historical local-path fields in nine figure-input tables before public release. |
| `docs/` | Reproducibility, availability, manual-curation, and audit documentation | Appropriate; update administrative placeholders. |
| `scripts/active_manuscript/` | Layer 2 manuscript-output implementation and repo-contained helpers | Required public release content; validated. |
| `scripts/full_pipeline/` | Layer 1 count-to-score computational provenance | Useful provenance content; not default runnable and documented as needing future path portability work. |
| `scripts/hpc_hypergator/` | Level 3 raw-data/HPC provenance | Optional public content only after account-path and email sanitization, or exclude from public source release. |
| `scripts/optional_internal/` | Internal readiness/legacy review utilities | Optional/internal; not invoked by default Layer 2 runner. |
| `scripts/utilities/` | Additional utilities and historical release notes | Needs manual triage; two utilities retain local absolute defaults and are not required by Layer 2. |
| `results_release_templates/` | Locally regenerated outputs and logs | Generated local material; currently ignored and should not be committed by default. |
| `local_test_logs/` | Local execution validation logs | Local-only and ignored. |

## B. Release-Blocking Issues

1. **No finalized license is present.** No `LICENSE` file or license-decision notice is currently included. Select and add an approved license before public publication.
2. **Citation and availability metadata are unfinished.** `CITATION.cff` contains repository URL, license, and release-date placeholders. Code/data availability drafts contain GitHub, reviewer-access, and archival DOI placeholders.
3. **Level 3 HPC provenance includes private/machine-specific identifiers.** `scripts/hpc_hypergator/` included institutional email references and account-specific cluster storage/environment paths. Replace with configurable placeholders or exclude these scripts from the public GitHub version.
4. **Released figure-input tables include historical local filesystem paths.** Nine files beneath `data/derived/figure_inputs/` contain local-source path strings. Replace only the provenance path fields with portable identifiers or clearly documented relative references before committing; numerical and scientific fields must remain unchanged.

## C. Minor Hygiene Issues

- The inspected folder is not currently a Git worktree, so `git status --short` and `git grep` cannot classify tracked changes here. Perform the final status review after copying into the existing working clone or initializing the intended repository.
- `config/config.yml` is present as an ignored local test configuration. It has been restored to `execute_active_scripts: false` and `run_internal_readiness: false`; do not commit it.
- The `Rscript` executable is installed but was not on `PATH` in the inspection shell. Reviewer instructions should state that R must be available on `PATH` or invoked by its installed executable path.
- `scripts/full_pipeline/` intentionally preserves historical local-root defaults for provenance and is correctly documented as not the default reviewer run. It should not be advertised as portable execution until configuration cleanup is completed.
- `scripts/utilities/project_tree.R` and `scripts/utilities/temporal_code.R` contain historical local-root assumptions and should be triaged as provenance/internal utilities or made portable before public presentation.
- A prior `results_release_templates/NAR_GB_readiness/` output folder remains locally from earlier work. Since `results_release_templates/` is ignored and NAR readiness is disabled by default, it should not enter the public source commit.
- `REVIEWER_REPRODUCTION_STATUS.md` is not present at repository root. It is optional, but a short public status document would help reviewers.

Safe hygiene applied during inspection:

- Removed local R session debris: root `.RData` and `scripts/full_pipeline/.Rhistory`.
- Reset ignored local `config/config.yml` from execution mode to checklist mode.

## D. Files To Commit

Commit after blockers are addressed:

- `README.md`, `.gitignore`, finalized `LICENSE`, finalized `CITATION.cff`, and `renv.lock`.
- `config/config_template.yml` and `config/README_config.md`.
- `run_all_manuscript_outputs_v6.R` and `run_full_pipeline_from_public_data_TEMPLATE.R`.
- `scripts/active_manuscript/`, including `lib/active_config.R`, `lib/figure_helpers_v6.R`, `lib/panel_builders_v6.R`, and `lib/merged_workflow_v6.R`.
- `scripts/full_pipeline/README_full_pipeline.md` and `scripts/full_pipeline/full_pipeline_script_inventory.tsv`; include additional Layer 1 scripts as provenance only with their documented limitations.
- `data/accessions/`, `data/curated_metadata/`, `data/split_designs/`, and sanitized `data/derived/` release inputs.
- Public-facing files under `docs/reproducibility/` and `docs/manual_curation/` after placeholder/path review.

Optional commit after sanitization:

- `scripts/hpc_hypergator/` as Level 3 provenance templates.
- `scripts/optional_internal/` only if the internal-status purpose is explicitly useful to reviewers.
- A curated set of expected-output references, if the release policy chooses to provide them outside ignored generated-output folders.

## E. Files To Remove Or Ignore

- `config/config.yml` (local-only; covered by `.gitignore`).
- `local_test_logs/` and all execution logs (covered by `.gitignore`).
- `results_release_templates/` by default; provide regenerated outputs in an archival/reviewer artifact instead if desired.
- Local R session state such as `.RData` and `.Rhistory` (removed during inspection and covered by `.gitignore`).
- Any future raw FASTQ/BAM files, caches, scratch/render folders, or private review exports.

## F. Path, Privacy, And Security Findings

- **Active Layer 2 code:** no dependency on external/local `v5_helpers.R`, archive scripts, HPC scripts, full-pipeline scripts, or legacy `revised_plots_v2/v3/v4` locations was found.
- **Historical/provenance code:** Layer 1 and utility scripts include historical local-root assumptions; these are not called by the default runner and are documented as nonportable provenance.
- **HPC privacy:** personal email and account-specific cluster paths remain in the Level 3 HPC script folder and require sanitization or exclusion.
- **Released data paths:** nine Layer 2 figure-input tables include historical local path text in provenance/source columns. This is portability leakage, not evidence of altered numeric results.
- **Credentials:** no passwords, API keys, bearer tokens, authorization secrets, or private keys were detected. Matches for the word `token` were ordinary split/group token field names.
- **Data scope:** no FASTQ, BAM, CRAM, or SAM files and no file over 100 MB were found. Public accession documentation exists with 16 dataset rows. Manual review should still confirm that all released metadata are public/de-identified and appropriate to redistribute.

## G. Layer 2 Validation Result

Dry-run validation:

- Command was executed with the installed R 4.3.3 `Rscript` executable because `Rscript` was not available on this shell's `PATH`.
- `execute_active_scripts: false`.
- Required steps ready: `manuscript_figures`, `priority3_gene_program_enrichment`, `supplementary_tables_s1_s5`.
- Optional step disabled: `internal_readiness`.
- No manuscript figure, enrichment, or supplementary table output was regenerated in dry-run mode.
- Checklist written to `results_release_templates/manifests/active_manuscript_step_checklist.tsv`.

Controlled Layer 2 validation:

- `manuscript_figures`: executed successfully.
- `priority3_gene_program_enrichment`: executed successfully.
- `supplementary_tables_s1_s5`: executed successfully.
- `internal_readiness`: skipped as `disabled_by_config`.
- Generated outputs were written beneath `results_release_templates/` only.
- No files in the external staging repository or original manuscript-output folder changed during this controlled run.
- The local config was restored to `execute_active_scripts: false` after validation.

## H. Documentation Assessment

The documentation correctly distinguishes the two-layer structure:

- `README.md` identifies Layer 2 as the default reviewer-facing quick run and states that it regenerates manuscript figures, supplementary tables, and Priority3 enrichment outputs from released derived inputs and scoring/evaluation tables.
- `README.md`, `config/README_config.md`, and the Layer 2 validation summary state that NAR G&B readiness is optional/internal and disabled by default.
- `scripts/full_pipeline/README_full_pipeline.md` documents Layer 1 from public gene-count matrices and verified metadata to frozen weights, sample-level scores, and Step09 evaluation tables without presenting it as the default reviewer run.
- `docs/reproducibility/frozen_gene_reconstruction_README.md` states that frozen-gene reconstruction is separate from the reviewer-facing run.
- `docs/reproducibility/metadata_and_split_provenance_README.md` documents the role of curated metadata and split definitions.

Priority3 gene-program enrichment is included as a required Layer 2 step and is explicitly checked by the runner.

## I. Recommended Final Commit Message

After resolving public blockers:

```text
Add IMRS two-layer reproducibility release candidate
```

## Recommended Next Actions

1. Select and add an approved open-source `LICENSE`.
2. Finalize `CITATION.cff` and replace availability-statement URL/DOI placeholders.
3. Sanitize or exclude `scripts/hpc_hypergator/` before a public push.
4. Sanitize local-source path fields in the nine released figure-input TSVs without changing scientific values, followed by a Layer 2 equivalence verification.
5. Transfer this cleaned content into a Git worktree and run `git status --short` before committing.
