# Final Pre-GitHub Summary

## What Is Ready

- The clean repository structure is present and previously tested for Level 1 reviewer-facing reproduction.
- `run_all_manuscript_outputs_v6.R` is present and configured to use staged Level 1 derived inputs.
- `config/config_template.yml` is present; local `config/config.yml` is ignored by `.gitignore`.
- `CITATION.cff` has been created with conservative placeholder-safe metadata.
- `LICENSE_DECISION_NEEDED.md` documents that a real license requires PI/institution/project approval.
- `renv.lock` has been created from detected repository dependencies without running analysis scripts.
- `data/derived`, `data/curated_metadata`, `data/split_designs`, `docs/reproducibility`, and `docs/manual_curation` are present.

## Manual Actions Still Needed Before Public Release

- Choose a real license and replace `LICENSE_PLACEHOLDER.txt` with `LICENSE` after approval.
- Update `CITATION.cff`: repository URL, final author list, final license, release date, and DOI/manuscript citation when available.
- Replace GitHub URL, reviewer-access, and Zenodo DOI placeholders in availability documents and generated candidate statements.
- Review `renv.lock` after `renv::restore(prompt = FALSE)` in a clean clone; local snapshot warned that transitive packages `fansi` and `plogr` were absent locally.
- Optionally remove ignored local test artifacts before packaging a final tarball; they are protected by `.gitignore`.

## Folder To Initialize As Git Repository

Initialize Git from the repository root containing this file: `<repository_root>`.

## Files Not To Commit

- `config/config.yml`
- `results_release_templates/`
- `results_release_templates/logs/`
- `local_test_logs/`
- `*.log`
- `_stdout.log` / `_stderr.log` files
- `.Rhistory`, `.RData`, `.Rproj.user/`
- generated local test TSV/TXT reports listed in `.gitignore`

## Recommended First Commit Message

`Initial IMRS reproducibility package`

## Reminder

Create the GitHub repository, make the first commit from this clean repo folder, push to GitHub, then create a release and Zenodo/archive DOI when ready. Do not run raw-data reconstruction, full framework reconstruction, HiPerGator/SLURM scripts, archive/legacy scripts, or manuscript analysis during this administrative commit preparation.
