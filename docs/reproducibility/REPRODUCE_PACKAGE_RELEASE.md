# Reproduce-package release

Prepared on 2026-10-06 from the latest available local `MANUSCRIPT_FINAL` package, based on `v6-reproducibility-release` commit `158c485b27fb6e7859a01df557fa42d75a2f987b`.

## Layout

The reference branch's root entry points and `config/`, `data/`, `docs/`, `scripts/`, and `publication_outputs/` directories are retained. The inherited Git history is preserved.

- `publication_outputs/figures/`: eight final composites in SVG, PNG, and PDF, copied without scientific changes.
- `publication_outputs/IMRS_Final_Figures_R_Assembled_PublicationReady.pptx`: supplied eight-slide R-assembled deck with editable locked captions.
- `data/derived/gene_program_enrichment/`: retained enrichment table and gene-mapping/background audits.
- `docs/reproducibility/manuscript_final/`: supplied figure QC and provenance records. Local package-root prefixes are normalized to repository-relative paths; these records describe the source package's earlier validation, not new statistical recomputation.
- `docs/reproducibility/reproduce_package_inventory.tsv`: included canonical Git-blob sizes and SHA-256 checksums, excluding this inventory itself and Git internals. These hashes are independent of Windows checkout line-ending conversion. Verify a committed checkout with `python scripts/reviewer/verify_release_inventory.py`.

Generated outputs remain under ignored `results_release_templates/`. Temporary files, build caches, duplicate archives, draft manuscripts, and older PowerPoint revisions are excluded. Original local package files are preserved separately.

## Scientific state and limitations

The final figures use 68 scored contrasts, grouped 40/14/11/3. The supplied corrected permutation tables contain 68 contrasts with 1,000 permutations each, 49 contrasts outside the strict 95% null interval, and 43 BH-significant two-sided results.

Some supplied supplementary outputs and upstream tables still reflect the historical 70-contrast state. In particular, the existing supplementary-table build summary, Supplementary Table S2, and several robustness rows are not fully synchronized with the final figures. The source package's central figure layer filters the two human HepG2 contrasts and recomputes selected figure summaries. Organizing this release does not silently rewrite those scientific source tables or assert complete manuscript/table synchronization.

The manual publication-ready PowerPoint is absent because the previously required exact source `Plots_V2(6).pptx` was unavailable. The supplied final R deck is included. Its Figure 4 uses the authoritative PNG fallback; the remaining composites use SVG.

## Reproduction

Restore the R environment recorded in `renv.lock`. The default `run_all_manuscript_outputs_v6.R` run is preflight-only. Enabling `execute_active_scripts: true` in a local configuration executes the figure, enrichment, and supplementary-table stages.

`RUN_FIGURES_v6.R` supports a fresh checkout by copying the committed S2 snapshot into the generated output folder when needed. It does not recompute enrichment. The full manuscript runner invokes the separate enrichment stage for recomputation.

The supplied final SVG/PDF/PNG snapshots and enrichment files are immutable reference outputs in the committed release; regeneration writes into `results_release_templates/`.

## Packaging validation on 2026-10-06

All 65 R scripts parsed successfully under installed R 4.3.3. The default reviewer preflight passed with all three required stages ready and the optional internal stage disabled. `RUN_FIGURES_v6.R` completed on the newly organized checkout with no skipped panels and no figure-generation warnings; the S2 snapshot fallback and label verification passed. R emitted locale-setting warnings during startup, which did not prevent these checks.

The 24 released figure files and final R PowerPoint matched the supplied originals byte-for-byte. The release validator checks the complete inventory, eight figure stems in three formats, eight PowerPoint slides, identical permutation source copies, and final permutation counts 68/49/43. Upstream score reconstruction and enrichment were not rerun during packaging; historical supplementary-table inconsistencies remain explicitly documented above.
