# Code and Data Availability Draft

## A. Pre-submission / Private-review Draft

Analysis code for the frozen, transfer-oriented IMRS scoring framework, reviewer-facing manuscript-output reproduction layer, retained-gene enrichment analysis, supplementary table generation, and reproducibility documentation will be made available at `TO_BE_ADDED_GITHUB_REPOSITORY_URL`. For peer review, a reviewer-access archive or private repository link can be provided at `TO_BE_ADDED_REVIEWER_ACCESS_LINK` if needed.

The release candidate includes staged derived inputs required by the Level 1 reviewer-facing runner, regenerated Supplementary Figure S2 outputs, regenerated Supplementary Tables S1-S5, and NAR G&B readiness materials. The default runner starts from staged derived inputs and does not run raw-data retrieval, full framework reconstruction, frozen-gene reconstruction, HiPerGator/SLURM scripts, or archive/legacy scripts.

## B. Final Publication Draft

Analysis code and reproducibility materials for the IMRS manuscript will be available at `TO_BE_ADDED_GITHUB_REPOSITORY_URL` and archived at `TO_BE_ADDED_ZENODO_DOI` or an equivalent archival DOI. Public dataset accessions are listed in the repository materials and manuscript supplementary tables. Released derived inputs include frozen IMRS gene-weight tables, metadata/provenance summaries, split-level manuscript inputs, robustness summaries, retained-gene enrichment results, and reviewer-facing regeneration outputs.

Raw sequencing files should generally be retrieved from public repositories by accession rather than redistributed, unless redistribution is explicitly allowed and practical. The repository documents three levels of reproducibility: Level 1 reviewer-facing manuscript-output reproduction from staged derived inputs, Level 2 analytic reconstruction from gene counts and curated metadata/split definitions, and Level 3 raw-data/HPC reconstruction from public accessions and alignment/counting workflows.

No new URLs or DOIs are invented in this draft. Replace placeholders only after the public GitHub repository and archival release are available.
