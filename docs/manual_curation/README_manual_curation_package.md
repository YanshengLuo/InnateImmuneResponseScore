# IMRS Manual Curation Package

Generated: 2026-05-23T17:33:37.913Z

## Purpose

This package documents manually defined and curated inputs used by the IMRS manuscript workflow. It focuses on metadata, split definitions, manuscript role assignments, publication-context mapping, boundary-context annotations, and manually placed caption/subtitle text. It is documentation only and does not rerun or alter IMRS scoring, validation, robustness, enrichment, figure, or supplementary-table outputs.

## Why Manual Curation Was Necessary

Several public transcriptomic datasets used in IMRS do not encode delivery group, control group, tissue, timepoint, or manuscript interpretation role in a uniformly machine-readable format. These fields were therefore curated into versioned tabular inputs with source evidence and rationale. Downstream scoring, validation, robustness, and figure-generation scripts consume these curated inputs without manual editing of numerical IMRS score values.

## What Was Curated

The curation package covers dataset role assignments, manuscript analysis groups, treatment/control labels, delivery-versus-control split definitions, tissue/timepoint annotations, delivery platform or modality labels, publication-context mappings, inclusion/exclusion or primary-claim handling, boundary categories for late or context-shifted settings, and manual figure-caption/subtitle text.

## Where Curated Files Are Located

Primary current sources include:

- <local_project_root>\audit\results\manuscript_dataset_role_table.tsv
- <local_project_root>\audit\results\supplement_dataset_split_provenance_v7.tsv
- <original_manuscript_output_folder>\supplementary_tables\Supplementary_Table_S3_late_context_shifted_boundary_audit.tsv
- <original_manuscript_output_folder>\v6_manual_caption_replacement_table.tsv
- <original_manuscript_output_folder>\Priority3_gene_program_enrichment\FigureS2_caption_candidate.txt
- <local_project_root>\audit\results\all_verified_metadata_folder_inventory.tsv
- <local_project_root>\audit\results\all_split_design_inventory.tsv

Verified metadata and split-design tables are located under `<local_project_root>/00_metadata/verified_metadata` and the mirrored HiPerGator source tree when present. Current manuscript-facing curated summaries are also represented in the v6 supplementary tables package.

## How Manual Curation Affects Reproducibility

Manual curation defines the metadata context consumed by automated scripts. It identifies which samples and groups form each delivery-versus-control split, how tissue/timepoint/platform labels are harmonized, and how datasets are assigned to manuscript analysis roles. These curated fields should be released as versioned tables so reviewers can inspect the decisions that precede automated scoring and plotting.

## Downstream Use

Downstream scripts consume curated tables for split-level scoring, transfer evaluation, role cleanup, boundary-context auditing, figure generation, and supplementary table packaging. The curated fields do not manually alter delivery-minus-control 螖IMRSz values, secondary AUC values, robustness metrics, or enrichment statistics.

## Boundary Context

Late or context-shifted settings are retained in the manuscript package and documented as boundary-setting evidence. They are not hidden or treated as unstructured negative evidence; they define the scope in which the frozen, transfer-oriented IMRS scoring framework should be interpreted.

## Captions and Manual Assembly

Current v6 caption handling is documented in `manual_caption_subtitle_audit.tsv`. No PowerPoint or Word caption source was found in the v6 folder during this packaging pass. Figure 1-5 and Supplementary Figure S1 caption updates remain manual placement items if the final manuscript/assembled figure source is outside v6. Supplementary Figure S2 has a caption candidate text file in the Priority3 enrichment folder.

## Files in This Package

- `manual_inputs_file_inventory.tsv`: inventory of curated/manual source files and generated tables documenting curation.
- `manual_curation_manifest.tsv`: row-level manifest of role, split, label, annotation, publication-context, boundary-context, and caption curation records.
- `dataset_role_curation_audit.tsv`: dataset-level role summary and review flags.
- `split_definition_curation_audit.tsv`: split-level delivery-versus-control definitions.
- `treatment_control_label_audit.tsv`: treatment and control label harmonization records.
- `boundary_context_curation_audit.tsv`: controlled audit for late or context-shifted settings.
- `publication_context_mapping_audit.tsv`: publication/provenance mapping summaries.
- `manual_caption_subtitle_audit.tsv`: manual caption/subtitle update records.
- `manual_curation_warnings.tsv`: missing metadata, ambiguous labels, and manual-caption follow-up warnings.

## No Analysis Changes

This package was generated from existing source tables and text files. It does not modify scoring, validation, robustness, enrichment, figure, supplementary table, or manuscript result values.

