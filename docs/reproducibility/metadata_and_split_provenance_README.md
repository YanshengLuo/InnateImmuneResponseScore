# Metadata And Split Provenance

The active reviewer-facing runner does not directly rebuild all metadata or split definitions from raw public metadata. Instead, it consumes staged derived manuscript inputs that already encode curated metadata, split definitions, manuscript role assignments, and boundary-context annotations.

The reviewer-facing runner consumes derived manuscript inputs that already encode curated metadata, split definitions, manuscript analysis roles, and boundary-context annotations. Curated metadata and split-definition files are provided separately so reviewers can inspect the decisions that precede scoring, validation, and supplementary table generation.

Curated metadata and split-definition files are staged for transparency and auditability under:

- `data_release_templates/curated_metadata`
- `data_release_templates/split_designs`
- `docs/manual_curation`

These files document how treatment/control groups, tissue and timepoint labels, manuscript analysis roles, and late or context-shifted boundary categories were defined before downstream scoring and manuscript table generation.

Full metadata-to-score reconstruction belongs to the Level 2/full-pipeline path. That path is separate from the Level 1 reviewer-facing reproduction layer and may require public data retrieval, additional configuration, and HiPerGator/HPC resources.

Manual curation defines metadata, split, role, and interpretation context. It does not manually alter delivery-minus-control ΔIMRSz values.
