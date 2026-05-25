# Frozen-Gene Reconstruction Provenance

The reviewer-facing manuscript reproduction path starts from:

- `data/derived/frozen_gene_weights.tsv`
- `data/derived/gene_power.tsv`

This is intentional because manuscript-output reproduction should not require rerunning the full anchor model construction. Those frozen derived files are the stable inputs used to regenerate manuscript-facing outputs such as Supplementary Figure S2 and the related supplementary tables.

The reviewer-facing manuscript reproduction path starts from frozen derived inputs to regenerate figures and supplementary tables. The full framework reconstruction path documents how frozen_gene_weights.tsv can be regenerated from upstream clean-gene and locked-anchor inputs. These are intentionally separated to keep the reviewer reproduction path stable while preserving provenance for the frozen IMRS coefficient table.

## Conceptual Workflow

The intended full reconstruction workflow is:

1. clean-gene / anchor-eligible inputs
2. locked-anchor support and filtering
3. frozen retained genes
4. frozen gene weights / coefficients
5. frozen IMRS scoring

The released `frozen_gene_weights.tsv` is the canonical frozen coefficient table used by manuscript figures, retained-gene enrichment, scoring, and transfer evaluation. Manual editing of frozen weights is not part of the manuscript workflow.

## Staged Provenance Files

Upstream reconstruction inputs are staged under:

- `data/derived/frozen_gene_reconstruction_inputs/`

Upstream reconstruction scripts are staged under:

- `scripts/full_pipeline/frozen_gene_reconstruction/`

These files document the reconstruction path. They are not called by `run_all_manuscript_outputs_v6.R` and are not part of the default reviewer-facing active manuscript run.
